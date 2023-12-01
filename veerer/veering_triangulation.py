r"""
Veering triangulations of surfaces.

An edge permutation is an involution (possibly with fixed points corresponding
to folded edges). It is said to be in canonical form if it starts with cycle
representatives.
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2018 Mark Bell
#                     2018-2023 Vincent Delecroix
#                     2018 Saul Schleimer
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************

import collections
import itertools
import numbers
from random import choice, shuffle
from array import array
import ppl

from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE, rich_to_bool

from .constants import *
from .permutation import *
from .misc import det2
from .triangulation import Triangulation
from .polyhedron import LinearExpressions, ConstraintSystem
from .polyhedron.linear_algebra import linear_form_project, linear_form_normalize


class VeeringTriangulation(Triangulation):
    r"""
    Veering triangulation.

    A *veering triangulation* is a triangulation of a surface together with
    a colouring of the edges in red or blue so that there is no monochromatic
    face and no monochromatic vertex.

    EXAMPLES::

        sage: from veerer import *  # random output due to deprecation warnings from realalg

    Built from an explicit triangulation (in cycle or list form) and a list of colours::

        sage: VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
        VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB")

    From a stratum::

        sage: from surface_dynamics import *  # optional - surface_dynamics

        sage: VeeringTriangulation.from_stratum(AbelianStratum(2))  # optional - surface_dynamics
        VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")

        sage: vt = VeeringTriangulation.from_stratum(QuadraticStratum({1:4}))  # optional - surface_dynamics
        sage: vt.stratum()                                                     # optional - surface_dynamics
        Q_2(1^4)

    From a flipper pseudo-Anosov map::

        sage: import flipper                # optional - flipper

        sage: T = flipper.load('S_2_1')     # optional - flipper
        sage: h = T.mapping_class('abcD')   # optional - flipper
        sage: h.is_pseudo_anosov()          # optional - flipper
        True
        sage: V = VeeringTriangulation.from_pseudo_anosov(h)  # optional - flipper
        sage: V # optional - flipper
        VeeringTriangulation("(0,~3,~1)(1,2,14)...(12,~14,~10)(~9,~4,~2)", "RBRBRRBRBBBBRBR")

    A triangulation with some purple::

        sage: VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBPBRBPRB")
        VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBPBRBPRB")
    """
    __slots__ = ['_colouring']

    def __init__(self, triangulation,  colouring, mutable=False, check=True):
        Triangulation.__init__(self, triangulation, mutable=mutable, check=False)
        n = self._n  # number of half edges (get initialized by the triangulation)

        if isinstance(colouring, str):
            colouring = [colour_from_char(c) for c in colouring]

        # set _colouring: half edge index --> {RED, BLUE}
        if len(colouring) == self.num_edges():
            colouring = [colouring[self._norm(i)] for i in range(n)]
        elif len(colouring) == n:
            colouring = colouring
        else:
            raise ValueError("'colouring' argument of invalid length")

        self._colouring = array('i', colouring)

        if check:
            self._check(ValueError)

    def _check(self, error=RuntimeError):
        """
        EXAMPLES::

            sage: from veerer import VeeringTriangulation
            sage: V = VeeringTriangulation("(0,1,2)", "GBR")
            Traceback (most recent call last):
            ...
            ValueError: invalid triangle (0, 1, 2) with colours (green, blue, red)
        """
        Triangulation._check(self, error)
        n = self.num_half_edges()
        ep = self._ep
        ev = self._vp
        cols = self._colouring
        if not isinstance(cols, array) or \
           len(cols) != n or \
           any(col not in COLOURS for col in cols) or \
           any(cols[e] != cols[ep[e]] for e in range(self._n)):
            raise error('bad colouring attribute')

        # faces must be of one of the following type (up to cyclic ordering)
        # non-dgenerate: BBR (BLUE), RRB (RED)
        # 1-degenerate: PBR (PURPLE), GRB (GREEN)
        # 2-degenerate: BPG (BLUE|PURPLE|GREEN), RGP (RED|GREEN|PURPLE)
        for a in range(n):
            col, a, b, c = self.triangle(a)
            good = False
            if col == BLUE:
                good = cols[a] == BLUE and cols[b] == BLUE and cols[c] == RED
            elif col == RED:
                good = cols[a] == RED and cols[b] == RED and cols[c] == BLUE
            elif col == PURPLE:
                good = cols[a] == PURPLE and cols[b] == BLUE and cols[c] == RED
            elif col == GREEN:
                good = cols[a] == GREEN and cols[b] == RED and cols[c] == BLUE
            elif col == BLUE | PURPLE | GREEN:
                good = cols[a] == BLUE and cols[b] == GREEN and cols[c] == PURPLE
            elif col == RED | PURPLE | GREEN:
                good = cols[a] == RED and cols[b] == PURPLE and cols[c] == GREEN
            if not good:
                raise error('invalid triangle ({}, {}, {}) with colours ({}, {}, {})'.format(a, b, c,
                    colour_to_string(cols[a]), colour_to_string(cols[b]), colour_to_string(cols[c])))

        # no monochromatic vertex
        for v in self.vertices():
            col = cols[v[0]]
            i = 1
            while i < len(v) and cols[v[i]] == col:
                i += 1
            if i == len(v):
                raise error('monochromatic vertex {} of colour {}'.format(v, colour_to_string(cols[v[0]])))

    def __getstate__(self):
        r"""
        TESTS::

            sage: from veerer import VeeringTriangulation
            sage: t = VeeringTriangulation("(0,1,2)", "BBR")
            sage: dumps(t)  # indirect doctest
            b'x\x9ck`J.KM-J-\xd2\x03Q\x99y\xe9\xf1%E\x99\x89y\xe9\xa59\x89%\x99\xf9y\\a\x10\xd1\x10\x14\xc1B\x06\xcd\xc6B\xc6\xd8B&\roFo&o\x06 \x04\xd1 \xc8\xd8\x99\x9a\xa4\x07\x005\xe2\x1bc'
        """
        a = list(self._fp)
        a.extend(self._ep)
        a.extend(self._colouring)
        a.append(self._mutable)
        return a

    def __setstate__(self, arg):
        r"""
        TESTS::

            sage: from veerer import VeeringTriangulation
            sage: t0 = VeeringTriangulation("(0,1,2)", "BBR", mutable=False)
            sage: t1 = VeeringTriangulation("(0,1,2)", "BBR", mutable=True)
            sage: s0 = loads(dumps(t0))  # indirect doctest
            sage: assert s0 == t0
            sage: s0._mutable
            False
            sage: s0._check()

            sage: s1 = loads(dumps(t1))  # indirect doctest
            sage: assert s1 == t1
            sage: s1._mutable
            True
            sage: s1._check()
        """
        n = (len(arg) - 1) // 3
        self._n = n
        self._fp = array('i', arg[:n])
        self._ep = array('i', arg[n:2*n])
        self._colouring = array('i', arg[2*n:3*n])
        self._mutable = arg[-1]
        self._vp = array('i', [-1] * n)
        for i in range(n):
            self._vp[self._fp[self._ep[i]]] = i

    def base_ring(self):
        from sage.rings.integer_ring import ZZ
        return ZZ

    def __hash__(self):
        r"""
        TESTS::

            sage: from veerer import *

            sage: triangulations = []
            sage: for cols in ["RRB", "RBR", "BRR", "BBR", "BRB", "RBB"]:
            ....:     t = VeeringTriangulation("(0,1,2)", cols)
            ....:     triangulations.append(t)

            sage: for i in range(len(triangulations)):
            ....:     for j in range(len(triangulations)):
            ....:         assert (triangulations[i] == triangulations[j]) == (i == j), (i, j)
            ....:         assert (triangulations[i] != triangulations[j]) == (i != j), (i, j)

            sage: hashes1 = {}
            sage: hashes2 = {}
            sage: for t in triangulations:
            ....:     h1 = hash(t) % (2 ** 16)
            ....:     h2 = (hash(t) >> 16) % (2 ** 16)
            ....:     if h1 in hashes1:
            ....:         print('collision 1: {} {}'.format(hashes1[h1], t))
            ....:     else:
            ....:         hashes1[h1] = t
            ....:     if h2 in hashes2:
            ....:         print('collision 2: {} {}'.format(hashes2[h2], t))
            ....:     else:
            ....:         hashes2[h2] = t
            sage: assert len(hashes1) == len(hashes2) == len(triangulations), (len(hashes1), len(hashes2), len(triangulations))

            sage: t = VeeringTriangulation("(0,1,2)", cols, mutable=True)
            sage: hash(t)
            Traceback (most recent call last):
            ...
            ValueError: mutable veering triangulation not hashable
        """
        if self._mutable:
            raise ValueError('mutable veering triangulation not hashable')

        x = 140737488617563
        x = ((x ^ hash(self._vp.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._ep.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._fp.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._colouring.tobytes())) * 2147483693) + 82520 + self._n + self._n

        return x

    def right_wedges(self, slope=VERTICAL):
        r"""
        Return the (vertical or horizontal) right sides of the wedges in this veering triangulation.

        A wedge is a pair of consecutive half edges at a vertex that jumps over a
        (vertical or horizontal) separatrix.

        INPUT:

        - ``slope`` -- either ``VERTICAL`` (default) or ``HORIZONTAL``

        EXAMPLES::

            sage: from veerer import VeeringTriangulation, VERTICAL, HORIZONTAL
            sage: vt = VeeringTriangulation("(0,1,2)(3,4,5)(6,7,8)(~8,~0,~7)(~6,~1,~5)(~4,~2,~3)", "RRBRRBRRB")
            sage: vt.right_wedges()
            [1, 4, 7, 10, 16, 13]
            sage: vt.right_wedges(HORIZONTAL)
            [2, 5, 8, 9, 12, 15]
        """
        if slope == VERTICAL:
            right_colour = RED
            left_colour = BLUE
        elif slope == HORIZONTAL:
            right_colour = BLUE
            left_colour = RED
        else:
            raise ValueError('invalid slope argument')
        wedges = []
        for face in self.faces():
            for i in range(3):
                if self.edge_colour(face[i]) == right_colour and self.edge_colour(face[(i + 1) % 3]) == left_colour:
                    break
            if self.edge_colour(face[i]) != right_colour or self.edge_colour(face[(i + 1) % 3]) != left_colour:
                raise ValueError('invalid colouring')
            wedges.append(face[i])
        return wedges

    def as_linear_family(self, mutable=False):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: vt.as_linear_family()
            VeeringTriangulationLinearFamily("(0,1,2)(~2,~0,~1)", "RRB", [(1, 0, -1), (0, 1, 1)])
        """
        from .linear_family import VeeringTriangulationLinearFamily
        P = self.train_track_linear_space()
        return VeeringTriangulationLinearFamily(self, self.train_track_linear_space().lines(), mutable=mutable)

    def triangle(self, a):
        r"""
        Return a quadruple ``(colour, e0, e1, e2)`` in canonical form for the triangle with half edge ``a``.

        The canonical form concerns the order of the colour as read counter-clockwise along the
        triangle:

        - ``RED`` triangle: ``(RED, RED, BLUE)``
        - ``BLUE`` triangle: ``(BLUE, BLUE, RED)``
        - ``PURPLE`` triangle: ``(PURPLE, BLUE, RED)``
        - ``GREEN`` triangle: ``(GREEN, RED, BLUE)``
        - ``RED|PURPLE|GREEN`` triangle: ``(RED, PURPLE, GREEN)``
        - ``BLUE|GREEN|PURPLE`` triangle: ``(BLUE, GREEN, PURPLE)``

        EXAMPLES::

            sage: from veerer import VeeringTriangulation, RED, BLUE, PURPLE, GREEN

            sage: V = VeeringTriangulation("(0,1,2)", "RRB")
            sage: assert V.triangle(0) == V.triangle(1) == V.triangle(2) == (RED, 0, 1, 2)
            sage: V = VeeringTriangulation("(0,1,2)", "BBR")
            sage: assert V.triangle(0) == V.triangle(1) == V.triangle(2) == (BLUE, 0, 1, 2)
            sage: V = VeeringTriangulation("(0,1,2)", "PBR")
            sage: assert V.triangle(0) == V.triangle(1) == V.triangle(2) == (PURPLE, 0, 1, 2)
            sage: V = VeeringTriangulation("(0,1,2)", "GRB")
            sage: assert V.triangle(0) == V.triangle(1) == V.triangle(2) == (GREEN, 0, 1, 2)
            sage: V = VeeringTriangulation("(0,1,2)", "RPG")
            sage: assert V.triangle(0) == V.triangle(1) == V.triangle(2) == (RED|PURPLE|GREEN, 0, 1, 2)
            sage: V = VeeringTriangulation("(0,1,2)", "BGP")
            sage: assert V.triangle(0) == V.triangle(1) == V.triangle(2) == (BLUE|GREEN|PURPLE, 0, 1, 2)
        """
        # non-dgenerate: BBR (BLUE), RRB (RED)
        # 1-degenerate: PBR (PURPLE), GRB (GREEN)
        # 2-degenerate: RGP (RED|GREEN|PURPLE), BPG (BLUE|GREEN|PURPLE)
        b = self._fp[a]
        c = self._fp[b]
        cols = self._colouring
        degenerate = PURPLE | GREEN
        standard = RED | BLUE
        ndegenerate = bool(cols[a] & degenerate) + bool(cols[b] & degenerate) + bool(cols[c] & degenerate)
        if ndegenerate == 0:
            if cols[a] == cols[b]:
                return (cols[a], a, b, c)
            elif cols[b] == cols[c]:
                return (cols[b], b, c, a)
            elif cols[c] == cols[a]:
                return (cols[c], c, a, b)
        elif ndegenerate == 1:
            if cols[a] & degenerate:
                return (cols[a], a, b, c)
            elif cols[b] & degenerate:
                return (cols[b], b, c, a)
            elif cols[c] & degenerate:
                return (cols[c], c, a, b)
        elif ndegenerate == 2:
            if cols[a] & standard:
                return (cols[a] | degenerate, a, b, c)
            elif cols[b] & standard:
                return (cols[b] | degenerate, b, c, a)
            elif cols[c] & standard:
                return (cols[c] | degenerate, c, a, b)
        return (None, None, None, None)

    @classmethod
    def from_pseudo_anosov(cls, h, mutable=False, check=True):
        r"""
        Construct the coloured triangulation of a pseudo-Anosov homeomorphism.

        EXAMPLES::

            sage: from flipper import *  # optional - flipper
            sage: from veerer import *

            sage: T = flipper.load('SB_4')  # optional - flipper
            sage: h = T.mapping_class('s_0S_1s_2S_3s_1S_2')  # optional - flipper
            sage: h.is_pseudo_anosov()  # optional - flipper
            True
            sage: VeeringTriangulation.from_pseudo_anosov(h) # optional - flipper
            VeeringTriangulation("(0,4,~1)(1,5,3)(2,~0,~3)(~5,~4,~2)", "RBBBBR")
        """
        FS = h.flat_structure()
        n = FS.triangulation.zeta

        X = {i.label: e.x for i,e in FS.edge_vectors.items()}
        Y = {i.label: e.y for i,e in FS.edge_vectors.items()}

        triangles = [[x.label for x in t] for t in FS.triangulation]
        colours = [RED if X[e]*Y[e] > 0 else BLUE for e in range(n)]
        return VeeringTriangulation(triangles, colours, mutable=mutable, check=check)

    @classmethod
    def from_square_tiled(cls, s, col=RED, mutable=False, check=True):
        r"""
        Build a veering triangulation from a square-tiled surface (or origami).

        INPUT:

        - ``s`` - a square-tiled surface

        - ``col`` - either ``RED`` or ``BLUE``

        EXAMPLES::

            sage: from surface_dynamics import Origami            # optional - surface_dynamics
            sage: from veerer import VeeringTriangulation

            sage: o = Origami('(1,2)', '(1,3)')                   # optional - surface_dynamics
            sage: T = VeeringTriangulation.from_square_tiled(o)   # optional - surface_dynamics
            sage: T                                               # optional - surface_dynamics
            VeeringTriangulation("(0,1,2)(3,4,5)(6,7,8)(~8,~0,~7)(~6,~1,~5)(~4,~2,~3)", "RRBRRBRRB")
            sage: o.stratum()                                     # optional - surface_dynamics
            H_2(2)
            sage: T.stratum()                                     # optional - surface_dynamics
            H_2(2)

        A one cylinder example in the odd component of H(4)::

            sage: o = Origami('(1,2,3,4,5)', '(1,4,3,5,2)')       # optional - surface_dynamics
            sage: T = VeeringTriangulation.from_square_tiled(o)   # optional - surface_dynamics
            sage: o.stratum()                                     # optional - surface_dynamics
            H_3(4)
            sage: T.stratum()                                     # optional - surface_dynamics
            H_3(4)
        """
        from .features import surface_dynamics_feature
        surface_dynamics_feature.require()

        from surface_dynamics.flat_surfaces.origamis.origami_dense import Origami_dense_pyx

        if col not in [BLUE, RED]:
            raise ValueError("col must be BLUE or RED")

        if not isinstance(s, Origami_dense_pyx):
            raise ValueError("input must be an origami")

        faces = []
        n = s.nb_squares()  # so 3n edges in the veering triangulation
                            # (6n half edges)
        r = s.r_tuple()
        u = s.u_tuple()
        ep = list(range(6*n-1, -1,- 1))
        fp = [None] * (6*n)
        N = 6*n - 1
        for i in range(0, n):
            ii = 3*i
            fp[ii] = ii+1
            fp[ii+1] = ii+2
            fp[ii+2] = ii

            j = N - (3*r[i] + 2)
            k = N - (3*u[i])
            l = N - (3*i+1)
            fp[j] = k
            fp[k] = l
            fp[l] = j

        colouring = [None] * (3*n)
        colouring[::3] = [RED]*n
        colouring[2::3] = [BLUE]*n
        colouring[1::3] = [col]*n
        colouring.extend(colouring[::-1])

        return cls.from_face_edge_perms(array('i', colouring),
                                        array('i', fp),
                                        array('i', ep), mutable=mutable, check=check)

    @classmethod
    def from_stratum(cls, c, folded_edges=False, mutable=False, check=True):
        r"""
        Return a Veering triangulation from either a stratum, a component
        of stratum or a cylinder diagram.

        EXAMPLES::

            sage: from surface_dynamics import *   # optional - surface_dynamics
            sage: from veerer import *

            sage: T = VeeringTriangulation.from_stratum(AbelianStratum(2))  # optional - surface_dynamics
            sage: T                                                         # optional - surface_dynamics
            VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: T.stratum()                                               # optional - surface_dynamics
            H_2(2)

            sage: Q = QuadraticStratum(9,-1)                                # optional - surface_dynamics

            sage: CTreg = VeeringTriangulation.from_stratum(Q.regular_component()) # optional - surface_dynamics
            sage: CTreg.stratum()  # optional - surface_dynamics
            Q_3(9, -1)

            sage: CTirr = VeeringTriangulation.from_stratum(Q.irregular_component())  # optional - surface_dynamics
            sage: CTirr.stratum()  # optional - surface_dynamics
            Q_3(9, -1)

        Some examples built from cylinder diagram::

            sage: c = QuadraticCylinderDiagram('(0)-(0)')    # optional - surface_dynamics
            sage: VeeringTriangulation.from_stratum(c)       # optional - surface_dynamics
            VeeringTriangulation("(0,2,~1)(1,~0,~2)", "RBB")

            sage: c = QuadraticCylinderDiagram('(0,0)-(1,1,2,2,3,3)')  # optional - surface_dynamics
            sage: VeeringTriangulation.from_stratum(c)                 # optional - surface_dynamics
            VeeringTriangulation("(0,10,~9)(1,~8,9)(2,~6,7)(3,~4,5)(4,~3,~11)(6,~2,~5)(8,~1,~7)(11,~10,~0)", "RRRRBBBBBBBB")

            sage: c = CylinderDiagram('(0,6,4,5)-(3,6,5) (1,3,2)-(0,1,4,2)')  # optional - surface_dynamics
            sage: CT = VeeringTriangulation.from_stratum(c)                   # optional - surface_dynamics
            sage: CT.stratum()                                                # optional - surface_dynamics
            H_4(6)
        """
        from .features import surface_dynamics_feature
        surface_dynamics_feature.require()

        # TODO: for now there is no account of possible folded edges
        # in the cylinder diagram. This has to be changed in
        # surface_dynamics...
        if folded_edges:
            raise NotImplementedError

        from surface_dynamics.flat_surfaces.strata import Stratum, StratumComponent
        from surface_dynamics.flat_surfaces.separatrix_diagram import \
                CylinderDiagram, QuadraticCylinderDiagram

        # make sure c is converted into a cylinder diagram
        if isinstance(c, Stratum):
            if not c.is_connected():
                print('Warning: the stratum %s is not connected' % c)
            c = c.one_component().one_cylinder_diagram()
        elif isinstance(c, StratumComponent):
            c = c.one_cylinder_diagram()
        elif not (isinstance(c, CylinderDiagram) or isinstance(c, QuadraticCylinderDiagram)):
            raise ValueError("c must either be a stratum or a component of a stratum"
                     " or a cylinder diagram")

        # translate the cylinder diagram c to a triangulation
        seen = [False] * c.nseps()
        triangles = []

        m = c.nseps()  # current counter
        nfolded = 0

        for bot,top in c.cylinders():
            # len(bot) + len(top)
            l = m + len(top) - 1
            oldi = None
            for i in bot:
                assert isinstance(i,int)
                if oldi == i and folded_edges:
                    nfolded += 1
                    continue
                if seen[i]:
                    i = ~i
                else:
                    seen[i] = True
                triangles.append((i,l+1,~l))
                l += 1
                oldi = i
            l = m + len(top) - 1
            for i in top:
                assert isinstance(i,int)
                if oldi == i and folded_edges:
                    nfolded += 1
                    continue
                if seen[i]:
                    i = ~i
                else:
                    seen[i] = True
                triangles.append((i,~(l-1),l))
                l -= 1
                oldi = i
            # last one got wrong
            i, j, k = triangles.pop()
            triangles.append((i,~(k+len(bot)+len(top)-1),k))

            m += len(bot) + len(top)

        nseps = c.nseps()

        colours = [RED] * nseps + [BLUE] * (2*nseps)
        return VeeringTriangulation(triangles, colours, mutable=mutable, check=check)

    @classmethod
    def from_face_edge_perms(self, colouring, fp, ep, vp=None, mutable=False, check=True):
        T = VeeringTriangulation.__new__(VeeringTriangulation)
        T._n = len(fp)
        T._fp = fp
        T._ep = ep

        if vp is None:
            n = T._n
            fp = T._fp
            ep = T._ep
            vp = array('i', [-1] * n)
            for i in range(n):
                vp[fp[ep[i]]] = i
        T._vp = vp
        T._colouring = colouring
        T._mutable = mutable

        if check:
            T._check(ValueError)

        return T

    @classmethod
    def from_string(cls, s, mutable=False, check=True):
        r"""
        Deserialization from the string ``s``.

        Such string is typically obtained from :meth:`to_string` or
        :meth:`iso_sig`. Not by calling ``str(triangulation)``.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation.from_string('RRB_120_012')
            sage: T.to_string()
            'RRB_120_012'

            sage: T = VeeringTriangulation.from_string('RRBBRR_120534_543210')
            sage: T.to_string()
            'RRBBRR_120534_543210'
            sage: VeeringTriangulation.from_string(T.to_string()) == T
            True
        """
        cols, fps, eps,  = s.split('_')
        n = len(cols)
        fp = perm_from_base64_str(fps, n)
        assert perm_base64_str(fp) == fps
        ep = perm_from_base64_str(eps, n)
        assert perm_base64_str(ep) == eps
        cols = array('i', [colour_from_char(c) for c in cols])
        return VeeringTriangulation.from_face_edge_perms(cols, fp, ep, mutable=mutable, check=check)

    def forgot_forward_flippable_colour(self, folded=True):
        r"""
        Make purple the colour of forward flippable edges.

        The result of this operation is entirely determined by the vertical
        train track (as only forward flippable edges have an undefined colour).

        EXAMPLES::

            sage: from veerer import VeeringTriangulation
            sage: t = VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBBBRBBRB", mutable=True)
            sage: t.forgot_forward_flippable_colour()
            sage: t
            VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBPBRBPRB")

            sage: t = VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBBBRBBRB", mutable=False)
            sage: t.forgot_forward_flippable_colour()
            Traceback (most recent call last):
            ...
            ValueError: immutable veering triangulation; use a mutable copy instead
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        ep = self._ep
        for e in self.forward_flippable_edges(folded=folded):
            self._colouring[e] = self._colouring[ep[e]] = PURPLE

    def forgot_backward_flippable_colour(self):
        r"""
        Make green the colour of backward flippable edges.

        The result of this operation is entirely determined by the horizontal
        train track (as only forward flippable edges have an undefined colour).

        EXAMPLES::

            sage: from veerer import VeeringTriangulation

            sage: t = VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBBBRBBRB", mutable=True)
            sage: t.forgot_backward_flippable_colour()
            sage: t
            VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RGBBRGBRB")

            sage: t = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB", mutable=True)
            sage: t.forgot_backward_flippable_colour()
            sage: t._check()

            sage: t = VeeringTriangulation("(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)", "RRRRRRBBBBBBBBBBBB", mutable=True)
            sage: t.forgot_backward_flippable_colour()
            sage: t._check()

            sage: t = VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBBBRBBRB", mutable=False)
            sage: t.forgot_backward_flippable_colour()
            Traceback (most recent call last):
            ...
            ValueError: immutable veering triangulation; use a mutable copy instead
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        ep = self._ep
        for e in self.backward_flippable_edges():
            self._colouring[e] = self._colouring[ep[e]] = GREEN

    def __eq__(self, other):
        r"""
        Equality test.
        """
        if type(self) != type(other):
            raise TypeError
        return self._n == other._n and self._fp == other._fp and self._ep == other._ep and self._colouring == other._colouring

    def __ne__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n != other._n or self._fp != other._fp or self._ep != other._ep or self._colouring != other._colouring

    def _richcmp_(self, other, op):
        c = (self._n > other._n) - (self._n < other._n)
        if c:
            return rich_to_bool(op, c)

        c = (self._colouring > other._colouring) - (self._colouring < other._colouring)
        if c:
            return rich_to_bool(op, c)

        c = (self._fp > other._fp) - (self._fp < other._fp)
        if c:
            return rich_to_bool(op, c)

        c = (self._ep > other._ep) - (self._ep < other._ep)
        return rich_to_bool(op, c)

    def to_core(self, slope=VERTICAL):
        r"""
        Change the colour of each forward (reps. backward) flippable edge in
        PURPLE if ``slope`` is ``VERTICAL`` (resp ``HORIZONTAL``)

        EXAMPLES::

            sage: from veerer import *

            sage: V = VeeringTriangulation("(0,1,2)", 'RRB')
            sage: V.to_core()
            sage: V
            VeeringTriangulation("(0,1,2)", "RPB")
            sage: V.forward_flippable_edges()
            [1]
        """
        if slope == VERTICAL:
            COLOR = PURPLE
        elif slope == HORIZONTAL:
            COLOR = GREEN
        else:
            raise ValueError("slope should either be HORIZONTAL or VERTICAL")

        n = self.num_edges()
        ep = self._ep
        for e in range(n):
            if ep[e] < e:
                raise ValueError("veering triangulation not in canonical form")
            if slope == VERTICAL:
                if self.is_forward_flippable(e, check=False):
                    self._colouring[e] = self._colouring[ep[e]] = COLOR
            else:
                if self.is_backward_flippable(e, check=False):
                    self._colouring[e] = self._colouring[ep[e]] = COLOR

    def copy(self, mutable=None):
        r"""
        Return a copy of this coloured triangulation.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], "RRB", mutable=True)
            sage: S1 = T.copy()
            sage: S2 = T.copy()
            sage: T == S1 == S2
            True
            sage: S1.flip(1,BLUE)
            sage: T == S1
            False
            sage: T == S2
            True

        TESTS::

            sage: from veerer import VeeringTriangulation
            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB", mutable=True)
            sage: U = T.copy(mutable=False)
            sage: _ = hash(U)
        """
        if mutable is None:
            mutable = self._mutable

        if not self._mutable and not mutable:
            # avoid copies of immutable objects
            if type(self) is VeeringTriangulation:
                return self
            else:
                T = VeeringTriangulation.__new__(VeeringTriangulation)
                T._n = self._n
                T._vp = self._vp
                T._ep = self._ep
                T._fp = self._fp
                T._colouring = self._colouring
                T._mutable = mutable

                return T
        else:
            T = VeeringTriangulation.__new__(VeeringTriangulation)
            T._n = self._n
            T._vp = self._vp[:]
            T._ep = self._ep[:]
            T._fp = self._fp[:]
            T._colouring = self._colouring[:]
            T._mutable = mutable

            return T

    def triangulation(self, mutable=False):
        r"""
        Return the underlying triangulation.

        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], "RRB")
            sage: T.triangulation()
            Triangulation("(0,1,2)(~2,~0,~1)")
        """
        return Triangulation.copy(self, mutable)

    def _colouring_string(self, short=False):
        n = self.num_half_edges()
        ep = self._ep
        if short:
            return ''.join(colour_to_char(self._colouring[e]) for e in range(n) if e <= ep[e])
        else:
            return ''.join(colour_to_char(self._colouring[e]) for e in range(n))

    def __str__(self):
        r"""
        TESTS::

            sage: from veerer import *

            sage: VeeringTriangulation("(0,1,2)", [RED, RED, BLUE])
            VeeringTriangulation("(0,1,2)", "RRB")
        """
        return "VeeringTriangulation(\"{}\", \"{}\")".format(
                perm_cycle_string(self._fp, False, self._n, self._ep),
                self._colouring_string(short=True))

    def __repr__(self):
        return str(self)

    # TODO: this duplicates the method edge_colour
    def colour(self, e):
        e = int(e)
        return self._colouring[e]

    def angles(self):
        r"""
        Return the list of angles (divided by \pi).

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(-3, 1, -1), (-2, 0, 2)], [RED, BLUE, BLUE])
            sage: T.angles()
            [2]

            sage: T = VeeringTriangulation([(-6, 2, -2), (-5, 1, 5), (-4, 0, 4), (-3, 3, -1)],
            ....:                           [RED, RED, BLUE, BLUE, BLUE, BLUE])
            sage: T.angles()
            [2, 2]

            sage: T = VeeringTriangulation([(-12, 4, -4), (-11, -1, 11), (-10, 0, 10),
            ....:                            (-9, 9, 1), (-8, 8, -2), (-7, 7, 2),
            ....:                            (-6, 6, -3), (-5, 5, 3)],
            ....:   [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE])
            sage: sorted(T.angles())
            [1, 1, 1, 1, 1, 3]

        Some examples with purple edges::

            sage: t = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB", mutable=True)
            sage: t.forgot_forward_flippable_colour()
            sage: t.angles()
            [6]

            sage: fp = "(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)"
            sage: cols = "RRRRRRBBBBBBBBBBBB"
            sage: t = VeeringTriangulation(fp, cols, mutable=True)
            sage: t.forgot_forward_flippable_colour()
            sage: t.angles()
            [3, 3, 3, 3]
        """
        n = self.num_half_edges()
        angles = []
        seen = [False] * n
        vp = self.vertex_permutation(copy=False)
        for e in range(n):
            if seen[e]: continue

            a = 0
            col = self._colouring[e]
            # NOTE: we count by multiples of pi/2 and then divide by 2
            while not seen[e]:
                seen[e] = True
                ee = vp[e]
                ccol = self._colouring[ee]
                a += bool(col != ccol and (((col & (BLUE|RED)) and (ccol & (BLUE|RED))) or
                                           (col & (GREEN|PURPLE))))
                e = ee
                col = ccol
            a += col != ccol

            angles.append(a // 2)

        return angles + [1] * self.num_folded_edges()

    def is_abelian(self, certificate=False):
        r"""
        Return whether this coloured triangulation is Abelian.

        EXAMPLES:

            sage: from veerer import *

            sage: T = VeeringTriangulation([(-3, 1, -1), (-2, 0, 2)], [RED, BLUE, BLUE])
            sage: T.is_abelian()
            True

            sage: T = VeeringTriangulation([(-6, 2, -2), (-5, 1, 5), (-4, 0, 4), (-3, 3, -1)],
            ....:       [RED, RED, BLUE, BLUE, BLUE, BLUE])
            sage: T.is_abelian()
            True
            sage: T.is_abelian(certificate=True)
            (True, [True, True, False, False, False, False, True, True, True, True, False, False])
            sage: T = VeeringTriangulation([(-12, 4, -4), (-11, -1, 11), (-10, 0, 10),
            ....:      (-9, 9, 1), (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)],
            ....:      [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE])
            sage: T.is_abelian()
            False

        Examples with purple edges::

            sage: t = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB", mutable=True)
            sage: t.forgot_forward_flippable_colour()
            sage: t.is_abelian()
            True

            sage: fp = "(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)"
            sage: cols = "RRRRRRBBBBBBBBBBBB"
            sage: t = VeeringTriangulation(fp, cols, mutable=True)
            sage: t.forgot_forward_flippable_colour()
            sage: t.is_abelian()
            False

        Examples with green edges::

            sage: fp = "(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)"
            sage: cols = "RRRBBBBBB"
            sage: t = VeeringTriangulation(fp, cols, mutable=True)
            sage: t.forgot_backward_flippable_colour()
            sage: t.is_abelian()
            True

            sage: fp = "(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)"
            sage: cols = "RRRRRRBBBBBBBBBBBB"
            sage: t = VeeringTriangulation(fp, cols, mutable=True)
            sage: t.forgot_backward_flippable_colour()
            sage: t.is_abelian()
            False
        """
        ep = self._ep
        vp = self._vp
        cols = self._colouring

        # Perform a BFS. We store an orientation for each edge that specifies
        # its orientation in R^2: True if x > 0 or (x == 0 and y > 0)
        #                         False if x < 0 or (x == 0 and y < 0)
        oris = [None] * self._n
        oris[0] = True
        q = [0]
        while q:
            e = e0 = q.pop()
            o = oris[e0]
            f = vp[e]
            while f != e0:
                if (cols[e] == RED and cols[f] == BLUE) or (cols[e] == GREEN):
                    o = not o

                oris[f] = o
                if oris[ep[f]] is None:
                    q.append(ep[f])
                    oris[ep[f]] = not o
                elif oris[ep[f]] == oris[f]:
                    return (False, None) if certificate else False
                e, f = f, vp[f]

            if (cols[e] == RED and cols[f] == BLUE) or cols[e] == GREEN:
                o = not o
            if oris[e0] != o:
                return (False, None) if certificate else False

        return (True, oris) if certificate else True

    def abelian_cover(self, mutable=None):
        r"""
        Return the orientation double cover of this veering triangulation.

        EXAMPLES::

            sage: from veerer import VeeringTriangulation

            sage: V = VeeringTriangulation("(0,1,2)", "BBR")
            sage: V.abelian_cover().stratum()                   # optional - surface_dynamics
            H_1(0)

            sage: V = VeeringTriangulation("(0,2,3)(1,4,~0)(5,6,~1)", "BRRBBBB")
            sage: A = V.abelian_cover()
            sage: A.stratum()                                   # optional - surface_dynamics
            H_2(2)

            sage: W = V.copy(mutable=True)
            sage: W.forgot_forward_flippable_colour()
            sage: B = A.copy(mutable=True)
            sage: B.forgot_forward_flippable_colour()
            sage: assert W.abelian_cover() == B

            sage: W = V.copy(mutable=True)
            sage: W.forgot_backward_flippable_colour()
            sage: B = A.copy(mutable=True)
            sage: B.forgot_backward_flippable_colour()
            sage: assert W.abelian_cover() == B
        """
        n = self._n
        ep = self._ep

        # edge permutation
        ep_cov = [None] * (2*n)
        for e in range(n):
            E = ep[e] + n
            ep_cov[e] = E
            ep_cov[E] = e

        # face permutation
        seen = [False] * n
        fp_cov = [None] * (2*n)
        for a in range(n):
            if seen[a]:
                continue
            col, a, b, c = self.triangle(a)
            seen[a] = seen[b] = seen[c] = True
            a1 = a
            a2 = a+n
            if col == PURPLE or col == BLUE:
                b1 = b+n; b2 = b
                c1 = c+n; c2 = c
            elif col == GREEN:
                b1 = b+n; b2 = b
                c1 = c; c2 = c+n
            elif col == RED:
                b1 = b+n; b2 = b
                c1 = c; c2 = c+n
            elif col == RED|PURPLE|GREEN:
                b1 = b+n; b2 = b
                c1 = c+n; c2 = c
            elif col == BLUE|PURPLE|GREEN:
                b1 = b; b2 = b+n
                c1 = c+n; c2 = c
            else:
                raise RuntimeError
            fp_cov[a1] = b1; fp_cov[b1] = c1; fp_cov[c1] = a1
            fp_cov[a2] = b2; fp_cov[b2] = c2; fp_cov[c2] = a2

        colouring_cov = self._colouring * 2
        # TODO: remove the check argument
        vt = self.from_face_edge_perms(colouring_cov, array('i', fp_cov), array('i', ep_cov), mutable=True, check=True)
        vt.relabel(vt._relabelling_from(0), check=False)
        if not mutable:
            vt.set_immutable()
        return vt

    def stratum(self):
        r"""
        Return the Abelian or quadratic stratum of this coloured triagulation.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.stratum()  # optional - surface_dynamics
            H_1(0)

            sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
            sage: cols = 'BRBRBBBRBR'
            sage: T = VeeringTriangulation(fp, cols)
            sage: T.stratum()  # optional - surface_dynamics
            Q_1(1^2, -1^2)

        Some examples with purple edges::

            sage: from surface_dynamics import AbelianStratum, QuadraticStratum  # optional - surface_dynamics
            sage: C = AbelianStratum(2)                                          # optional - surface_dynamics
            sage: t = VeeringTriangulation.from_stratum(C, mutable=True)         # optional - surface_dynamics
            sage: t.forgot_forward_flippable_colour()                            # optional - surface_dynamics
            sage: t.stratum()                                                    # optional - surface_dynamics
            H_2(2)

            sage: C = QuadraticStratum(1,1,1,1)                                  # optional - surface_dynamics
            sage: t = VeeringTriangulation.from_stratum(C, mutable=True)         # optional - surface_dynamics
            sage: t.forgot_forward_flippable_colour()                            # optional - surface_dynamics
            sage: t.stratum()                                                    # optional - surface_dynamics
            Q_2(1^4)
        """
        from .features import surface_dynamics_feature
        surface_dynamics_feature.require()

        A = self.angles()
        if any(a%2 for a in A) or not self.is_abelian():
            from surface_dynamics import QuadraticStratum
            return QuadraticStratum([(a-2) for a in A])
        else:
            from surface_dynamics import AbelianStratum
            return AbelianStratum([(a-2)/2 for a in A])

    def dimension(self):
        r"""
        Return the dimension of the ambient stratum of Abelian or quadratic differential.
        """
        # each folded edge gives a pole
        return 2 * self.genus() - 2 + self.num_vertices() + self.num_folded_edges() + self.is_abelian()

    stratum_dimension = dimension

    def colours_about_edge(self, e):
        e = int(e)
        return [self._colouring[f] for f in self.square_about_edge(e)]

    def alternating_square(self, e):
        r"""
        Return whether there is an alternating square around the edge ``e``.
        """
        e = int(e)
        colours = self.colours_about_edge(e)
        if any(colours[f] == GREEN or colours[f] == PURPLE for f in range(4)):
            return False
        return all(colours[f] != colours[(f+1) % 4] for f in range(4))

    # TODO: change the names
    # TODO: compute it combinatorially, not via polytope computation
    # (these are the cycles and, if not orientable, the dumbbells)
    def vertex_cycles(self, slope=VERTICAL, backend=None):
        r"""
        Return the rays of the train-track polytope.

        EXAMPLES::

            sage: from veerer import VeeringTriangulation, HORIZONTAL, VERTICAL
            sage: V = VeeringTriangulation([(0,1,2), (-1,-2,-3)], "RRB")
            sage: sorted(V.vertex_cycles())
            [[0, 1, 1], [1, 1, 0]]
            sage: sorted(V.vertex_cycles(HORIZONTAL))
            [[1, 0, 1], [1, 1, 0]]
        """
        return self.train_track_polytope(slope, backend=backend).rays()

    def branches(self, slope=VERTICAL):
        r"""
        Return a 3-tuple made of the lists of respectively the small, mixed
        and large edges of the underlying train track.

        INPUT:

        - ``slope`` (optional) - either ``HORIZONTAL`` or ``VERTICAL``

        EXAMPLES::

            sage: from veerer import *
            sage: vt = VeeringTriangulation.from_string('RBBBBRBRBRBBBBRBRBRBBBBR_k509clabdjgfie876m4321nh_nmlkjihgfedcba9876543210')
            sage: vt.branches()
            ([0, 1, 5, 7, 9, 10], [2, 3, 8, 11], [4, 6])
            sage: vt.branches(HORIZONTAL)
            ([0, 4, 5, 6, 7, 9], [2, 3, 8, 11], [1, 10])
        """
        n = self.num_half_edges()
        ne = self.num_edges()
        fp = self.face_permutation(copy=False)
        ep = self.edge_permutation(copy=False)
        w = [0] * ne
        seen = [False] * n
        if slope == VERTICAL:
            left = BLUE
            right = RED
        elif slope == HORIZONTAL:
            left = RED
            right = BLUE
        else:
            raise ValueError('slope must be HORIZONTAL or VERTICAL')
        for e in range(n):
            if seen[e]:
                continue

            # go to the right transition
            a = e
            while self.colour(a) != left or self.colour(fp[a]) != right:
                a = fp[a]
            b = fp[a]
            c = fp[b]
            w[self._norm(c)] += 1  # half large
            seen[a] = seen[b] = seen[c] = True

        small = []
        mixed = []
        large = []
        for e in range(ne):
            if ep[e] == e:  # folded edge
                w[e] *= 2
            if w[e] == 0:
                small.append(e)
            elif w[e] == 1:
                mixed.append(e)
            elif w[e] == 2:
                large.append(e)
            else:
                raise RuntimeError

        return (small, mixed, large)


    def is_flippable(self, e, check=True):
        r"""
        Return whether the edge ``e`` can be flipped.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.is_flippable(0)
            True
            sage: T.is_flippable(1)
            True
            sage: T.is_flippable(2)
            False
        """
        if check:
            e = self._check_half_edge(e)
        return Triangulation.is_flippable(self, e) and self.alternating_square(e)

    def is_forward_flippable(self, e, check=True):
        r"""
        Return whether one can perform a forward flip to ``e``.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: [T.is_forward_flippable(e) for e in range(3)]
            [False, True, False]

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [GREEN, RED, BLUE])
            sage: [T.is_forward_flippable(e) for e in range(3)]
            [False, True, True]

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [PURPLE, BLUE, RED])
            sage: [T.is_forward_flippable(e) for e in range(3)]
            [True, False, False]
        """
        if check:
            e = self._check_half_edge(e)
        if self._colouring[e] == GREEN or not Triangulation.is_flippable(self, e):
            return False
        if self._colouring[e] == PURPLE:
            return True
        ca, cb, cc, cd = self.colours_about_edge(e)
        return bool(ca & (BLUE | GREEN)) and bool(cb & (RED | GREEN)) and bool(cc & (BLUE | GREEN)) and bool(cd & (RED | GREEN))

    def is_backward_flippable(self, e, check=True):
        r"""
        Return whether one can perform a backward flip to ``e``.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: [T.is_backward_flippable(e) for e in range(3)]
            [True, False, False]

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [GREEN, RED, BLUE])
            sage: [T.is_backward_flippable(e) for e in range(3)]
            [True, False, False]

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [PURPLE, BLUE, RED])
            sage: [T.is_backward_flippable(e) for e in range(3)]
            [False, True, True]
        """
        if check:
            e = self._check_half_edge(e)
        if self._colouring[e] == PURPLE or not Triangulation.is_flippable(self, e):
            return False
        if self._colouring[e] == GREEN:
            return True
        ca, cb, cc, cd = self.colours_about_edge(e)
        return bool(ca & (RED | PURPLE)) and bool(cb & (BLUE | PURPLE)) and bool(cc & (RED | PURPLE)) and bool(cd & (BLUE | PURPLE))

    def forward_flippable_edges(self, folded=True):
        r"""
        Return the set of forward flippable edges.

        INPUT:

        - ``folded`` - boolean (default ``True``) - whether to include folded edges

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.forward_flippable_edges()
            [1]

            sage: T = VeeringTriangulation("(0,1,2)", [RED, RED, BLUE])
            sage: T.forward_flippable_edges()
            [1]

        TESTS::

            sage: from veerer.permutation import perm_random
            sage: from veerer.veering_triangulation import VeeringTriangulation
            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], "RRB", mutable=True)
            sage: for _ in range(10):
            ....:     rel = perm_random(6)
            ....:     T.relabel(rel)
            ....:     assert len(T.forward_flippable_edges()) == 1
        """
        ep = self._ep
        n = self._n
        if folded:
            return [e for e in range(n) if e <= ep[e] and self.is_forward_flippable(e, check=False)]
        return [e for e in range(n) if e < ep[e] and self.is_forward_flippable(e, check=False)]

    def backward_flippable_edges(self, folded=True):
        r"""
        Return the list of backward flippable edges.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.backward_flippable_edges()
            [0]

            sage: T = VeeringTriangulation("(0,1,2)", [RED, RED, BLUE])
            sage: T.backward_flippable_edges()
            [0]
        """
        ep = self._ep
        n = self._n
        if folded:
            return [e for e in range(n) if e <= ep[e] and self.is_backward_flippable(e, check=False)]
        return [e for e in range(n) if e <= ep[e] and self.is_backward_flippable(e, check=False)]

    def purple_edges(self, folded=True):
        r"""
        EXAMPLES::

            sage: from veerer import VeeringTriangulation
            sage: t = VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBPBRBPRB")
            sage: t.purple_edges()
            [2, 6]
        """
        ep = self._ep
        n = self._n
        if folded:
            return [e for e in range(n) if e <= ep[e] and self._colouring[e] == PURPLE]
        return [e for e in range(n) if e < ep[e] and self._colouring[e] == PURPLE]

    def mostly_sloped_edges(self, slope):
        if slope == HORIZONTAL:
            return self.forward_flippable_edges()
        elif slope == VERTICAL:
            return self.backward_flippable_edges()
        else:
            raise ValueError

    def relabel(self, p, check=True):
        r"""
        Relabel inplace this veering triangulation according to the permutation ``p``.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RBB", mutable=True)
            sage: T.relabel([0,1,3,2,5,4])
            sage: T
            VeeringTriangulation("(0,1,~2)(2,~0,~1)", "RBB")
            sage: T._check()

        Composing relabellings and permutation composition::

            sage: from veerer.permutation import perm_compose, perm_random_centralizer
            sage: fp = "(0,16,~15)(1,19,~18)(2,22,~21)(3,21,~20)(4,20,~19)(5,23,~22)(6,18,~17)(7,17,~16)(8,~1,~23)(9,~2,~8)(10,~3,~9)(11,~4,~10)(12,~5,~11)(13,~6,~12)(14,~7,~13)(15,~0,~14)"
            sage: cols = "RRRRRRRRBBBBBBBBBBBBBBBB"
            sage: T0 = VeeringTriangulation(fp, cols)
            sage: for _ in range(10):
            ....:     p1 = perm_random_centralizer(T0.edge_permutation(copy=False))
            ....:     p2 = perm_random_centralizer(T0.edge_permutation(copy=False))
            ....:     T1 = T0.copy(mutable=True)
            ....:     T1.relabel(p1)
            ....:     T1.relabel(p2)
            ....:     T2 = T0.copy(mutable=True)
            ....:     T2.relabel(perm_compose(p1, p2))
            ....:     assert T1  == T2

        TESTS:

        This example used to be wrong::

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE], mutable=True)
            sage: T.relabel([1,5,0,2,4,3])
            sage: T.edge_colour(0) == BLUE
            True
            sage: T.edge_colour(1) == RED
            True
            sage: T._check()

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE], mutable=True)
            sage: from veerer.permutation import perm_random
            sage: for _ in range(10):
            ....:     r = perm_random(6)
            ....:     T.relabel(r)
            ....:     T._check()
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        n = self._n
        ep = self._ep
        if check and not perm_check(p, n):
            p = perm_init(p, n, ep)
            if not perm_check(p, n):
                raise ValueError('invalid relabelling permutation')

        Triangulation.relabel(self, p, check=False)
        perm_on_list(p, self._colouring)

    def _automorphism_good_starts(self, has_purple=None, has_green=None):
        r"""
        Start at a RED before a BLUE.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)", "RRB")
            sage: len(T._automorphism_good_starts())
            1
            sage: t = VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBPBRBPRB")
            sage: len(t._automorphism_good_starts())
            2

            sage: T = VeeringTriangulation("(0,1,2)(3,4,~1)(5,6,~4)", "RBGGRBG")
            sage: T._automorphism_good_starts()
            [2, 3, 6]
        """
        starts = []
        best_word = None

        n = self._n
        vp = self._vp
        cols = self._colouring

        if has_purple is None:
            has_purple = False
            for e in range(n):
                if self._colouring[e] == PURPLE:
                    has_purple = True
        if has_green is None:
            has_green = False
            for e in range(n):
                if self._colouring[e] == GREEN:
                    has_green = True

        # first run: compare edges using colours and vertex permutation orbits
        # FIX: Take advantage of folded edges, if they are present?
        for e in range(n):
            f = vp[e]

            # if there is any purple start from there
            # otherwise start from red -> blue transition
            if has_purple:
                if cols[e] != PURPLE:
                    continue
            elif has_green:
                if cols[e] != GREEN:
                    continue
            elif cols[e] != RED or cols[f] != BLUE:
                continue

            # build the word
            w = []
            g = f
            while True:

                # go along the v permutation and write a word in BLUE/RED
                # run through all the BLUE
                n = 0
                while cols[g] == BLUE:
                    g = vp[g]
                    n += 1
                w.append(n)
                # possibly run through a GREEN
                if cols[g] == GREEN:
                    w.append(-1)
                    g = vp[g]
                # run through all the RED
                while cols[g] == RED:
                    g = vp[g]
                    n += 1
                w.append(n)
                # possibly run through a PURPLE
                if cols[g] == PURPLE:
                    w.append(-1)
                    g = vp[g]

                if g == f or \
                   (best_word is not None and len(w) > len(best_word)):
                    break

            if best_word is None or \
               (len(w) < len(best_word) or (len(w) == len(best_word) and w < best_word)):
                starts = [e]
                best_word = w
            elif w == best_word:
                starts.append(e)

        assert starts

        if len(starts) == 1:
            return starts

        # try to break ties using face permutation orbits
        # (because the start edge e is RED before BLUE or PURPLE before BLUE (around the vertex)
        # the only colour we do not know is fp[e]
        fp = self._fp
        best_word = None
        for e in starts:
            w = cols[fp[e]]
            if best_word is None or w < best_word:
                new_starts = [e]
                best_word = w
            elif w == best_word:
                new_starts.append(e)

#        if len(new_starts) < len(starts):
#            print('1 - new_starts got better: was %s and is now %s' %(len(starts), len(new_starts)))
        starts = new_starts
#        if len(starts) == 1:
#            return starts
#
#        # try to break using vef orbits
#        ep = self._triangulation.edge_permutation(copy=False)
#        best_word = None
#        for e in starts:
#            f = vp[ep[fp[e]]]
#            w = [cols[e]]
#            while f != e:
#                w.append(cols[f])
#                f = vp[ep[fp[f]]]
#
#            if best_word is None or w < best_word:
#                new_starts = [e]
#                best_word = w
#            elif w == best_word:
#                new_starts.append(e)
#
##        if len(new_starts) < len(starts):
##            print('2 - new_starts got better: was %s and is now %s' %(len(starts), len(new_starts)))
#        starts = new_starts

        return starts

    def edge_colour(self, e):
        return self._colouring[e]

    def set_edge_colour(self, e, col):
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        if self._colouring[e] != PURPLE and self._colouring[e] != GREEN:
            raise ValueError("only PURPLE and GREEN edges could be changed colours")
        if col != BLUE and col != RED:
            raise ValueError("the new colour 'col' must be RED or BLUE")
        E = self._ep[e]
        self._colouring[e] = self._colouring[E] = col

    def set_random_colours(self):
        r"""
        Set random colours to the GREEN and PURPLE edges.
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        ep = self._ep
        recolour = []
        cols = [BLUE, RED]
        for e, col in enumerate(self._colouring):
            if col == GREEN or col == PURPLE:
                E = ep[e]
                shuffle(cols)
                oldcol = self._colouring[e]
                self._colouring[e] = self._colouring[E] = cols[0]
                if not self.edge_has_curve(e):
                    self._colouring[e] = self._colouring[E] = cols[1]
                    assert self.edge_has_curve(e)
                recolour.append((e, oldcol))
                if e != E:
                    recolour.append((E, oldcol))

    def set_colour(self, col):
        r"""
        Set all GREEN or PURPLE edges to ``col``.

        The colour ``col`` must be RED or BLUE.
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        if col != RED and col != BLUE:
            raise ValueError("'col' must be RED or BLUE")

        for e in range(self._n):
            if self._colouring[e] == GREEN or self._colouring[e] == PURPLE:
                self._colouring[e] = col

    def best_relabelling(self, all=False):
        r"""
        Return a pair ``(r, (cols, fp, ep))`` where the triple ``(cols, fp,
        ep)`` is the data of the canonical labelling of this veering
        triangulation and ``r`` is a relabelling that brings this triangulation
        to the canonical one.

        EXAMPLES::

            sage: from veerer import *
            sage: from veerer.permutation import perm_random_centralizer
            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: V = VeeringTriangulation(fp, cols, mutable=True)
            sage: r, (cols, fp, ep) = V.best_relabelling()
            sage: for _ in range(10):
            ....:     p = perm_random_centralizer(V.edge_permutation(copy=False))
            ....:     V.relabel(p)
            ....:     r2, (cols2, fp2, ep2) = V.best_relabelling()
            ....:     assert cols2 == cols and fp2 == fp and ep2 == ep
        """
        best = None
        if all:
            relabellings = []

        for start_edge in self._automorphism_good_starts():
            relabelling = self._relabelling_from(start_edge)

            # relabelled data
            fp = perm_conjugate(self._fp, relabelling)
            ep = perm_conjugate(self._ep, relabelling)
            cols = self._colouring[:]
            perm_on_list(relabelling, cols)

            T = (cols, fp, ep)
            if best is None or T < best:
                best = T
                best_relabelling = relabelling
                if all:
                    del relabellings[:]
                    relabellings.append(relabelling)
            elif all and T == best:
                relabellings.append(relabelling)

        return (relabellings, best) if all else (best_relabelling, best)

    def rotate(self):
        r"""
        Apply the pi/2 rotation.

        EXAMPLES::

            sage: from veerer import *

            sage: faces = "(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)"
            sage: cols = [BLUE,RED,RED,BLUE,RED,RED]
            sage: T = VeeringTriangulation(faces, cols)
            sage: T.rotate()
            sage: T
            VeeringTriangulation("(0,1,2)(3,4,5)(~5,~3,~1)(~4,~2,~0)", "RBBRBB")
            sage: T._check()

            sage: fp = "(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)"
            sage: cols = "RRRRRRBBBBBBBBBBBB"
            sage: T0 = VeeringTriangulation(fp, cols)
            sage: T = T0.copy(mutable=True)
            sage: T.rotate()
            sage: T.conjugate()
            sage: S = T0.copy(mutable=True)
            sage: S.conjugate()
            sage: S.rotate()
            sage: S == T
            True

        Check that PURPLE edges are mapped to GREEN::

            sage: T = VeeringTriangulation("(0,1,2)(3,4,5)(~5,~3,~1)(~4,~2,~0)", "BRPBRP", mutable=True)
            sage: T.rotate()
            sage: T
            VeeringTriangulation("(0,1,2)(3,4,5)(~5,~3,~1)(~4,~2,~0)", "RBGRBG")
            sage: T.rotate()
            sage: T
            VeeringTriangulation("(0,1,2)(3,4,5)(~5,~3,~1)(~4,~2,~0)", "BRPBRP")
        """
        for i, col in enumerate(self._colouring):
            if col == RED:
                self._colouring[i] = BLUE
            elif col == BLUE:
                self._colouring[i] = RED
            elif col == PURPLE:
                self._colouring[i] = GREEN
            else:
                assert col == GREEN
                self._colouring[i] = PURPLE

    def conjugate(self):
        r"""
        Conjugate this triangulation.

        EXAMPLES::

            sage: from veerer import *

            sage: faces = "(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)"
            sage: cols = [BLUE,RED,RED,BLUE,RED,RED]
            sage: T = VeeringTriangulation(faces, cols, mutable=True)
            sage: T.conjugate()
            sage: T
            VeeringTriangulation("(0,2,4)(1,3,5)(~5,~4,~3)(~2,~1,~0)", "RBBRBB")
            sage: T._check()

            sage: T = VeeringTriangulation(faces, cols, mutable=False)
            sage: T.conjugate()
            Traceback (most recent call last):
            ...
            ValueError: immutable veering triangulation; use a mutable copy instead
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        Triangulation.conjugate(self)
        transp = {RED: BLUE, BLUE: RED, GREEN: GREEN, PURPLE: PURPLE}
        self._colouring = array('i', [transp[col] for col in self._colouring])

    # TODO: finish this!!
    def automorphism_quotient(self, aut):
        r"""
        Return the canonical triangulation of the quotient.

        (vertex fixed, edge fixed, faces fixed)

        EXAMPLES::

            sage: from veerer import *

            sage: p = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = 'BRBBBRRBBBBR'
            sage: T = VeeringTriangulation(p, cols)
        """
        raise NotImplementedError
        nb_verts = 0  # vertices stabilized
        nb_faces = 0  # faces stabilized
        nb_edges = 0  # edges stabilized

        n = self._n # ???
        ep = self._triangulation.edge_permutation()

        # edges
        nb_edges = sum(aut[e] == ep[e] or aut[e] == e for e in range(n))

        # vertices
        verts, edge_to_vert = perm_cycles(self._vp)
        for i,v in enumerate(verts):
            if edge_to_vert[aut[v[0]]] == i:
                # deg does not change
                nb_edges += 1
            else:
                # deg is wrapped around
                pass

        # faces
        faces, edge_to_face = perm_cycles(self._fp)
        for i,f in enumerate(faces):
            nb_faces += edge_to_face[aut[f[0]]] == i

        return (nb_verts, nb_edges, nb_faces)

    def to_string(self):
        r"""
        Serialization to string.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)", "RRB")
            sage: T.to_string()
            'RRB_120_012'

            sage: from surface_dynamics import *  # optional - surface_dynamics
            sage: T = VeeringTriangulation.from_stratum(QuadraticStratum({1:20}))  # optional - surface_dynamics
            sage: s = T.to_string()  # optional - surface_dynamics
            sage: TT = VeeringTriangulation.from_string(s)  # optional - surface_dynamics
            sage: T == TT  # optional - surface_dynamics
            True
        """
        colours = self._colouring_string(short=False)
        fp = perm_base64_str(self._fp)
        ep = perm_base64_str(self._ep)
        return colours + '_' + fp + '_' + ep

    def _non_isom_easy(self, other):
        return (Triangulation._non_isom_easy(self, other) or
                self.num_folded_edges() != other.num_folded_edges() or
                self._colouring.count(RED) != other._colouring.count(RED) or
                self._colouring.count(BLUE) != other._colouring.count(BLUE) or
                self._colouring.count(PURPLE) != other._colouring.count(PURPLE))

    def flip(self, e, col, Lx=None, Gx=None, reduced=None, check=True):
        r"""
        Flip an edge inplace.

        INPUT:

        - ``e`` - edge number

        - ``col`` - colour of the edge after the flip (ie either ``RED``, ``BLUE`` or ``GREEN``)

        - ``Lx`` - (optional) - matrix whose rows are equations in a linear subspace
          that has to be carried around

        - ``Gx`` - (optional) - matrix whose rows are generators of a linear subspace
          that has to be carried around

        - ``reduced`` -- whether to change the colours of the neighbors edges to ``PURPLE``
          if they become forward flippable.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB", mutable=True)
            sage: T.flip(1, RED)
            sage: T
            VeeringTriangulation("(0,~2,1)(2,~1,~0)", "RRB")
            sage: T.flip(0, RED)
            sage: T
            VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB")
            sage: T.flip(1, BLUE)
            sage: T
            VeeringTriangulation("(0,~2,1)(2,~1,~0)", "RBB")
            sage: T.flip(2, BLUE)
            sage: T
            VeeringTriangulation("(0,~1,~2)(1,2,~0)", "RBB")

        The same flip sequence with reduced veering triangulations (forward flippable
        edges in ``PURPLE``)::

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB", mutable=True)
            sage: T.forgot_forward_flippable_colour()
            sage: T
            VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RPB")
            sage: T.flip(1, RED)
            sage: T
            VeeringTriangulation("(0,~2,1)(2,~1,~0)", "PRB")
            sage: T.flip(0, RED)
            sage: T
            VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RPB")
            sage: T.flip(1, BLUE)
            sage: T
            VeeringTriangulation("(0,~2,1)(2,~1,~0)", "RBP")
            sage: T.flip(2, BLUE)
            sage: T
            VeeringTriangulation("(0,~1,~2)(1,2,~0)", "RPB")

        Some examples involving linear subspaces::

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: T = T.copy(mutable=True)
            sage: Gx = matrix(ZZ, [s, t])
            sage: T.flip(3, 2, Gx=Gx)
            sage: T.flip(4, 2, Gx=Gx)
            sage: T.flip(5, 2, Gx=Gx)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(0), VERTICAL)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(1), VERTICAL)

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(2, 3, 4, 5, 1, 2)
            sage: T = T.copy(mutable=True)
            sage: Gx = matrix(ZZ, [s, t])
            sage: flip_sequence = [(3, 2), (4, 1), (5, 2), (6 , 2), (5, 1), (1, 1), (5, 1)]
            sage: for e, col in flip_sequence:
            ....:     T.flip(e, col, Gx=Gx)
            ....:     T._set_switch_conditions(T._tt_check, Gx.row(0), VERTICAL)
            ....:     T._set_switch_conditions(T._tt_check, Gx.row(1), VERTICAL)
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        if check:
            if col != BLUE and col != RED and col != GREEN:
                raise ValueError("'col' must be BLUE, RED or GREEN")

            e = self._check_half_edge(e)
            if not self.is_forward_flippable(e, check=False):
                raise ValueError("half-edge e={} is not forward flippable".format(e))

        if reduced is None:
            reduced = self._colouring[e] == PURPLE

        if Lx is not None:
            raise NotImplementedError("not implemented for linear equations")
        if Gx is not None:
            a, b, c, d = self.square_about_edge(e)
            e = self._norm(e)
            a = self._norm(a)
            d = self._norm(d)

            b = self._norm(b)
            c = self._norm(c)
            assert Gx.column(e) == Gx.column(a) + Gx.column(b) == Gx.column(c) + Gx.column(d)
            if col == RED:
                # ve <- vd - va
                # e-th column becomes d-th column minus a-th column
                Gx.add_multiple_of_column(e, e, -1)
                Gx.add_multiple_of_column(e, d, +1)
                Gx.add_multiple_of_column(e, a, -1)
            elif col == BLUE:
                # ve <- va - vd
                Gx.add_multiple_of_column(e, e, -1)
                Gx.add_multiple_of_column(e, d, -1)
                Gx.add_multiple_of_column(e, a, +1)
            else:
                raise NotImplementedError('GREEN not implemented with linear subspace')

        # flip and set colour
        E = self._ep[e]
        Triangulation.flip(self, e, check=False)
        self._colouring[e] = self._colouring[E] = col

        if reduced:
            ep = self._ep
            a, b, c, d = self.square_about_edge(e)
            assert self._colouring[a] & (RED | GREEN), (a, colour_to_string(self._colouring[a]))
            assert self._colouring[b] & (BLUE | GREEN), (b, colour_to_string(self._colouring[b]))
            assert self._colouring[c] & (RED | GREEN), (c, colour_to_string(self._colouring[c]))
            assert self._colouring[d] & (BLUE | GREEN), (d, colour_to_string(self._colouring[d]))

            # assertions to be removed
            assert not self.is_forward_flippable(e)

            if col == BLUE:
                if self.is_forward_flippable(b, check=False):
                    self._colouring[b] = PURPLE
                    self._colouring[ep[b]] = PURPLE
                if d != ep[b] and self.is_forward_flippable(d, check=False):
                    self._colouring[d] = PURPLE
                    self._colouring[ep[d]] = PURPLE
                assert not self.is_forward_flippable(a)
                assert not self.is_forward_flippable(c)
            elif col == RED:
                if self.is_forward_flippable(a, check=False):
                    self._colouring[a] = PURPLE
                    self._colouring[ep[a]] = PURPLE
                if c != ep[a] and self.is_forward_flippable(c, check=False):
                    self._colouring[c] = PURPLE
                    self._colouring[ep[c]] = PURPLE
                assert not self.is_forward_flippable(b)
                assert not self.is_forward_flippable(d)
            else:
                assert col == GREEN
                # should we put all edges in a cylinder PURPLE?

    def cylinders(self, col):
        r"""
        Return the list of cylinders of colour ``col``.

        Each cylinder is given as a quadruple ``(edges, rbdry, lbdry, half)`` where

        - ``edges`` are the edges crossed by the core curve of the annulus

        - ``rbdry`` and ``lbdry`` are respectively the list of edges on the right
          and left boundaries

        - ``half`` is a boolean that is ``True`` when it is cut by two folded edges

        EXAMPLES::

            sage: from veerer import VeeringTriangulation, RED, BLUE

        The torus::

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.cylinders(RED)
            [([0, 4], [2], [3], False)]
            sage: T.cylinders(BLUE)
            []

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "BRB")
            sage: T.cylinders(BLUE)
            [([2, 5], [1], [4], False)]
            sage: T.cylinders(RED)
            []

        Some examples in Q(4, -1^8)::

            sage: T = VeeringTriangulation("(0,1,2)", "BBR")
            sage: T.cylinders(BLUE)
            [([0, 1], [], [2], True)]
            sage: T.cylinders(RED)
            []

            sage: fp = "(5,4,7)(~5,3,10)(1,~0,8)(~1,~4,11)(2,6,9)(~2,0,12)"
            sage: cols = "BBBBBBBRRRRRR"
            sage: T = VeeringTriangulation(fp, cols, mutable=True)
            sage: T.cylinders(BLUE)
            [([6, 2, 0, 1, 14, 5, 3], [9, 8, 7], [12, 11, 10], True)]
            sage: T.cylinders(RED)
            []
            sage: T.forgot_forward_flippable_colour()
            sage: T.cylinders(RED)
            [([12, 15, 9], [6], [0], True),
             ([8, 1, 11], [14], [17], True),
             ([10, 13, 7], [4], [3], True)]

            sage: fp = "(0,12,3)(1,2,11)(4,6,7)(5,~2,10)(8,~5,~4)(9,~1,~0)"
            sage: cols = "BBRRRBBBBRRRB"
            sage: T = VeeringTriangulation(fp, cols)
            sage: T.cylinders(BLUE)
            [([6, 7], [], [4], True)]
            sage: T.cylinders(RED)
            [([10, 15, 11], [5], [1], True)]

            sage: fp = '(0,~5,4)(1,~3,2)(3,5,~4)(6,~1,~0)'
            sage: cols = 'RBBRBRB'
            sage: T = VeeringTriangulation(fp, cols)
            sage: T.cylinders(BLUE)
            [([2, 1, 6], [11], [9], True)]
            sage: T.cylinders(RED)
            []

            sage: T = VeeringTriangulation.from_string("RBRBBBRRBRBR_908274a531b6_ba9385764210")
            sage: T.cylinders(BLUE)
            [([5, 4, 3], [], [7, 2], True)]

        Torus with PURPLE edge::

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "PBR", mutable=True)
            sage: T.cylinders(RED)
            [([2, 5], [1], [4], False)]
            sage: T.cylinders(BLUE)
            [([0, 4], [2], [3], False)]

            sage: R = T.copy()
            sage: R.set_colour(RED)
            sage: R.cylinders(RED)
            [([2, 5], [1], [4], False)]

            sage: B = T.copy()
            sage: B.set_colour(BLUE)
            sage: B.cylinders(BLUE)
            [([0, 4], [2], [3], False)]

        Longer examples::

            sage: fp = "(~0,5,1)(~1,6,2)(~2,3,9)(~3,4,8)(~4,7,0)"
            sage: cols = "BBBBBRRRRR"
            sage: T = VeeringTriangulation(fp, cols, mutable=True)
            sage: T.cylinders(BLUE)
            [([0, 1, 2, 3, 4], [7, 5, 6], [9, 8], False)]
            sage: T.forgot_forward_flippable_colour()
            sage: T.cylinders(BLUE)
            [([0, 1, 2, 3, 4], [7, 5, 6], [9, 8], False)]
            sage: T.set_colour(BLUE)
            sage: T.rotate()
            sage: T.cylinders(RED)
            [([0, 1, 2, 3, 4], [7, 5, 6], [9, 8], False)]

            sage: V = VeeringTriangulation("(0,3,8)(1,~4,2)(4,5,~11)(6,~5,7)(9,~8,11)(10,~2,~1)(~10,~7,~9)(~6,~0,~3)", "PBPBRPRBRBRB")
            sage: V.cylinders(BLUE)
            [([0, 20], [8], [17], False),
             ([2, 22], [19], [10], False),
             ([5, 7, 14, 11], [4, 15], [6, 13], False)]

            sage: V = VeeringTriangulation("(0,12,~11)(1,13,~12)(2,3,14)(4,10,9)(5,~10,11)(6,~5,~17)(7,~0,~6)(8,~4,~7)(15,~13,~14)(16,17,~2)(~16,~3,~15)(~9,~8,~1)", "RRRBRRBBBBPBBBRRPB")
            sage: V.cylinders(BLUE)
            []
        """
        if col != RED and col != BLUE:
            raise ValueError("'col' must be RED or BLUE")

        opcol = RED if col == BLUE else BLUE

        n = self._n
        fp = self._fp
        ep = self._ep
        cols = self._colouring

        cylinders = []

        seen = [False] * n
        for a in range(n):
            if seen[a]:
                continue

            # triangle (a,b,c)
            b = fp[a]
            c = fp[b]
            if seen[b] or seen[c]:
                continue

            if (cols[a] == opcol) + (cols[b] == opcol) + (cols[c] == opcol) > 1:
                seen[a] = seen[b] = seen[c] = True
                continue

            # We normalize the edges (a,b,c) so that a is the one going forward
            # (each triangle in a cylinder has two "doors" and one "boundary"). We
            # choose a traversal order in the triangle and hence have two types of
            # triangles.
            #
            #  right triangle       left triangle
            #
            #       x               x---------x
            #      / \               \    b  /
            #     /b  \ -->           \c    /  -->
            #    /    a\               \  a/
            #   /  c    \               \ /
            #  x---------x               x
            #
            # We make it so we start with a right triangle.
            if cols[a] == opcol:
                a, b, c = b, c, a
            elif cols[b] == opcol:
                a, b, c = c, a, b
            assert cols[a] == col or cols[a] == PURPLE
            assert cols[b] == col or cols[b] == PURPLE
            assert cols[c] == opcol
            a0, b0, c0 = a, b, c
            RIGHT = 0
            LEFT = 1
            typ = RIGHT

            cc = []    # cycle of half-edges inside the cylinder
            rbdry = [] # right boundary
            lbdry = [] # left boundary
            half_turn = False # whether the cylinder is a folded cylinder
            cyl = True
            while not seen[a]:
                assert cols[a] == col or cols[a] == PURPLE, (typ, a, colour_to_string(cols[a]), b, colour_to_string(cols[b]), c, colour_to_string(cols[c]))
                seen[a] = seen[ep[a]] = True
                cc.append(a)
                if typ == RIGHT:
                    assert cols[c] == opcol, (typ, a, colour_to_string(cols[a]), b, colour_to_string(cols[b]), c, colour_to_string(cols[c]))
                    assert cols[b] == col or cols[b] == PURPLE, (typ, a, colour_to_string(cols[a]), b, colour_to_string(cols[b]), c, colour_to_string(cols[c]))
                    rbdry.append(c)
                else:
                    assert cols[b] == opcol, (typ, a, colour_to_string(cols[a]), b, colour_to_string(cols[b]), c, colour_to_string(cols[c]))
                    assert cols[c] == col or cols[c] == PURPLE, (typ, a, colour_to_string(cols[a]), b, colour_to_string(cols[b]), c, colour_to_string(cols[c]))
                    lbdry.append(b)

                if a == ep[a]:
                    # found a folded edge... we can not continue in this
                    # direction. Either stop or start again from the other door.
                    if half_turn:
                        cyl = True
                        break
                    half_turn = True
                    cc = [ep[x] for x in reversed(cc)]
                    lbdry.reverse()
                    rbdry.reverse()
                    rbdry, lbdry = lbdry, rbdry
                    lbdry.pop()  # a0 will be added again at the beginning of next loop
                    assert cols[c0] == opcol
                    a, b, c = b0, c0, a0
                    typ = LEFT
                    continue

                # go to triangle across the door
                a = ep[a]
                b = fp[a]
                c = fp[b]
                assert cols[a] == col or cols[a] == PURPLE, (self, a, colour_to_string(cols[a]), b, colour_to_string(cols[b]), c, colour_to_string(cols[c]))
                if cols[b] == col or cols[b] == PURPLE:
                    assert cols[c] == opcol, (self, a, colour_to_string(cols[a]), b, colour_to_string(cols[b]), c, colour_to_string(cols[c]))
                    # LEFT type, next door is  b
                    a, b, c = b, c, a
                    typ = LEFT
                elif cols[c] == col or cols[c] == PURPLE:
                    assert cols[b] == opcol, (self, a, colour_to_string(cols[a]), b, colour_to_string(cols[b]), c, colour_to_string(cols[c]))
                    # RIGHT type, next door is c
                    a, b, c = c, a, b
                    typ = RIGHT
                else:
                    cyl = False
                    break

            if cyl and (a == a0 or half_turn):
                cylinders.append((cc, rbdry, lbdry, half_turn))

        return cylinders

    def dehn_twists(self, col):
        r"""
        Return the list of Dehn twists along the cylinders of colour ``col``.

        EXAMPLES::

            sage: from veerer import VeeringTriangulation, BLUE, RED

            sage: T = VeeringTriangulation("(0,~4,5)(1,~0,6)(2,8,~1)(3,~7,~2)(4,~3,7)(~8,9,~11)(10,~5,~9)(11,~6,~10)", "BBBBBRRRRBBB")
            sage: b1, b2 = T.dehn_twists(BLUE)
            sage: b1
            VeeringFlipSequence(VeeringTriangulation("(0,~4,5)...(11,~6,~10)", "BBBBBRRRRBBB"), "1B 0B 4B 2B 1B 0B", "(0,3,1,4,~2)(2,~0,~3,~1,~4)(5)(6)(7)(8)(9)(10)(11)(~11)(~10)(~9)(~8)(~7)(~6)(~5)")
            sage: b2
            VeeringFlipSequence(VeeringTriangulation("(0,~4,5)...(11,~6,~10)", "BBBBBRRRRBBB"), "9B 10B", "(0)(1)(2)(3)(4)(5)(6)(7)(8)(9,~10,11)(10,~11,~9)(~8)(~7)(~6)(~5)(~4)(~3)(~2)(~1)(~0)")

            sage: T.rotate()
            sage: r1, r2 = T.dehn_twists(RED)
            sage: r1
            VeeringFlipSequence(VeeringTriangulation("(0,~4,5)...(11,~6,~10)", "RRRRRBBBBRRR"), "3R 4R 0R 2R 3R 4R", "(0,2,4,1,3)(5)(6)(7)(8)(9)(10)(11)(~11)(~10)(~9)(~8)(~7)(~6)(~5)(~4,~1,~3,~0,~2)")
            sage: r2
            VeeringFlipSequence(VeeringTriangulation("(0,~4,5)...(11,~6,~10)", "RRRRRBBBBRRR"), "11R 10R", "(0)(1)(2)(3)(4)(5)(6)(7)(8)(9,11,10)(~11,~10,~9)(~8)(~7)(~6)(~5)(~4)(~3)(~2)(~1)(~0)")

        A (purple) square-tiled surface corresponds to a Penner system. A
        product associated to the Dehn twists is of pseudo-Anosov type if and
        only if all twists appear at least once::

            sage: T = VeeringTriangulation("(0,~2,1)(2,~11,~3)(3,10,~4)(4,~15,~5)(5,14,~6)(6,~10,~7)(7,9,~8)(8,17,~9)(11,15,~12)(12,~17,~13)(13,~16,~14)(16,~1,~0)", "PRBPRBPRBPBRBPRPBR")
            sage: b1,b2 = T.dehn_twists(BLUE)
            sage: r1,r2,r3 = T.dehn_twists(RED)

            sage: (r1 * r3 * b1 * b2).is_pseudo_anosov()
            False
            sage: (r1 * r3 * b1 * b2 * r2).is_pseudo_anosov()
            True
        """
        if col != RED and col != BLUE:
            raise ValueError("'col' must be RED or BLUE")

        from .flip_sequence import VeeringFlipSequence
        twists = []

        opcol = BLUE if col == RED else RED
        vp = self._vp
        ep = self._ep
        cols = self._colouring
        for mid, rbdry, lbdry, half in self.cylinders(col):
            if half:
                raise NotImplementedError

            # count packets
            edges = []
            packets = []
            p = []
            for r in lbdry:
                assert cols[r] == opcol
                s = vp[r]
                del p[:]
                while cols[s] != opcol:
                    p.append(s)
                    s = vp[s]
                edges.extend(p)
                packets.append(len(p))

            # for as many times as there are saddle connection
            # in the bdry twist
            F = VeeringFlipSequence(self)
            m = len(lbdry)
            n = len(edges)
            flipsmod2 = [0] * n
            for shift in range(m):
                j = shift if col == BLUE else -shift
                for i in range(m):
                    K = range(packets[i] - 1, 0, -1) if col == BLUE else range(0, packets[i] - 1, 1)
                    for k in K:
                        l = (j + k) % n
                        F.append_flip(edges[l], col)
                        flipsmod2[l] = 1 - flipsmod2[l]
                    j += packets[i]
            r = perm_id(self._n)
            for i in range(n):
                j = (i-m)%n if col == BLUE else (i+m)%n
                e = edges[i]
                f = edges[j]
                if col == BLUE and flipsmod2[i]:
                    r[e] = ep[f]
                    r[ep[e]] = f
                else:
                    r[e] = f
                    r[ep[e]] = ep[f]
            F.append_relabelling(r)

            # TODO: remove assertion check
            assert F.start() == F.end()

            twists.append(F)

        return twists

    def is_cylindrical(self, col=None):
        r"""
        Return whether this veering triangulation is cylindrical.

        A Veering triangulation is blue cylindrical (resp red cylindrical) if
        all its triangles have two blue edges (resp red). Here a purple edge
        counts for both blue and red.

        It is purple cylindrical if all its triangle are adjacent to a purple
        edge. This is equivalent to say that both the red and blue colouring
        of the purple edges are respectively red and blue cylindrical

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB", mutable=True)
            sage: T.is_cylindrical()
            True

            sage: T.forgot_forward_flippable_colour()
            sage: T.is_cylindrical(PURPLE)
            True
            sage: T.is_cylindrical(RED)
            True
            sage: T.is_cylindrical(BLUE)
            True

            sage: T = VeeringTriangulation("(0,3,4)(1,~3,5)(2,6,~4)", "PPPBRRB")
            sage: T.is_cylindrical(PURPLE)
            True
        """
        if col is None:
            return self.is_cylindrical(PURPLE) or self.is_cylindrical(BLUE) or self.is_cylindrical(RED)
        elif col != BLUE and col != RED and col != PURPLE:
            raise ValueError("'col' must be one of BLUE, RED or PURPLE")

        n = self._n
        fp = self._fp
        cols = self._colouring
        seen = [False] * n
        for a in range(n):
            if seen[a]:
                continue
            b = fp[a]
            c = fp[b]
            if col == PURPLE:
                if cols[a] != PURPLE and cols[b] != PURPLE and cols[c] != PURPLE:
                    return False
            else:
                if (cols[a] == PURPLE) + (cols[b] == PURPLE) + (cols[c] == PURPLE) + \
                   (cols[a] == col) + (cols[b] == col) + (cols[c] == col) != 2:
                    return False
        return True

    def is_quadrangulable(self):
        k = 0
        n = self.num_edges()
        ep = self._ep
        for e in range(n):
            if not self.is_forward_flippable(e, check=False):
                continue

            k += 1 if ep[e] == e else 2

        return k == self.num_faces()

    def is_square_tiled(self, col=PURPLE):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RBR")
            sage: T.is_square_tiled(BLUE)
            False
            sage: T.is_square_tiled(RED)
            True
            sage: T.is_cylindrical()
            True

            sage: T = VeeringTriangulation("(0,1,4)(~1,2,5)(~2,~4,3)(~3,~5,~0)", "RRRRBB")
            sage: T.is_square_tiled()
            False
            sage: T.is_cylindrical()
            True

            sage: T = VeeringTriangulation("(~0,1,4)(~1,~4,2)(~2,3,5)(~3,6,0)", "RRRRBBB")
            sage: T.is_square_tiled(RED)
            True
            sage: T.is_square_tiled(BLUE)
            False
            sage: T.is_cylindrical()
            True

            sage: T = VeeringTriangulation("(0,1,2)", "RRB")
            sage: T.is_square_tiled(RED)
            True
            sage: T.is_square_tiled(BLUE)
            False
            sage: T.is_cylindrical()
            True

            sage: T = VeeringTriangulation("(1,~3,2)(~2,0,3)", "RRRB")
            sage: T.is_square_tiled(RED)
            True
            sage: T.is_square_tiled(BLUE)
            False
            sage: T.is_cylindrical()
            True
        """
        k = 0
        n = self.num_edges()
        ep = self._ep
        for e in range(n):
            if not self.is_forward_flippable(e, check=False):
                continue
            if self._colouring[e] != col:
                return False

            k += 1 if ep[e] == e else 2

        return k == self.num_faces()

    def properties_code(self):
        r"""
        Return an integer code that gathers boolean properties of this veering
        triangulation.

        EXAMPLES::

            sage: from veerer import *
            sage: from veerer.constants import properties_to_string
            sage: T = VeeringTriangulation("(0,1,8)(2,~7,~1)(3,~0,~2)(4,~5,~3)(5,6,~4)(7,~8,~6)", "BRRRRBRBR")
            sage: T.properties_code()
            17
            sage: properties_to_string(81)
            'red geometric'
        """
        from .constants import BLUE, RED, SQUARETILED, QUADRANGULABLE, GEOMETRIC

        code = 0
        if self.is_square_tiled(RED):
            code |= SQUARETILED
            code |= RED
        if self.is_square_tiled(BLUE):
            code |= SQUARETILED
            code |= BLUE
        if self.is_quadrangulable():
            code |= QUADRANGULABLE
        if self.is_cylindrical(RED):
            code |= RED
        if self.is_cylindrical(BLUE):
            code |= BLUE
        if self.is_geometric():
            code |= GEOMETRIC

        if code & BLUE and code & RED:
            raise RuntimeError("found a blue and red triangulations!")
        if code & SQUARETILED:
            if not code & BLUE and not code & RED:
                raise RuntimeError("square-tiled should be coloured")
            if not code & GEOMETRIC:
                raise RuntimeError("square-tiled should be geometric")

        return code

    def flip_back(self, e, col, Lx=None, Gx=None, check=True):
        r"""
        Flip backward an edge in place

        EXAMPLES::

            sage: from veerer import *

            sage: T0 = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T = T0.copy(mutable=True)
            sage: T.flip(1, RED)
            sage: T.flip(0, RED)
            sage: T.flip_back(0, RED)
            sage: T.flip_back(1, RED)
            sage: T == T0
            True

            sage: T.flip(1, BLUE)
            sage: T.flip(2, BLUE)
            sage: T.flip_back(2, BLUE)
            sage: T.flip_back(1, RED)
            sage: T == T0
            True

        Some examples involving linear subspaces::

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: T = T.copy(mutable=True)
            sage: Gx = matrix(ZZ, [s, t])
            sage: Gx.echelon_form()
            [1 0 0 1 1 1 1]
            [0 1 1 1 1 1 0]
            sage: T.flip(3, 2, Gx=Gx)
            sage: T.flip(4, 2, Gx=Gx)
            sage: T.flip(5, 2, Gx=Gx)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(0), VERTICAL)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(1), VERTICAL)
            sage: Gx.echelon_form()
            [ 1  0  0  1  1  1  1]
            [ 0  1  1 -1 -1 -1  0]
            sage: T.flip_back(5, 2, Gx=Gx)
            sage: T.flip_back(4, 2, Gx=Gx)
            sage: T.flip_back(3, 2, Gx=Gx)
            sage: Gx.echelon_form()
            [1 0 0 1 1 1 1]
            [0 1 1 1 1 1 0]
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        if check:
            if col != BLUE and col != RED and col != PURPLE:
                raise ValueError("'col' must be BLUE, RED or PURPLE")

            e = self._check_half_edge(e)
            if not self.is_backward_flippable(e, check=False):
                raise ValueError('half-edge e={} is not backward flippable'.format(e))

        if Lx is not None:
            raise NotImplementedError("not implemented for linear equations")

        E = self._ep[e]

        Triangulation.flip_back(self, e, check=False)
        old_col = self._colouring[e]
        self._colouring[e] = self._colouring[E] = col

        if Gx is not None:
            a, b, c, d = self.square_about_edge(e)
            e = self._norm(e)
            a = self._norm(a)
            b = self._norm(b)
            c = self._norm(c)
            d = self._norm(d)
            Gx.add_multiple_of_column(e, e, -1)
            Gx.add_multiple_of_column(e, a, +1)
            Gx.add_multiple_of_column(e, b, +1)

    def _set_switch_conditions(self, insert, x, slope=VERTICAL):
        r"""
        These are the linear parts of the train-track equations

        INPUT:

        - ``insert`` - a function to be called for each equation

        - ``x`` - variable factory (the variable for edge ``e`` is constructed
          via ``x[e]``)

        - ``slope`` - (default ``VERTICAL``) the slope of the train-track
          ``HORIZONTAL`` or ``VERTICAL``

        EXAMPLES::

            sage: from veerer import VeeringTriangulation
            sage: import ppl
            sage: vt = VeeringTriangulation("(0,1,2)(3,4,5)(6,7,8)(~8,~0,~7)(~6,~1,~5)(~4,~2,~3)", "RRBRRBRRB")
            sage: cs = ppl.Constraint_System()
            sage: x = [ppl.Variable(e) for e in range(vt.num_edges())]
            sage: vt._set_switch_conditions(cs.insert, x)
            sage: for g in cs:
            ....:     print(vector(ZZ, g.coefficients()))
            (1, -1, 1, 0, 0, 0, 0, 0, 0)
            (0, 0, 0, 1, -1, 1, 0, 0, 0)
            (0, 0, 0, 0, 0, 0, 1, -1, 1)
            (1, 0, 0, 0, 0, 0, 0, -1, 1)
            (0, 1, 0, 0, 0, -1, -1, 0, 0)
            (0, 0, 1, 1, -1, 0, 0, 0, 0)
        """
        if slope == VERTICAL:
            LAR = PURPLE
            POS = BLUE
            NEG = RED
            ZERO = GREEN
        elif slope == HORIZONTAL:
            LAR = GREEN
            POS = RED
            NEG = BLUE
            ZERO = PURPLE
        else:
            raise ValueError('bad slope parameter')

        for (i,j,k) in self.faces():
            i = self._norm(i)
            ci = self._colouring[i]
            j = self._norm(j)
            cj = self._colouring[j]
            k = self._norm(k)
            ck = self._colouring[k]

            if ci == ZERO and cj == NEG and ck == POS:
                # i is degenerate
                insert(x[i] == 0)
                insert(x[j] == x[k])
            elif cj == ZERO and ck == NEG and ci == POS:
                # j is degenerate
                insert(x[j] == 0)
                insert(x[k] == x[i])
            elif ck == ZERO and ci == NEG and cj == POS:
                # k is degenerate
                insert(x[k] == 0)
                insert(x[i] == x[j])
            elif ck == LAR or (ci == POS and cj == NEG):
                # k is large
                insert(x[k] == x[i] + x[j])
            elif ci == LAR or (cj == POS and ck == NEG):
                # i is large
                insert(x[i] == x[j] + x[k])
            elif cj == LAR or (ck == POS and ci == NEG):
                # j is large
                insert(x[j] == x[k] + x[i])
                l,s1,s2 = j,k,i
            else:
                raise ValueError('can not determine the nature of triangle (%s, %s, %s) with colors (%s, %s, %s) in %s direction' %
                                 (self._edge_rep(i), self._edge_rep(j), self._edge_rep(k),
                                  colour_to_string(ci), colour_to_string(cj), colour_to_string(ck),
                                  'horizontal' if slope == HORIZONTAL else 'vertical'))

    @staticmethod
    def _tt_check(x):
        if not x:
            raise AssertionError("does not satisfy train-track constraints")

    def train_track_switch_constraints(self, slope=VERTICAL):
        r"""
        Return the linear constraints of the train-track.

        EXAMPLES::

            sage: from veerer import VeeringTriangulation
            sage: fp = "(0,3,8)(~0,5,6)(~3,4,2)(~4,1,7)"
            sage: cols = "BBBRRRRRR"
            sage: vt = VeeringTriangulation(fp, cols)
            sage: cs = vt.train_track_switch_constraints()
            sage: cs
            {x0 + x3 - x8 == 0, x1 - x4 + x7 == 0, x2 + x3 - x4 == 0, x0 + x5 - x6 == 0}
        """
        from sage.rings.integer_ring import ZZ
        from .polyhedron.linear_expression import LinearExpressions
        L = LinearExpressions(ZZ)
        cs = ConstraintSystem()
        ne = self.num_edges()
        variables = [L.variable(e) for e in range(ne)]
        self._set_switch_conditions(cs.insert, variables, slope)
        return cs

    def _set_train_track_constraints(self, insert, x, slope, low_bound, allow_degenerations):
        r"""
        Sets the equations and inequations for train tracks.

        INPUT:

        - ``insert`` - a function to be called for each equation or inequation

        - ``x`` - variable factory (the variable for edge ``e`` is constructed
          via ``x[e]``)

        - ``slope`` - the slope of the train-track ``HORIZONTAL`` or ``VERTICAL``

        - ``low_bound`` - boolean - whether to set lower bounds to 0 or 1

        - ``allow_degenerations`` - boolean - allow to ignore the lower bound
          in appropriate situations

        EXAMPLES::

            sage: from veerer import *

            sage: T  = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, 0, False)
            sage: l
            [x0 == x1 + x2, x0 == x1 + x2, x0 >= 0, x1 >= 0, x2 >= 0]

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, 1, False)
            sage: l
            [x0 == x1 + x2, x0 == x1 + x2, x0 >= 1, x1 >= 1, x2 >= 1]

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, 3, False)
            sage: l
            [x0 == x1 + x2, x0 == x1 + x2, x0 >= 3, x1 >= 3, x2 >= 3]

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, 2, True)
            sage: l
            [x0 == x1 + x2, x0 == x1 + x2, x0 >= 2, x1 >= 0, x2 >= 2]

        This can also be used to check that a given vector satisfies the train-track
        equations::

            sage: T._set_train_track_constraints(T._tt_check, [2,1,1], HORIZONTAL, False, False)
            sage: T._set_train_track_constraints(T._tt_check, [1,1,1], HORIZONTAL, False, False)
            Traceback (most recent call last):
            ...
            AssertionError: does not satisfy train-track constraints

        Check equations with folded edges (that are "counted twice")::

            sage: T = VeeringTriangulation("(0,2,3)(~0,1,4)(~1,5,6)", [BLUE, RED, RED, BLUE, BLUE, BLUE, BLUE])
            sage: T._set_train_track_constraints(T._tt_check, [0,1,1,1,1,1,0], VERTICAL, False, False)
            sage: T._set_train_track_constraints(T._tt_check, [1,2,3,4,3,7,5], VERTICAL, False, False)
        """
        if slope == VERTICAL:
            POS = BLUE
            NEG = RED
            ZERO = GREEN
        elif slope == HORIZONTAL:
            POS = RED
            NEG = BLUE
            ZERO = PURPLE
        else:
            raise ValueError('bad slope parameter')

        low_bound = max(0, int(low_bound))
        ne = self.num_edges()
        ep = self._ep

        # switch
        self._set_switch_conditions(insert, x, slope)

        # non-negativity
        for e in range(ne):
            if ep[e] != e and ep[e] < ne:
                raise ValueError('edge permutation not in standard form')
            if self._colouring[e] == ZERO:
                # already done in switch conditions: insert(x[e] == 0)
                pass
            elif not low_bound or \
                (allow_degenerations and \
                 ((slope == HORIZONTAL and self.is_forward_flippable(e, check=False)) or \
                 (slope == VERTICAL and self.is_backward_flippable(e, check=False)))):
                insert(x[e] >= 0)
            else:
                insert(x[e] >= low_bound)

    def _set_geometric_constraints(self, insert, x, y, hw_bound=0):
        r"""
        Set the geometric constraints.

        INPUT:

        - ``insert`` - function

        - ``x``, ``y`` - variables

        - ``hw_bound`` - a nonegative number

        EXAMPLES::

            sage: from veerer import *

            sage: T  = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: y = [SR.var("y0"), SR.var("y1"), SR.var("y2")]
            sage: T._set_geometric_constraints(l.append, x, y)
            sage: l
            [x1 <= y0 + y2, y0 <= x1 + x2]
        """
        hw_bound = max(0, int(hw_bound))
        for e in self.forward_flippable_edges():
            a, b, c, d = self.square_about_edge(e)
            insert(x[self._norm(e)] <= y[self._norm(a)] + y[self._norm(d)] - hw_bound)
        for e in self.backward_flippable_edges():
            a, b, c, d = self.square_about_edge(e)
            insert(y[self._norm(e)] <= x[self._norm(a)] + x[self._norm(d)] - hw_bound)

    def _set_balance_constraints(self, insert, x, slope, homogeneous):
        r"""
        Linear constraints for the balanced polytope.
        """
        if homogeneous:
            if slope == VERTICAL:
                for eff, ebf in itertools.product(self.forward_flippable_edges(), self.backward_flippable_edges()):
                    a, b, c, d = self.square_about_edge(ebf)
                    insert(x[self._norm(b)] + x[self._norm(c)] >= x[self._norm(eff)])
            elif slope == HORIZONTAL:
                for eff, ebf in itertools.product(self.forward_flippable_edges(), self.backward_flippable_edges()):
                    a, b, c, d = self.square_about_edge(eff)
                    insert(x[self._norm(b)] + x[self._norm(c)] >= x[self._norm(ebf)])
            else:
                raise ValueError("slope must be HORIZONTAL or VERTICAL")
        else:
            if slope == VERTICAL:
                for e in self.forward_flippable_edges():
                    insert(x[self._norm(e)] <= 1)
                for e in self.backward_flippable_edges():
                    a, b, c, d = self.square_about_edge(e)
                    insert(x[self._norm(b)] + x[self._norm(c)] >= 1)
            elif slope == HORIZONTAL:
                for e in self.forward_flippable_edges():
                    a, b, c, d = self.square_about_edge(e)
                    insert(x[self._norm(b)] + x[self._norm(c)] >= 1)
                for e in self.backward_flippable_edges():
                    insert(x[self._norm(e)] <= 1)
            else:
                raise ValueError("slope must be HORIZONTAL or VERTICAL")

    def train_track_linear_space(self, slope=VERTICAL, backend=None):
        r"""
        Return the polytope determined by the switch equations (a linear subspace)

        INPUT:

        - ``slope`` - the slope for the train track (``HORIZONTAL`` or ``VERTICAL``)

        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB")
            sage: T.train_track_linear_space()
            Cone of dimension 2 in ambient dimension 3 made of 0 facets (backend=ppl)
            sage: T.train_track_linear_space(backend='sage')
            Cone of dimension 2 in ambient dimension 3 made of 0 facets (backend=sage)
        """
        from sage.rings.integer_ring import ZZ
        cs = ConstraintSystem()
        L = LinearExpressions(ZZ)
        ne = self.num_edges()
        x = [L.variable(e) for e in range(ne)]
        self._set_switch_conditions(cs.insert, x, slope)
        return cs.cone(backend)

    def train_track_polytope(self, slope=VERTICAL, low_bound=0, backend=None):
        r"""
        Return the polytope determined by the constraints.

        INPUT:

        - ``slope`` - the slope for the train track

        - ``low_bound`` - integer - optional lower bound for the lengths
          (default to 0)

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2),(-1,-2,-3)], [RED, RED, BLUE])
            sage: P = T.train_track_polytope(VERTICAL)
            sage: P
            Cone of dimension 2 in ambient dimension 3 made of 2 facets (backend=ppl)
            sage: sorted(P.rays())
            [[0, 1, 1], [1, 1, 0]]

            sage: P = T.train_track_polytope(VERTICAL, low_bound=3)  # not tested
            sage: P.generators()  # not tested
            Generator_System {ray(1, 1, 0), ray(0, 1, 1), point(3/1, 6/1, 3/1)}

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [GREEN, RED, BLUE])
            sage: sorted(T.train_track_polytope(VERTICAL).rays())
            [[0, 1, 1]]
            sage: sorted(T.train_track_polytope(HORIZONTAL).rays())
            [[1, 0, 1], [1, 1, 0]]

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [PURPLE, BLUE, RED])
            sage: sorted(T.train_track_polytope(VERTICAL).rays())
            [[1, 0, 1], [1, 1, 0]]
            sage: sorted(T.train_track_polytope(HORIZONTAL).rays())
            [[0, 1, 1]]

        One can also use other backends::

            sage: sorted(T.train_track_polytope(VERTICAL, backend='sage').rays())
            [[1, 0, 1], [1, 1, 0]]
            sage: sorted(T.train_track_polytope(HORIZONTAL, backend='sage').rays())
            [[0, 1, 1]]
        """
        from sage.rings.integer_ring import ZZ
        L = LinearExpressions(ZZ)
        cs = ConstraintSystem()
        ne = self.num_edges()
        variables = [L.variable(e) for e in range(ne)]
        self._set_train_track_constraints(cs.insert, variables, slope, low_bound, False)
        return cs.cone(backend)

    def train_track_min_solution(self, slope=VERTICAL, allow_degenerations=False):
        r"""
        Return the minimal integral point satisfying the constraints.

        INPUT:

        - ``slope`` - the slope of the train track

        - ``allow_degenerations`` - boolean - if ``True`` then allow certain
          degenerations to occur.

        OUTPUT: a point from ppl

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2),(-1,-2,-3)], [RED, RED, BLUE])
            sage: T.train_track_min_solution(VERTICAL)
            point(1/1, 2/1, 1/1)
            sage: T.train_track_min_solution(VERTICAL, allow_degenerations=True)
            point(0/1, 1/1, 1/1)

            sage: T.train_track_min_solution(HORIZONTAL)
            point(2/1, 1/1, 1/1)
            sage: T.train_track_min_solution(HORIZONTAL, allow_degenerations=True)
            point(1/1, 0/1, 1/1)
        """
        n = self.num_edges()
        M = ppl.MIP_Problem(n)

        x = [ppl.Variable(e) for e in range(n)]
        M.set_objective_function(-sum(x))
        self._set_train_track_constraints(M.add_constraint, x, slope, 1, allow_degenerations)
        return M.optimizing_point()

    def geometric_polytope(self, x_low_bound=0, y_low_bound=0, hw_bound=0, backend=None):
        r"""
        Return the geometric polytope of this veering triangulation.

        The geometric polytope is the polytope of length and heights data that
        corresponds to L-infinity Delaunay triangulations.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.geometric_polytope()
            Cone of dimension 4 in ambient dimension 6 made of 6 facets (backend=ppl)
            sage: T.geometric_polytope(x_low_bound=1, y_low_bound=1, hw_bound=1)  # not tested

            sage: T.geometric_polytope(backend='sage')
            Cone of dimension 4 in ambient dimension 6 made of 6 facets (backend=sage)
        """
        from sage.rings.integer_ring import ZZ
        L = LinearExpressions(ZZ)
        ne = self.num_edges()
        x = [L.variable(e) for e in range(ne)]
        y = [L.variable(ne+e) for e in range(ne)]
        cs = ConstraintSystem()
        self._set_train_track_constraints(cs.insert, x, VERTICAL, x_low_bound, False)
        self._set_train_track_constraints(cs.insert, y, HORIZONTAL, y_low_bound, False)
        self._set_geometric_constraints(cs.insert, x, y, hw_bound=hw_bound)
        return cs.cone(backend)

    def _complexify_generators(self, Gx):
        r"""
        Given ``Gx`` a matrix whose rows are admissible lengths return the
        corresponding admissible heights.

        EXAMPLES::

            sage: from veerer import *

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: Gx = matrix(QQ, 2, [s, t])
            sage: Gy = T._complexify_generators(Gx)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(0), VERTICAL)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(1), VERTICAL)
            sage: T._set_switch_conditions(T._tt_check, Gy.row(0), HORIZONTAL)
            sage: T._set_switch_conditions(T._tt_check, Gy.row(1), HORIZONTAL)

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(2, 3, 4, 5, 1, 2)
            sage: Gx = matrix(QQ, 2, [s, t])
            sage: Gy = T._complexify_generators(Gx)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(0), VERTICAL)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(1), VERTICAL)
            sage: T._set_switch_conditions(T._tt_check, Gy.row(0), HORIZONTAL)
            sage: T._set_switch_conditions(T._tt_check, Gy.row(1), HORIZONTAL)
        """
        ne = self.num_edges()
        ep = self._ep
        if Gx.ncols() != ne:
            raise ValueError
        Gy = Gx.__copy__()
        for j in range(ne):
            if ep[j] < j:
                raise ValueError("triangulation not in standard form")
            if self._colouring[j] == BLUE:
                for i in range(Gy.nrows()):
                    Gy[i, j] *= -1
        return Gy

    def _complexify_equations(self, Lx):
        r"""
        Given ``Lx`` a matrix whose right kernel is some subspace of
        admissible lengths, return the admissible heights.

        EXAMPLES::

            sage: from veerer import *

            sage: V = VectorSpace(QQ, 7)

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: Gx = matrix(QQ, 2, [s, t])
            sage: Gy = T._complexify_generators(Gx)
            sage: Lx = Gx.right_kernel_matrix()
            sage: V1 = V.subspace(T._complexify_equations(Lx))
            sage: V2 = V.subspace(Gy.right_kernel_matrix())
            sage: assert V1 == V2

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(2, 3, 4, 5, 1, 2)
            sage: Gx = matrix(QQ, 2, [s, t])
            sage: Gy = T._complexify_generators(Gx)
            sage: Lx = Gx.right_kernel_matrix()
            sage: V1 = V.subspace(T._complexify_equations(Lx))
            sage: V2 = V.subspace(Gy.right_kernel_matrix())
            sage: assert V1 == V2
        """
        ne = self.num_edges()
        ep = self._ep
        if Lx.ncols() != ne:
            raise ValueError
        Ly = Lx.__copy__()
        for j in range(ne):
            if ep[j] < j:
                raise ValueError("triangulation not in standard form")
            if self._colouring[j] == BLUE:
                for i in range(Ly.nrows()):
                    Ly[i, j] *= -1
        return Ly

    def _flat_structure_from_train_track_lengths(self, VH, VV, base_ring=None, mutable=False, check=True):
        r"""
        Return a flat structure from two vectors ``VH`` and ``VV``
        satisfying the train track equations.
        """
        from sage.modules.free_module import VectorSpace

        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ

        assert len(VH) == len(VV) == self.num_edges()
        assert all(x >=0 for x in VH)
        assert all(x >= 0 for x in VV)

        self._set_train_track_constraints(self._tt_check, VH, HORIZONTAL, False, False)
        self._set_train_track_constraints(self._tt_check, VV, VERTICAL, False, False)

        V = VectorSpace(base_ring, 2)
        vectors = [V((x, y if self._colouring[i] == RED else -y)) for \
                   i,(x,y) in enumerate(zip(VV, VH))]
        m = self.num_edges()
        n = self.num_half_edges()
        ep = self._ep
        vectors.extend(vectors[ep[e]] for e in range(m,n))

        # get correct signs for each triangle
        for i,j,k in self.faces():
            if det2(vectors[i], vectors[j]) < 0:
                vectors[j] = -vectors[j]
            if det2(vectors[j], vectors[k]) < 0:
                vectors[k] = -vectors[k]
            if vectors[i] + vectors[j] + vectors[k]:
                raise RuntimeError('bad vectors for %s:\n vec[%s] = %s\n vec[%s] = %s\n vec[%s] = %s' \
                                   % (self.to_string(), self._edge_rep(i), vectors[i], self._edge_rep(j), \
                                      vectors[j], self._edge_rep(k), vectors[k]))

            if det2(vectors[k], vectors[i]) < 0:
                raise RuntimeError

        from .flat_structure import FlatVeeringTriangulation
        return FlatVeeringTriangulation(self, vectors, mutable=mutable, check=check)

    def flat_structure_middle(self, backend=None):
        r"""
        Return a flat structure with this Veering triangulation.

        Note that this triangulation must be core. The point is chosen
        by taking the interior point of the polytope obtained by
        summing each ray.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)", "RRB")
            sage: T.flat_structure_middle()
            FlatVeeringTriangulation(Triangulation("(0,1,2)"), [(1, 2), (-2, -1), (1, -1)])

            sage: from surface_dynamics import *              # optional - surface_dynamics
            sage: Q = QuadraticStratum({1:4, -1:4})           # optional - surface_dynamics
            sage: CT = VeeringTriangulation.from_stratum(Q)   # optional - surface_dynamics
            sage: CT.flat_structure_middle()                  # optional - surface_dynamics
            FlatVeeringTriangulation(Triangulation("(0,18,~17)(1,20,~19)...(19,~18,~0)"), [(3, 3), (1, 1), ..., (-3, -3)])

        TESTS::

            sage: from veerer import *
            sage: t = VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBPBRBPRB")
        """
        n = self.num_edges()

        PH = self.train_track_polytope(HORIZONTAL, backend=backend)
        PV = self.train_track_polytope(VERTICAL, backend=backend)

        # pick sum of rays
        VH = PH.rays()
        VH = [sum(v[i] for v in VH) for i in range(n)]
        VV = PV.rays()
        VV = [sum(v[i] for v in VV) for i in range(n)]

        return self._flat_structure_from_train_track_lengths(VH, VV)

    def flat_structure_min(self, allow_degenerations=False):
        r"""
        Return a flat structure with this Veering triangulation.

        Note that this triangulation must be core. The point is chosen
        by taking the minimum integral point in the cone.

        EXAMPLES::

            sage: from veerer import *

            sage: from surface_dynamics import *             # optional - surface_dynamics
            sage: Q = QuadraticStratum({1:4, -1:4})          # optional - surface_dynamics
            sage: CT = VeeringTriangulation.from_stratum(Q)  # optional - surface_dynamics
            sage: CT.flat_structure_min()                    # optional - surface_dynamics
            FlatVeeringTriangulation(Triangulation("(0,18,~17)(1,20,~19)...(19,~18,~0)"), [(3, 3), (1, 1), ..., (-3, -3)])

        By allowing degenerations you can get a simpler solution but
        with some of the edges horizontal or vertical::

            sage: F = CT.flat_structure_min(True)                 # optional - surface_dynamics
            sage: F                                               # optional - surface_dynamics
            FlatVeeringTriangulation(Triangulation("(0,18,~17)(1,20,~19)...(19,~18,~0)"), [(3, 3), (1, 1), ..., (-3, -3)])
            sage: F.to_veering_triangulation()                    # optional - surface_dynamics
            VeeringTriangulation("(0,18,~17)(1,20,~19)...(19,~18,~0)", "RRRRRRRRBBBBBBBBBGBBBBBP")
        """
        VH = self.train_track_min_solution(HORIZONTAL, allow_degenerations=allow_degenerations)
        VV = self.train_track_min_solution(VERTICAL, allow_degenerations=allow_degenerations)

        assert VH.is_point()
        assert VV.is_point()

        from sage.rings.rational import Rational
        assert VH.divisor() == 1
        VH = [Rational(c) for c in VH.coefficients()]

        assert VV.divisor() == 1
        VV = [Rational(c) for c in VV.coefficients()]

        return self._flat_structure_from_train_track_lengths(VH, VV)

    def flat_structure_geometric_middle(self, backend=None):
        r"""
        Return a geometric flat structure obtained by averaging the
        vertices of the geometric polytope.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.flat_structure_geometric_middle()
            FlatVeeringTriangulation(Triangulation("(0,1,2)(~2,~0,~1)"), [(4, 9), (-9, -4), (5, -5), (-5, 5), (9, 4), (-4, -9)])
        """
        ne = self.num_edges()
        r = self.geometric_polytope(backend=backend).rays()
        VV = [sum(v[i] for v in r) for i in range(ne)]
        VH = [sum(v[ne + i] for v in r) for i in range(ne)]

        return self._flat_structure_from_train_track_lengths(VH, VV)

    def zippered_rectangles(self, x, y, base_ring=None, check=True):
        r"""
        Return the zippered rectangle surface associated to the holonomy data ``x`` and ``y``.

        This construction only works for Abelian differentials.

        For each wedge with an upward separatrix we consider a horizontal
        segment based at the bottom of the wedge as wide as the triangle. Then
        we consider downward vertical separatrices and extend them until they
        touch the union of these horizontal segments.

        INPUT:

        - ``x``, ``y`` -- list of positive numbers (as many as edges in this veering
          triangulation)

        - ``base_ring`` -- an optional base ring to build the surface

        - ``check`` -- optional boolean (default ``True``)

        EXAMPLES::

            sage: from veerer import VeeringTriangulation, BLUE, RED

            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RBB")
            sage: x = [1, 2, 1]
            sage: y = [1, 1, 2]
            sage: vt.zippered_rectangles(x, y)  # optional: sage_flatsurf
            Translation Surface built from a square and a rectangle

            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,4)(~2,5,3)(~3,~4,~5)", "RBBRBR")
            sage: x = [1,2,1,2,1,1]
            sage: y = [1,1,2,1,2,3]
            sage: vt.zippered_rectangles(x, y)  # optional: sage_flatsurf
            Translation Surface built from 3 squares and a rectangle

            sage: vt = VeeringTriangulation("(0,~2,1)(2,~8,~3)(3,~7,~4)(4,6,~5)(5,8,~6)(7,~1,~0)", "PRBPRBPBR")
            sage: R0, R1 = vt.dehn_twists(RED)
            sage: B0, B1 = vt.dehn_twists(BLUE)
            sage: f = B0 **2 * R0 **3 * B1 * R1 ** 5
            sage: a, x, y = f.self_similar_widths_and_heights()
            sage: S = vt.zippered_rectangles(x, y)  # optional: sage_flatsurf
            sage: S  # optional: sage_flatsurf
            Translation Surface built from 6 rectangles

        We now check that labelling of the rectangles in ``S`` coincide with
        the order of faces in the veering triangulation::

            sage: from veerer import HORIZONTAL
            sage: for i, r in enumerate(vt.right_wedges(HORIZONTAL)):  # optional: sage_flatsurf
            ....:     e0, e1, e2, e3 = S.polygon(i).erase_marked_vertices().edges()
            ....:     w = e0[0]
            ....:     h = e1[1]
            ....:     l = vt.previous_in_face(r)
            ....:     nr = vt._norm(r)
            ....:     nl = vt._norm(l)
            ....:     assert (w == x[nr] and h == y[nl]) or (w == x[nl] and h == y[nr])

        TESTS::

            sage: from veerer import VeeringTriangulation, BLUE, RED
            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RBB")
            sage: x = [1, 2, 1]
            sage: y = [1, 1, 2]
            sage: vt.zippered_rectangles(x, y, base_ring=AA)  # optional: sage_flatsurf
            Translation Surface built from a square and a rectangle
        """
        ans, edge_orientations = self.is_abelian(certificate=True)
        if not ans:
            raise ValueError('the construction is only valid for Abelian differentials')

        if check:
            x, y = self._check_xy(x, y)
            if base_ring is None:
                from sage.structure.sequence import Sequence
                base_ring = Sequence(list(x) + list(y)).universe()
        elif base_ring is None:
            base_ring = self.base_ring()

        colouring = self.colouring_from_xy(x, y, check=False)
        for col0, col1 in zip(colouring, self._colouring):
            if col1 == GREEN or col1 == PURPLE:
                continue
            if col0 != col1:
                raise ValueError('invalid monodromy data')

        ep = self._ep

        # Each separatrix is in between two consecutive edges of
        # distinct colours around a vertex
        # We label separatrices as follows
        #
        #            4i+2
        #             ^
        #             |
        #             |
        #           c2|c1
        # 4i+3  <-----x----> 4i+1
        #             |c0
        #             |
        #             |
        #             v
        #             4i

        # compute for each half-edge the label of the next separatrix
        left_wedges = []
        right_wedges = []
        next_separatrix = [-1] * self._n
        previous_separatrix = [-1] * self._n
        i = 0
        for vertex in self.vertices():
            # find a blue half-edge after a down separatrix, that is
            # consecutive a, b at a vertex such that
            # edge_orientations[a] == False and edge_orientations[b] == True
            a = vertex[0]
            b = vertex[1]
            while not edge_orientations[a] or edge_orientations[b]:
                a = b
                b = self.next_at_vertex(a)
            assert colouring[a] == RED and colouring[b] == BLUE, (a, colouring[a], colouring[b])
            while previous_separatrix[b] == -1:
                if i % 2 == 0:
                    assert colouring[a] == RED and colouring[b] == BLUE, (i, a, colouring[a], b, colouring[b])
                else:
                    assert colouring[a] == BLUE and colouring[b] == RED, (i, a, colouring[a], b, colouring[b])

                right_wedges.append(a)
                left_wedges.append(b)
                a = b
                b = self.next_at_vertex(a)
                previous_separatrix[a] = i
                next_separatrix[a] = i + 1
                while colouring[a] == colouring[b]:
                    previous_separatrix[b] = i
                    next_separatrix[b] = i + 1
                    a = b
                    b = self.next_at_vertex(a)
                # colour change
                i += 1

            # correct the last quadrant
            assert colouring[a] == RED and colouring[b] == BLUE
            j = previous_separatrix[b]
            while colouring[a] == RED:
                next_separatrix[a] = j
                a = self.previous_at_vertex(a)

            # check that we did a multiple of 2pi
            assert i % 4 == 0

        assert all(x != -1 for x in next_separatrix)
        assert all(x != -1 for x in previous_separatrix)
        assert len(right_wedges) == len(left_wedges) == 2 * self.num_faces()
        assert all(colouring[right_wedges[i]] == RED for i in range(0, 2 * self.num_faces(), 2))
        assert all(colouring[right_wedges[i]] == BLUE for i in range(1, 2 * self.num_faces(), 2))
        assert all(colouring[left_wedges[i]] == BLUE for i in range(0, 2 * self.num_faces(), 2))
        assert all(colouring[left_wedges[i]] == RED for i in range(1, 2 * self.num_faces(), 2))

        # In order to have a consistent labelling between the triangles as provided by self.faces()
        # and the rectangles we compute the face index associated to each half-edge and use it
        # later to order the rectangles
        half_edge_face_index = [-1] * self._n
        for i, face in enumerate(self.faces()):
            for e in face:
                half_edge_face_index[e] = i

        # There are as many rectangles in the Markov partitions as triangles in
        # the veering triangulation
        # Each left/right/bottom separatrix has a certain number of cut points
        # (on both sides)
        rectangles = [None] * self.num_faces()
        for i in range(1, 2 * self.num_faces(), 2):
            l = left_wedges[i]
            L = ep[l]
            r = right_wedges[i]
            R = ep[r]
            e = self.next_in_face(r)
            E = ep[e]
            assert self.next_in_face(e) == L

            # The rectangles are always built starting from the bottom left corner
            # as in the following picture
            #
            #   o----------------o
            #   |  p5         p4 |
            #   |p6            p3|
            #   |                |
            #   |p7            p2|
            #   | p0          p1 |
            #   o----------------o
            if i % 4 == 1:
                # right rectangle
                # bottom side
                p0 = (next_separatrix[R], RIGHT, x[r])
                if x[l] < x[r]:
                    # small case: x[l] = x[r] - x[e]
                    p1 = (next_separatrix[R], RIGHT, x[e])
                else:
                    # big case: x[l] = x[r] + x[e]
                    p1 = (previous_separatrix[e], LEFT, x[e])
                # right side
                p2 = (next_separatrix[L], RIGHT, y[e])
                p3 = (next_separatrix[L], RIGHT, y[l])
                # top side
                p4 = (next_separatrix[r], RIGHT, x[l])
                p5 = (next_separatrix[r], RIGHT, 0)
                # left side
                p6 = (previous_separatrix[r], LEFT, 0)
                p7 = (previous_separatrix[r], LEFT, y[r])

            else:
                # bottom side
                if x[r] < x[l]:
                    # small case: x[r] = x[l] - x[e]
                    p0 = (previous_separatrix[L], LEFT, x[e])
                else:
                    # big case: x[r] = x[l] + x[e]
                    p0 = (next_separatrix[E], RIGHT, x[e])
                p1 = (previous_separatrix[L], LEFT, x[l])
                # right side
                p2 = (next_separatrix[l], RIGHT, y[l])
                p3 = (next_separatrix[l], RIGHT, 0)
                # top side
                p4 = (next_separatrix[r], LEFT, 0)
                p5 = (next_separatrix[r], LEFT, x[r])
                # left side
                p6 = (previous_separatrix[R], LEFT, y[r])
                p7 = (previous_separatrix[R], LEFT, y[e])

            j = half_edge_face_index[r]
            rectangles[j] = (p0, p1, p2, p3, p4, p5, p6, p7)

        from .tatami_decomposition import tatami_decomposition
        return tatami_decomposition(rectangles, base_ring)

    def is_core(self, backend=None):
        r"""
        Test whether this coloured triangulation is core.

        EXAMPLES::

            sage: from veerer import *

            sage: triangles = [(-24, -2, -23), (-22, 2, 22), (-21, 3, 21), (-20, 4, 20),
            ....:              (-19, 1, 19), (-18, 6, 18), (-17, 7, 17), (-16, 16, -1),
            ....:              (-15, -8, -14), (-13, 13, -7), (-12, 12, -6), (-11, 11, -5),
            ....:              (-10, -4, 8), (-9, -3, 23), (0, 15, 14), (5, 10, 9)]
            sage: colours = [RED, RED, RED, RED, RED, RED, RED, RED, BLUE, BLUE, BLUE,
            ....:            BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE,
            ....:            BLUE, BLUE, BLUE]

            sage: T = VeeringTriangulation(triangles, colours)
            sage: T.is_core()
            True

            sage: U = T.copy(mutable=True)
            sage: U.flip(10, BLUE)
            sage: U.is_core()
            True

            sage: U = T.copy(mutable=True)
            sage: U.flip(10, RED)
            sage: U.is_core()
            False

        Examples involving linear subspaces::

            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: vt.as_linear_family().is_core()
            True
            sage: VeeringTriangulationLinearFamily(vt, [1, 0, -1]).is_core()
            False
        """
        if any(c == PURPLE or c == GREEN for c in self._colouring):
            raise ValueError('core not implemented with PURPLE or GREEN colour')

        # TODO: could use v.has_curve for every edge?
        # In theory LP should be much faster but in practice (in small dimensions)
        # polytope is much better
        d = self.dimension()
        return self.train_track_polytope(HORIZONTAL, backend=backend).affine_dimension() == d and \
               self.train_track_polytope(VERTICAL, backend=backend).affine_dimension() == d

    def is_geometric(self, backend=None):
        r"""
        Test whether this coloured triangulation is geometric.

        EXAMPLES::

            sage: from veerer import VeeringTriangulation

            sage: vt1 = VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBBBRBBRB")
            sage: vt2 = VeeringTriangulation("(0,~8,~3)(1,6,~2)(2,~1,~0)(3,7,~4)(4,8,~5)(5,~7,~6)", "RBBBRBRBB")
            sage: vt1.is_geometric()
            True
            sage: vt2.is_geometric()
            False

        An example in genus 2 involving a linear subspace::

            sage: from veerer import VeeringTriangulations, VeeringTriangulationLinearFamily
            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: f = VeeringTriangulationLinearFamily(T, [s, t])
            sage: f.is_geometric()
            True
        """
        dim = self.dimension()
        P = self.geometric_polytope(backend=backend)
        Pdim = P.affine_dimension()
        assert Pdim <= 2 * dim
        return Pdim == 2 * dim

    def balanced_polytope(self, slope=VERTICAL, homogeneous=False, backend=None):
        r"""
        Return the set of balanced coordinates for this veering triangulation

        INPUT:

        - ``slope`` - either ``VERTICAL`` (default) or ``HORIZONTAL``

        - ``homogeneous`` - boolean. Default to ``False``.
        """
        from sage.rings.rational_field import QQ
        ne = self.num_edges()
        cs = ConstraintSystem()
        L = LinearExpressions(QQ)
        x = [L.variable(e) for e in range(ne)]
        self._set_train_track_constraints(cs.insert, x, slope, False, False)
        self._set_balance_constraints(cs.insert, x, slope, homogeneous)
        return cs.polyhedron(backend)

    def is_balanced(self, slope=VERTICAL, backend=None):
        r"""
        Check balanceness

        EXAMPLES::

            sage: from veerer import VeeringTriangulation

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RBR")
            sage: T.is_balanced()  # not tested
            True

            sage: T = VeeringTriangulation("(0,~3,2)(1,4,~2)(3,5,~4)(~5,~1,~0)", "RBBBRB")
            sage: T.is_balanced()  # not tested
            False

            sage: T = VeeringTriangulation("(0,1,8)(2,~7,~1)(3,~0,~2)(4,~5,~3)(5,6,~4)(7,~8,~6)", "BRRRRBRBR")
            sage: T.is_balanced()  # not tested
            False
            sage: T.rotate()
            sage: T.is_balanced()  # not tested
            True
        """
        return self.balanced_polytope(backend=backend).affine_dimension() == self.dimension()

    def edge_has_curve(self, e, verbose=False):
        r"""
        INPUT:

        - e - edge label

        OUTPUT: boolean - whether there is a curve which crosses the associated
        (dual) train-track in the correct direction.

        EXAMPLES::

            sage: from veerer import *

            sage: t = [(-6, -4, -5), (-3, -1, 3), (-2, 0, 4), (1, 5, 2)]
            sage: c = [RED, RED, BLUE, RED, BLUE, RED]
            sage: T0 = VeeringTriangulation(t, c)
            sage: T0.is_core()
            True

        Flipping edge 3 in RED is fine (it remains a core triangulation)::

            sage: T1 = T0.copy(mutable=True)
            sage: T1.flip(3, RED)
            sage: T1.edge_has_curve(3)
            True
            sage: T1.is_core()
            True

        However, flipping edge 3 in BLUE leads to a non-core triangulation::

            sage: T2 = T0.copy(mutable=True)
            sage: T2.flip(3, BLUE)
            sage: T2.edge_has_curve(3)
            False
            sage: T2.is_core()
            False

        Equivantly, the train track polytope is degenerate::

            sage: P1 = T1.train_track_polytope(VERTICAL)
            sage: P1.affine_dimension()
            3
            sage: P2 = T2.train_track_polytope(VERTICAL)
            sage: P2.affine_dimension()
            2
        """
        # TODO: we should only searching for vertex cylces; i.e. not allow more
        # than two pairs (i, ~i) to be both seen (barbell are fine but not more)

        # TODO: we want the._colouring to also take care of negative edges
        # (bad alternative: take "norm" function from flipper)
        colouring = self._colouring
        edge_rep = self._edge_rep
        ep = self._ep

        if verbose:
            print('[edge_has_curve] checking edge %s with colour %s' % (edge_rep(e), colouring[e]))

        a, b, c, d = self.square_about_edge(e)
        if colouring[a] == BLUE or colouring[b] == RED:
            assert colouring[c] == BLUE or colouring[d] == RED
            POS, NEG = BLUE, RED
            if verbose:
                print('[edge_has_curve] checking HORIZONTAL track')
        else:
            assert colouring[a] == RED or colouring[b] == BLUE
            assert colouring[c] == RED or colouring[d] == BLUE
            POS, NEG = RED, BLUE
            if verbose:
                print('[edge_has_curve] checking VERTICAL track')

        # check alternating condition
        assert colouring[e] == BLUE or colouring[e] == RED
        assert colouring[a] != colouring[b], (a, b, colouring[a], colouring[b])
        assert colouring[c] != colouring[d], (c, d, colouring[c], colouring[d])

        if colouring[e] == NEG:
            start = b
            end = ep[d]
            if verbose:
                print('[edge_has_curve] try to find path from b=%s to ~d=%s' %
                         (edge_rep(start), edge_rep(end)))
        else:
            start = a
            end = ep[c]
            if verbose:
                print('[edge_has_curve] try to find path from a=%s to ~c=%s' %
                         (edge_rep(start), edge_rep(end)))

        if start == end:
            return True

        n = self._n
        fp = self._fp
        seen = [False] * n
        seen[start] = True
        q = collections.deque()
        q.append(start)
        while q:
            f = q.popleft()
            if verbose:
                print('[edge_has_curve] crossing %s' % edge_rep(f))

            # here we set r and s so that the triangle is (r, s, ~f)
            r = fp[ep[f]]
            s = fp[r]
            if verbose:
                print('[edge_has_curve] switch with r=%s s=%s' % (edge_rep(r), edge_rep(s)))
            if not seen[r] and not (colouring[r] == POS and colouring[ep[f]] == NEG):
                if r == end:
                    if verbose:
                        print('[edge_has_curve] done at %s' % edge_rep(r))
                    return True
                seen[r] = True
                q.append(r)
                if verbose:
                    print('[edge_has_curve] adding %s on top of the queue' % edge_rep(r))
            if not seen[s] and not (colouring[s] == NEG and colouring[ep[f]] == POS):
                if s == end:
                    if verbose:
                        print('[edge_has_curve] done at %s' % edge_rep(s))
                    return True
                seen[s] = True
                q.append(s)
                if verbose:
                    print('[edge_has_curve] adding %s on top of the queue' % edge_rep(s))

        return False

    def geometric_flips(self, backend=None):
        r"""
        Return the list of geometric flips.

        A flip, a rather a list of flips, is geometric if it arises generically
        as a flip of L^oo-Delaunay triangulations along Teichmueller geodesics.
        Each geometric flip corresponds to a facet of the geometric polytope.

        OUTPUT: a list of pairs ``(edge_number, new_colour)``

        EXAMPLES::

            sage: from veerer import *

            sage: vt = VeeringTriangulation("(0,2,3)(1,4,~0)(5,6,~1)", "BRRBBBB")
            sage: sorted(vt.geometric_flips())
            [([3], 1), ([3], 2), ([4], 1), ([4], 2), ([5], 1), ([5], 2)]
            sage: sorted(vt.geometric_flips(backend='sage'))
            [([3], 1), ([3], 2), ([4], 1), ([4], 2), ([5], 1), ([5], 2)]

        L-shaped square tiled surface with 3 squares (given as a sphere with
        3 triangles). It has two geometric neighbors corresponding to simultaneous
        flipping of the diagonals 3, 4 and 5::

            sage: from veerer import *
            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: f = VeeringTriangulationLinearFamily(T, [s, t])
            sage: sorted(f.geometric_flips(backend='ppl'))
            [([3, 4, 5], 1), ([3, 4, 5], 2)]
            sage: sorted(f.geometric_flips(backend='sage'))
            [([3, 4, 5], 1), ([3, 4, 5], 2)]
            sage: sorted(f.geometric_flips(backend='normaliz-QQ'))  # optional - pynormaliz
            [([3, 4, 5], 1), ([3, 4, 5], 2)]

        To be compared with the geometric flips in the ambient stratum::

            sage: sorted(T.geometric_flips())
            [([3], 1), ([3], 2), ([4], 1), ([4], 2), ([5], 1), ([5], 2)]
            sage: sorted(T.as_linear_family().geometric_flips())
            [([3], 1), ([3], 2), ([4], 1), ([4], 2), ([5], 1), ([5], 2)]


        A more complicated example in which edge 4 have a forced colour after
        flip and where the flippable edges 0 and 3 are not part of any geometric
        flips::

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(2, 3, 5, 2, 1, 1)
            sage: f = VeeringTriangulationLinearFamily(T, [s, t])
            sage: T.flippable_edges()
            [0, 3, 4, 5, 6]
            sage: sorted(f.geometric_flips(backend='ppl'))
            [([4], 2), ([5], 1), ([5], 2)]
            sage: sorted(f.geometric_flips(backend='sage'))
            [([4], 2), ([5], 1), ([5], 2)]
            sage: sorted(f.geometric_flips(backend='normaliz-QQ'))  # optional - pynormaliz
            [([4], 2), ([5], 1), ([5], 2)]

        TESTS::

            sage: from veerer import VeeringTriangulation
            sage: fp = "(0,~8,~7)(1,3,~2)(2,7,~3)(4,6,~5)(5,8,~6)(~4,~1,~0)"
            sage: cols = "RBRRRRBBR"
            sage: vt = VeeringTriangulation(fp, cols)
            sage: sorted(vt.geometric_flips())
            [([2], 1), ([2], 2), ([4, 8], 1), ([4, 8], 2)]
            sage: sorted(vt.as_linear_family().geometric_flips())
            [([2], 1), ([2], 2), ([4, 8], 1), ([4, 8], 2)]
        """
        from sage.matrix.constructor import matrix

        dim = self.dimension()
        base_ring = self.base_ring()
        ne = ambient_dim = self.num_edges()
        L = LinearExpressions(base_ring)
        x = [L.variable(e) for e in range(ne)]
        y = [L.variable(ne + e) for e in range(ne)]
        P = self.geometric_polytope(backend=backend)
        if P.affine_dimension() != 2 * dim:
            raise ValueError('not geometric P.dimension() = {} while 2 * dim = {}'.format(P.affine_dimension(), 2 * dim))

        # compute the Delaunay flip facets and the associated
        # subset of edges
        eqns = matrix(base_ring, P.eqns())
        delaunay_facets = {}
        for e in self.forward_flippable_edges():
            a, b, c, d = self.square_about_edge(e)
            constraint = x[self._norm(e)] == y[self._norm(a)] + y[self._norm(d)]
            constraint = constraint.coefficients(dim=2*ne, homogeneous=True)
            linear_form_project(eqns, constraint)
            linear_form_normalize(base_ring, constraint)
            constraint = tuple(constraint)
            if constraint in delaunay_facets:
                delaunay_facets[constraint].append(e)
            else:
                delaunay_facets[constraint] = [e]

        # determine the possible colours
        neighbours = []
        for ieq, edges in delaunay_facets.items():
            # build the facet
            F = P.add_constraint(L(ieq) == 0)
            if not F.affine_dimension() == 2 * dim - 1:
                continue

            # test each edge colour conditions
            # NOTE: all simultaneous flips must be of the same colour
            assert all(self._colouring[e] == self._colouring[edges[0]] for e in edges)
            # NOTE: the equations for the different edges are all equivalent, it
            # is hence enough to use the first edge
            a, b, c, d = self.square_about_edge(edges[0])
            Fred = F.add_constraint(x[self._norm(a)] <= x[self._norm(d)])
            if Fred.affine_dimension() == 2 * dim - 1:
                neighbours.append((edges, RED))
                Fblue = F.add_constraint(x[self._norm(a)] >= x[self._norm(d)])
                if Fblue.affine_dimension() == 2 * dim - 1:
                    neighbours.append((edges, BLUE))
            else:
                neighbours.append((edges, BLUE))

        return neighbours

    def random_forward_flip_sequence(self, length=1, relabel=False):
        V = self.copy()
        cols = [RED, BLUE]
        flips = []
        for _ in range(length):
            e = choice(V.forward_flippable_edges())
            col = choice(cols)

            # TODO: this is a bit annoying. There should be a method to
            # test what are the valid colouring
            V.flip(e, col, reduced=False)
            if not V.edge_has_curve(e):
                col = BLUE if col == RED else RED
            V.flip_back(e, PURPLE)
            V.flip(e, col)
            flips.append((e, col))

        if relabel:
            relabelling = perm_random_centralizer(self._ep)
        else:
            relabelling = perm_id(self._n)

        from .flip_sequence import VeeringFlipSequence
        return VeeringFlipSequence(self, flips, relabelling)

    # TODO: this will not work with purple edges
    def random_forward_flip(self, repeat=1):
        r"""
        Apply a forward flip randomly among the ones that keeps the triangulation core.

        INPUT:

        - ``repeat`` - integer (default 1) - if provided make ``repeat`` flips instead of 1.
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation; use a mutable copy instead')

        cols = [RED, BLUE]
        for _ in range(repeat):
            e = choice(self.forward_flippable_edges())
            old_col = self._colouring[e]
            shuffle(cols)
            for c in cols:
                self.flip(e, c)
                if self.edge_has_curve(e):
                    break
                else:
                    self.flip_back(e, old_col)


class VeeringTriangulations:
    @staticmethod
    def L_shaped_surface(a1, a2, b1, b2, t1=0, t2=0):
        r"""
        Return the quotient of the L-shaped surface.

        The corresponding surface belongs to the quadratic stratum Q(1, -1^5) and
        its orientation double cover in H(2).

        EXAMPLES::

            sage: from veerer import *
            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: T
            VeeringTriangulation("(0,2,3)(1,4,~0)(5,6,~1)", "BRRBBBB")
            sage: s
            (0, 1, 1, 1, 1, 1, 0)
            sage: t
            (1, 0, 0, 1, 1, 1, 1)
            sage: T._set_switch_conditions(T._tt_check, s, VERTICAL)
            sage: T._set_switch_conditions(T._tt_check, t, VERTICAL)

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(2, 3, 4, 5, 1, 2)
            sage: T._set_switch_conditions(T._tt_check, s, VERTICAL)
            sage: T._set_switch_conditions(T._tt_check, t, VERTICAL)
        """
        # Return the (quotient by the hyperelliptic involution of the) L-shaped surface
        # together with the equations of the GL2R deformation
        #

        #       +----+
        #    ^  |    |
        #    |   |    |
        # b2 |    |    |
        #    |     |    |
        #    v      |    |
        #            +    +---------+
        #    ^       |              |
        #    |       |              |
        # b1 |       |              |
        #    |        |              |
        #    v        |              |
        #             +----+---------+
        #             <---> <------->
        #               a1      a2
        #
        #     <----><->
        #       t2  t1

        if not isinstance(a1, numbers.Integral) or \
           not isinstance(a2, numbers.Integral) or \
           not isinstance(b1, numbers.Integral) or \
           not isinstance(b2, numbers.Integral) or \
           not isinstance(t1, numbers.Integral) or \
           not isinstance(t2, numbers.Integral):
            raise TypeError("a1, a2, b1, b2, t1, t2 must be integers")
        a1 = int(a1)
        a2 = int(a2)
        b1 = int(b1)
        b2 = int(b2)
        if a1 <= 0 or a2 <= 0 or b1 <= 0 or b2 <= 0 or t1 < 0 or t2 < 0:
            raise ValueError("a1, a2, b1, b2 must be positive and t1, t2 non-negative")

        T = VeeringTriangulation("(0,2,3)(~0,1,4)(~1,5,6)", [BLUE, RED, RED, BLUE, BLUE, BLUE, BLUE])
        s = (t1, a1, a2, a2 + t1, a1 + t1, a1 + t2, t2)
        t = (b1, 0, 0, b1, b1, b2, b2)

        return T, s, t

    @staticmethod
    def ngon(n):
        #TODO: here we should include the Teichmueller curve constraints
        # (over the real part of cyclotomic fields)
        n = int(n)
        assert(n > 4 and n % 2 == 0)
        m = (n - 2) // 2
        T = [(i, 2*m+i, ~(i+1)) for i in range(m)] + \
            [(~0, ~(2*m), ~(m+1))] + \
            [(m+i, ~(2*m+i), ~(m+i+1)) for i in range(1, m-1)] + \
            [(2*m-1, ~(3*m-1), m)]
        colouring = [RED] * (2*m) + [BLUE] * m
        return VeeringTriangulation(T, colouring)
