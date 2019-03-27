r"""
Veering triangulations of surfaces.
"""
from __future__ import print_function, absolute_import

import numbers

from collections import deque
from itertools import product

from .constants import *
from .permutation import *
from .misc import det2
from .triangulation import Triangulation

from .env import curver, sage, surface_dynamics, ppl, flipper, random, require_package

from random import choice, shuffle

def ppl_cone_to_hashable(P):
    r"""
    EXAMPLES::

        sage: from veerer.veering_triangulation import VeeringTriangulation, ppl_cone_to_hashable, ppl_cone_from_hashable
        sage: from surface_dynamics import *

        sage: T = VeeringTriangulation.from_stratum(AbelianStratum(2))
        sage: P = T.geometric_polytope()
        sage: h = ppl_cone_to_hashable(P)
        sage: P == ppl_cone_from_hashable(h)
        True
    """
    eqns = []
    ieqs = []
    for constraint in P.minimized_constraints():
        if constraint.inhomogeneous_term():
            raise ValueError('not a cone')
        if constraint.is_equality():
            eqns.append(tuple(constraint.coefficients()))
        elif constraint.is_inequality():
            ieqs.append(tuple(constraint.coefficients()))
        else:
            raise RuntimeError
    eqns.sort()
    ieqs.sort()
    return (P.space_dimension(), tuple(eqns), tuple(ieqs))

def ppl_cone_from_hashable(args):
    d, eqns, ieqs = args
    P = ppl.C_Polyhedron(d)
    for constraint in eqns:
        P.add_constraint(sum(coeff * ppl.Variable(i) for i,coeff in enumerate(constraint)) == 0)
    for constraint in ieqs:
        P.add_constraint(sum(coeff * ppl.Variable(i) for i,coeff in enumerate(constraint)) >= 0)
    return P

class VeeringTriangulation(Triangulation):
    r"""
    Veering triangulation.

    A *veering triangulation* is a triangulation of a surface together with
    a colouring of the edges in red or blue so that there is no monochromatic
    face and no monochromatic vertex.

    EXAMPLES::

        sage: from veerer import *

    Built from an explicit triangulation (in cycle or list form) and a list of colors::

        sage: VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
        VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB")

    From a stratum::

        sage: from surface_dynamics import *

        sage: VeeringTriangulation.from_stratum(AbelianStratum(2))
        VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")

        sage: VeeringTriangulation.from_stratum(QuadraticStratum({1:4}))
        VeeringTriangulation("(0,12,~11)(1,13,~12)...(16,~15,~1)", "RRRRRRBBBBBBBBBBBB")

    From a flipper pseudo-Anosov map::

        sage: import flipper

        sage: T = flipper.load('S_2_1')
        sage: h = T.mapping_class('abcD')
        sage: h.is_pseudo_anosov()
        True
        sage: VeeringTriangulation.from_pseudo_anosov(h)
        VeeringTriangulation("(0,~3,~1)...(12,~14,~10)(~2,~9,~4)", "RBRBRRBRBBBBRBR")
    """
    __slots__ = ['_colouring']

    def __init__(self, triangulation,  colouring, check=True):
        Triangulation.__init__(self, triangulation, check=False)

        if isinstance(colouring, str):
            colouring = [RED if c == 'R' else BLUE for c in colouring]

        n = self._n
        # A list : edge_indices --> {Red, Blue}
        if len(colouring) == self.num_edges():
            colouring = [colouring[self._norm(i)] for i in range(n)]
        elif len(colouring) == n:
            colouring = colouring
        else:
            raise ValueError('wrong colouring argument')

        self._colouring = array('l', colouring)

        if check:
            self._check(ValueError)

    def _check(self, error=RuntimeError):
        Triangulation._check(self, error)
        n = self.num_half_edges()
        ep = self._ep
        ev = self._vp
        cols = self._colouring
        if not isinstance(cols, array) or \
           len(cols) != n or \
           any(col not in COLOURS for col in cols) or \
           any(cols[e] != cols[ep[e]] for e in range(n)):
            raise error('bad colouring attribute')

        # no monochromatic face
        allowed = [set([BLUE,RED]),
                   set([PURPLE,GREEN,BLUE]),
                   set([PURPLE,GREEN,RED]),
                   set([PURPLE,BLUE,RED]),
                   set([GREEN,BLUE,RED])]
        for f in self.faces():
            fcols = set(cols[e] for e in f)
            if not (PURPLE in fcols and GREEN in fcols) and \
               not (RED in fcols and BLUE in fcols):
                raise error('monochromatic face {}'.format(f))

        # no purple -> {blue,purple} or green -> {red,green} around a vertex
        for e in range(n):
            c0 = cols[e]
            c1 = cols[ev[e]]
            if c0 == PURPLE and c1 == PURPLE:
                raise error('two consecutive purple edges in a vertex at {}'.format(e))
            if c0 == GREEN and c1 == GREEN:
                raise error('two consecutive green edges in a vertex at {}'.format(e))
            if c0 == PURPLE and c1 == BLUE:
                raise error('blue after purple in a vertex at {}'.format(e))
            if c0 == GREEN and c1 == RED:
                raise error('red after green in a vertex at {}'.format(e))

        # no monochromatic vertex
        for v in self.vertices():
            if len(set(self._colouring[e] for e in v)) == 1:
                raise error('monochromatic vertex {}'.format(v))

    @classmethod
    def from_pseudo_anosov(cls, h):
        r"""
        Construct the coloured triangulation of a pseudo-Anosov homeomorphism.

        EXAMPLES::

            sage: from flipper import *
            sage: from veerer import *

            sage: T = flipper.load('SB_4')
            sage: h = T.mapping_class('s_0S_1s_2S_3s_1S_2')
            sage: h.is_pseudo_anosov()
            True
            sage: VeeringTriangulation.from_pseudo_anosov(h)
            VeeringTriangulation("(0,3,~2)(1,~0,~4)(2,5,4)(~1,~5,~3)", "BBBRRB")
        """
        FS = h.flat_structure()
        n = FS.triangulation.zeta

        X = {i.label: e.x for i,e in FS.edge_vectors.iteritems()}
        Y = {i.label: e.y for i,e in FS.edge_vectors.iteritems()}

        triangles = [[x.label for x in t] for t in FS.triangulation]
        colours = [RED if X[e]*Y[e] > 0 else BLUE for e in range(n)]
        return VeeringTriangulation(triangles, colours)

    @classmethod
    def from_square_tiled(cls, s, col=RED):
        r"""
        Build a veering triangulation from a square-tiled surface (or origami).

        INPUT:

        - ``s`` - a square-tiled surface

        - ``col`` - either ``RED`` or ``BLUE``

        EXAMPLES::

            sage: from surface_dynamics import Origami
            sage: from veerer import VeeringTriangulation

            sage: o = Origami('(1,2)', '(1,3)')
            sage: T = VeeringTriangulation.from_square_tiled(o)
            sage: T
            VeeringTriangulation("(0,1,2)(3,4,5)(6,7,8)(~8,~0,~7)(~6,~1,~5)(~4,~2,~3)", "RBBRBBRBB")
            sage: o.stratum()
            H_2(2)
            sage: T.stratum()
            H_2(2)

        A one cylinder example in the odd component of H(4)::

            sage: o = Origami('(1,2,3,4,5)', '(1,4,3,5,2)')
            sage: T = VeeringTriangulation.from_square_tiled(o)
            sage: o.stratum()
            H_3(4)
            sage: T.stratum()
            H_3(4)
        """
        require_package('surface_dynamics', 'from_square_tiled')
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

        return cls.from_face_edge_perms(array('l', colouring),
                                        array('l', fp),
                                        array('l', ep))

    @classmethod
    def from_stratum(cls, c, folded_edges=False):
        r"""
        Return a Veering triangulation from either a stratum, a component
        of stratum or a cylinder diagram.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: from veerer import *

            sage: T = VeeringTriangulation.from_stratum(AbelianStratum(2))
            sage: T
            VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: T.stratum()
            H_2(2)

            sage: Q = QuadraticStratum(9,-1)
            sage: CTreg = VeeringTriangulation.from_stratum(Q.regular_component())
            sage: CTreg.stratum()
            Q_3(9, -1)

            sage: CTirr = VeeringTriangulation.from_stratum(Q.irregular_component())
            sage: CTirr.stratum()
            Q_3(9, -1)

        Some examples built from cylinder diagram::

            sage: c = QuadraticCylinderDiagram('(0)-(0)')
            sage: VeeringTriangulation.from_stratum(c)
            VeeringTriangulation("(0,2,~1)(1,~0,~2)", "RBB")

            sage: c = QuadraticCylinderDiagram('(0,0)-(1,1,2,2,3,3)')
            sage: VeeringTriangulation.from_stratum(c)
            VeeringTriangulation("(0,10,~9)(1,~8,9)(2,~6,7)(3,~4,5)(4,~3,~11)(6,~2,~5)(8,~1,~7)(11,~10,~0)", "RRRRBBBBBBBB")

            sage: c = CylinderDiagram('(0,6,4,5)-(3,6,5) (1,3,2)-(0,1,4,2)')
            sage: CT = VeeringTriangulation.from_stratum(c)
            sage: CT.stratum()
            H_4(6)
        """
        require_package('surface_dynamics', 'from_stratum')

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

        colors = [RED] * nseps + [BLUE] * (2*nseps)
        return VeeringTriangulation(triangles, colors)

    @classmethod
    def from_face_edge_perms(self, colouring, fp, ep, vp=None, check=True):
        T = VeeringTriangulation.__new__(VeeringTriangulation)
        T._n = len(fp)
        T._fp = fp
        T._ep = ep

        if vp is None:
            n = T._n
            fp = T._fp
            ep = T._ep
            vp = array('l', [-1] * n)
            for i in range(n):
                vp[fp[ep[i]]] = i
        T._vp = vp
        T._colouring = colouring

        if check:
            T._check(ValueError)

        return T

    @classmethod
    def from_string(cls, s, check=True):
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
        cols = array('l', [colour_from_char(c) for c in cols])
        return VeeringTriangulation.from_face_edge_perms(cols, fp, ep, check=check)

    def to_flipper(self):
        r"""
        Return the corresponding flipper triangulation

        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: T.to_flipper()
            [(~2, ~0, ~1), (0, 1, 2)]
        """
        require_package('flipper', 'to_flipper')
        ep = self._ep
        F = []
        for f in self.faces():
            face = []
            for e in f:
                if ep[e] == e:
                    raise ValueError("flipper do not accept folded edges")
                if ep[e] < e:
                    face.append(~int(ep[e]))
                else:
                    face.append(int(e))
            F.append(tuple(face))

        return flipper.create_triangulation(F)

    def to_curver(self):
        r"""
        Return the corresponding curver triangulation

        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: T.to_curver()
            [(~2, ~0, ~1), (0, 1, 2)]
        """
        require_package('curver', 'to_curver')
        ep = self._ep
        F = []
        for f in self.faces():
            face = []
            for e in f:
                if ep[e] == e:
                    raise ValueError("curver do not accept folded edges")
                if ep[e] < e:
                    face.append(~int(ep[e]))
                else:
                    face.append(int(e))
            F.append(tuple(face))

        return curver.create_triangulation(F)

    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError
        return Triangulation.__eq__(self, other) and self._colouring == other._colouring

    def __ne__(self, other):
        r"""
        TESTS::

            sage: from veerer import *

            sage: VeeringTriangulation("(0,1,2)", 'RRB') != VeeringTriangulation("(0,1,2)", 'RRB')
            False
            sage: VeeringTriangulation("(0,1,2)", 'RRB') != VeeringTriangulation("(0,1,2)", 'RBR')
            True

        """
        if type(self) != type(other):
            raise TypeError
        return Triangulation.__ne__(self, other) or self._colouring != other._colouring

    def to_core(self, slope=VERTICAL):
        r"""
        Change the color of each forward (reps. backward) flippable edge in
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
                if self.is_forward_flippable(e):
                    self._colouring[e] = self._colouring[ep[e]] = COLOR
            else:
                if self.is_backward_flippable(e):
                    self._colouring[e] = self._colouring[ep[e]] = COLOR

    def copy(self):
        r"""
        Return a copy of this coloured triangulation

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], "RRB")
            sage: S1 = T.copy()
            sage: S2 = T.copy()
            sage: T == S1 == S2
            True
            sage: S1.flip(1,BLUE)
            sage: T == S1
            False
            sage: T == S2
            True
        """
        T = VeeringTriangulation.__new__(VeeringTriangulation)
        T._n = self._n
        T._vp = self._vp[:]
        T._ep = self._ep[:]
        T._fp = self._fp[:]
        T._colouring = self._colouring[:]
        return T

    def _faces_string(self):
        return ''.join('(' + ','.join(self._edge_rep(e) for e in f) + ')' for f in self.faces())

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
                self._faces_string(), self._colouring_string(short=True))

    def __repr__(self):
        return str(self)

    def colour(self, e):
        e = int(e)
        return self._colouring[e]

    def angles(self):
        r"""
        Return the list of angles (divided by \pi).

        Note that the vertices of the triangulation are labeled. The
        angles are given in the same order.

        # ??? Confused - there is no labelling yet.

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
        """
        n = self.num_half_edges()
        angles = []
        seen = [False] * n
        vp = self.vertex_permutation(copy=False)
        for e in range(n):
            if seen[e]: continue

            a = 0
            col = self._colouring[e]
            while not seen[e]:
                seen[e] = True
                ee = vp[e]
                ccol = self._colouring[ee]
                a += col != ccol
                e = ee
                col = ccol
            a += col != ccol

            angles.append(a/2)

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
            (True, [True, True, True, True, True, True, False, False, False, False, False, False])

            sage: T = VeeringTriangulation([(-12, 4, -4), (-11, -1, 11), (-10, 0, 10),
            ....:      (-9, 9, 1), (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)],
            ....:      [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE])
            sage: T.is_abelian()
            False
        """
        if self.num_folded_edges() > 0:
            return False

        ep = self._ep
        vp = self._vp

        # Perform BFS to check that we can coherently orient the
        # imaginary part of each edge
        oris = [None] * self._n
        oris[0] = True   # = going up
        oris[ep[0]] = False # = going down
        q = [0, ep[0]]
        while q:
            e0 = q.pop()
            o = oris[e0]

            e = e0
            f = vp[e]
            while f != e0:
                if self._colouring[e] == BLUE and self._colouring[f] == RED:
                    o = not o
                oris[f] = o
                if oris[ep[f]] is None:
                    q.append(ep[f])
                    oris[ep[f]]= not o
                elif oris[ep[f]] == oris[f]:
                    return (False, None) if certificate else False
                e,f = f, vp[f]

        return (True, oris) if certificate else True

    def stratum(self):
        r"""
        Return the Abelian or quadratic stratum of this coloured triagulation.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.stratum()
            H_1(0)

            sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
            sage: cols = 'BRBRBBBRBR'
            sage: T = VeeringTriangulation(fp, cols)
            sage: T.stratum()
            Q_1(1^2, -1^2)
        """
        require_package('surface_dynamics', 'stratum')

        A = self.angles()
        if any(a%2 for a in A) or not self.is_abelian():
            from surface_dynamics import QuadraticStratum
            return QuadraticStratum([(a-2) for a in A])
        else:
            from surface_dynamics import AbelianStratum
            return AbelianStratum([(a-2)/2 for a in A])

    def stratum_dimension(self):
        # each folded edge gives a pole
        dim1 = 2*self.genus() - 2 \
               + self.num_vertices() \
               + self.num_folded_edges() \
               + (1 if self.is_abelian() else 0)
        dim2 = self.stratum().dimension()
        assert dim1 == dim2
        return dim1

    def colours_about_edge(self, e):
        e = int(e)
        return [self._colouring[f] for f in self.square_about_edge(e)]

    def alternating_square(self, e):
        e = int(e)
        colours = self.colours_about_edge(e)
        if any(colours[f] == GREEN or colours[f] == PURPLE for f in range(4)):
            print("Warning: alternating_square is not carefullly defined with GREEN/PURPLE edges")
        return all(colours[f] != colours[(f+1) % 4] for f in range(4))

    def branches(self, slope=VERTICAL):
        r"""
        Return a 3-uplet made of the lists of respectively the small, mixed
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


    def is_flippable(self, e):
        r"""
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
        e = int(e)
        return Triangulation.is_flippable(self, e) and self.alternating_square(e)

    def flippable_edges(self):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.flippable_edges()
            [0, 1]
        """
        n = self._n
        ep = self._ep
        return [e for e in range(n) if e <= ep[e] and self.is_flippable(e)]

    def is_forward_flippable(self, e):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.is_forward_flippable(0)
            False
            sage: T.is_forward_flippable(1)
            True
            sage: T.is_forward_flippable(2)
            False
        """
        return Triangulation.is_flippable(self, e) and self.colours_about_edge(e) == [BLUE, RED, BLUE, RED]

    def forward_flippable_edges(self):
        r"""
        Return the set of forward flippable edges.

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
            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], "RRB")
            sage: for _ in range(10):
            ....:     rel = perm_random(6)
            ....:     T.relabel(rel)
            ....:     assert len(T.forward_flippable_edges()) == 1
        """
        ep = self._ep
        n = self._n
        return [e for e in range(n) if e <= ep[e] and self.is_forward_flippable(e)]

    def is_backward_flippable(self, e):
        r"""
        Test whether the edge ``e`` is backward flippable.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.is_backward_flippable(0)
            True
            sage: T.is_backward_flippable(1)
            False
            sage: T.is_backward_flippable(2)
            False
        """
        return Triangulation.is_flippable(self, e) and self.colours_about_edge(e) == [RED, BLUE, RED, BLUE]

    def backward_flippable_edges(self):
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
        return [e for e in range(n) if e <= ep[e] and self.is_backward_flippable(e)]

    def mostly_sloped_edges(self, slope):
        if slope == HORIZONTAL:
            return self.forward_flippable_edges()
        elif slope == VERTICAL:
            return self.backward_flippable_edges()
        else:
            raise ValueError

    def relabel(self, p, Lx=None, Gx=None):
        r"""
        Relabel inplace this veering triangulation according to the permutation ``p``.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, BLUE, BLUE])
            sage: T.relabel([0,1,3,2,5,4])
            sage: T
            VeeringTriangulation("(0,1,~2)(2,~0,~1)", "RBB")
            sage: T._check()

        Relabeling the subspace as well::

            sage: from veerer.permutation import perm_random_centralizer
            sage: T, s, t = VeeringTriangulations.L_shaped_surface(2, 3, 4, 5, 1, 2)
            sage: Gx = matrix(ZZ, [s, t])
            sage: for _ in range(10):
            ....:     p = perm_random_centralizer(T.edge_permutation(copy=False))
            ....:     T.relabel(p, Gx=Gx)
            ....:     T._set_switch_conditions(T._tt_check, Gx.row(0), VERTICAL)
            ....:     T._set_switch_conditions(T._tt_check, Gx.row(1), VERTICAL)

        Composing relabelings and permutation composition::

            sage: from veerer.permutation import perm_compose
            sage: from surface_dynamics import *
            sage: T0 = VeeringTriangulation.from_stratum(AbelianStratum(1,1,1,1))
            sage: for _ in range(10):
            ....:     p1 = perm_random_centralizer(T0.edge_permutation(copy=False))
            ....:     p2 = perm_random_centralizer(T0.edge_permutation(copy=False))
            ....:     T1 = T0.copy()
            ....:     T1.relabel(p1); T1.relabel(p2)
            ....:     T2 = T0.copy()
            ....:     T2.relabel(perm_compose(p1, p2))
            ....:     assert T1  == T2

        TESTS:

        This example used to be wrong::

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.relabel([1,5,0,2,4,3])
            sage: T.edge_colour(0) == BLUE
            True
            sage: T.edge_colour(1) == RED
            True
            sage: T._check()

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: from veerer.permutation import perm_random
            sage: for _ in range(10):
            ....:     r = perm_random(6)
            ....:     T.relabel(r)
            ....:     T._check()
        """
        n = self._n
        ep = self._ep
        if not perm_check(p, n):
            p = perm_init(p)
            if not perm_check(p, n):
                raise ValueError('invalid relabeling permutation')

        if Lx:
            raise NotImplementedError
        if Gx:
            seen = [False] * self.num_edges()
            for c in perm_cycles(p):
                c0 = self._norm(c[0])
                seen[c0] = True
                for i in range(1,len(c)):
                    ci = self._norm(c[i])
                    if seen[ci]:
                        break
                    seen[ci] = True
                    Gx.swap_columns(c0, ci)

        Triangulation.relabel(self, p)
        self._colouring = perm_on_array(p, self._colouring)

    def _automorphism_good_starts(self):
        r"""
        Start at a RED before a BLUE.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)", "RRB")
            sage: len(T._automorphism_good_starts())
            1
        """
        starts = []
        best_word = None

        n = self._n
        vp = self._vp
        cols = self._colouring

        # first run: compare edges using colors and vertex permutation orbits
        # FIX: Take advantage of folded edges, if they are present?
        for e in range(n):
            if cols[e] != RED:
                continue
            f = vp[e]
            if cols[f] != BLUE:
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
                # run through all the RED
                while cols[g] == RED:
                    g = vp[g]
                    n += 1
                w.append(n)

                if g == f or \
                   (best_word is not None and len(w) > len(best_word)):
                    break

            if best_word is None or \
               (len(w) < len(best_word) or (len(w) == len(best_word) and w < best_word)):
                starts = [e]
                best_word = w
            elif w == best_word:
                starts.append(e)

        if len(starts) == 1:
            return starts

        # try to break ties using face permutation orbits
        # (because the start edge e is RED before BLUE (on the vertex)
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

    def automorphisms(self):
        r"""
        Return the list of automorphisms of this veering triangulation.

        The output is a list of arrays that are permutations acting on the set
        of half edges.

        EXAMPLES::

            sage: from veerer import *

        An example with 4 symmetries in genus 2::

            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: V = VeeringTriangulation(fp, cols)
            sage: A = V.automorphisms()
            sage: len(A)
            4
            sage: S = V.copy()
            sage: for a in A:
            ....:     S.relabel(a)
            ....:     assert S == V
        """
        best = None
        best_relabellings = []
        for start_edge in self._automorphism_good_starts():
            relabelling = self._relabelling_from(start_edge)

            fp = perm_conjugate(self._fp, relabelling)
            ep = perm_conjugate(self._ep, relabelling)
            cols = perm_on_array(relabelling, self._colouring)

            T = (cols, fp, ep)
            if best is None or T == best:
                best_relabellings.append(relabelling)
                best = T
            elif T < best:
                del best_relabellings[:]
                best_relabellings.append(relabelling)
                best = T

        p0 = perm_invert(best_relabellings[0])
        return [perm_compose(p, p0) for p in best_relabellings]

    def edge_colour(self, e):
        return self._colouring[e]

    def best_relabelling(self, all=False):
        r"""
        """
        best = None
        if all:
            relabellings = []
        for start_edge in self._automorphism_good_starts():
            relabelling = self._relabelling_from(start_edge)

            # relabelled data
            fp = perm_conjugate(self._fp, relabelling)
            ep = perm_conjugate(self._ep, relabelling)
            cols = perm_on_array(relabelling, self._colouring)

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
            VeeringTriangulation("(0,1,2)(3,4,5)(~5,~3,~1)(~2,~0,~4)", "RBBRBB")
            sage: T._check()

            sage: from surface_dynamics import *
            sage: T0 = VeeringTriangulation.from_stratum(QuadraticStratum(1,1,1,1))
            sage: T = T0.copy()
            sage: T.rotate()
            sage: T.conjugate()
            sage: S = T0.copy()
            sage: S.conjugate()
            sage: S.rotate()
            sage: S == T
            True
        """
        self._colouring = array('l', [RED if x == BLUE else BLUE for x in self._colouring])

    def conjugate(self):
        r"""
        Conjugate this triangulation.

        EXAMPLES::

            sage: from veerer import *

            sage: faces = "(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)"
            sage: cols = [BLUE,RED,RED,BLUE,RED,RED]
            sage: T = VeeringTriangulation(faces, cols)
            sage: T.conjugate()
            sage: T
            VeeringTriangulation("(0,2,4)(1,3,5)(~5,~4,~3)(~1,~0,~2)", "RBBRBB")
            sage: T._check()
        """
        Triangulation.conjugate(self)
        self._colouring = array('l', [RED if x == BLUE else BLUE for x in self._colouring])

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
            sage: from surface_dynamics import *

            sage: T = VeeringTriangulation("(0,1,2)", "RRB")
            sage: T.to_string()
            'RRB_120_012'

            sage: T = VeeringTriangulation.from_stratum(QuadraticStratum({1:20}))
            sage: s = T.to_string()
            sage: TT = VeeringTriangulation.from_string(s)
            sage: T == TT
            True
        """
        colours = self._colouring_string(short=False)
        fp = perm_base64_str(self._fp)
        ep = perm_base64_str(self._ep)
        return colours + '_' + fp + '_' + ep

    def iso_sig(self, Lx=None, Gx=None):
        r"""
        Return a canonical string ("isomorphism signature").

        INPUT:

        - ``Lx`` - (optional) a matrix whose rows are the linear equations
          for the admissible train-track lengths.

        - ``Gx`` - (optional) a matrix whose rows are the generators of admissible
          train-track lengths

        EXAMPLES::

            sage: from veerer import *

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = VeeringTriangulation(t, cols)
            sage: T.iso_sig()
            'RBRBRBBBRBBBBBBRBBBRBRBR_7ad9e8fbhjl0mkig654321nc_nmlkjihgfedcba9876543210'

        If we relabel the triangulation, the isomorphic signature does not change::

            sage: from veerer.permutation import perm_random
            sage: p = perm_random(24)
            sage: T.relabel(p)
            sage: T.iso_sig()
            'RBRBRBBBRBBBBBBRBBBRBRBR_7ad9e8fbhjl0mkig654321nc_nmlkjihgfedcba9876543210'

        An isomorphic triangulation can be reconstructed from the isomorphic
        signature via::

            sage: s = T.iso_sig()
            sage: T2 = VeeringTriangulation.from_string(s)
            sage: T == T2
            False
            sage: T.is_isomorphic_to(T2)
            True

        TESTS::

            sage: from veerer.veering_triangulation import VeeringTriangulation
            sage: from veerer.permutation import perm_random

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = VeeringTriangulation(t, cols)
            sage: iso_sig = T.iso_sig()
            sage: for _ in range(10):
            ....:     p = perm_random(24)
            ....:     T.relabel(p)
            ....:     assert T.iso_sig() == iso_sig
        """
        n = self._n

        if Lx:
            R, (cols, fp, ep) = self.best_relabelling(True)
            raise NotImplementedError("not implemented for linear equations")
        if Gx:
            raise NotImplementedError("not implemented for generators")
        else:
            r, (cols, fp, ep) = self.best_relabelling()
            cols = ''.join(colour_to_char(col) for col in cols)

            fp = perm_base64_str(fp)
            ep = perm_base64_str(ep)

            return cols + '_' + fp + '_' + ep

    def set_canonical_labels(self, Lx=None, Gx=None):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = VeeringTriangulation(t, cols)
            sage: T.set_canonical_labels()
            sage: s0 = T.to_string()
            sage: T.set_canonical_labels()
            sage: assert s0 == T.to_string()
            sage: from veerer.permutation import perm_random_centralizer
            sage: for _ in range(10):
            ....:     p = perm_random_centralizer(T.edge_permutation(copy=False))
            ....:     T.relabel(p)
            ....:     T.set_canonical_labels()
            ....:     assert s0 == T.to_string()
        """
        if Lx:
            raise NotImplementedError
        elif Gx:
            R, (cols, fp, ep) = self.best_relabelling(all=True)
            self.relabel(R[0], Gx=Gx)
            if len(R) == 1:
                return
            for i in range(2, len(R)):
                pass
        else:
            r, (cols, fp, ep) = self.best_relabelling()
            self.relabel(r)

    def canonical(self):
        r"""
        Return a canonically labeled isomorphic coloured triangulation.
        """
        raise ValueError

    def is_isomorphic_to(self, other):
        r"""
        Test whether this triangulation is isomorphic to ``other``.
        """
        return self.iso_sig() == other.iso_sig()

# TODO?
#    def isometries_to(self, other):
#        return [isom for isom in self._triangulation.isometries_to(other._triangulation) if all(self._colouring[i] == other._colouring[isom.index_map[i]] for i in range(self.zeta))]

# TODO?
#    def self_isometries(self):
#        return self.isometries_to(self)

    def flip(self, e, col, Lx=None, Gx=None):
        r"""
        Flip an edge inplace.

        INPUT:

        - ``e`` - edge number

        - ``col`` - color of the edge after the flip (ie either ``RED`` or ``BLUE``)

        - ``Lx`` - (optional) - matrix whose rows are equations in a linear subspace
          that has to be carried around

        - ``Gx`` - (optional) - matrix whose rows are generators of a linear subspace
          that has to be carried around

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.flip(1, RED); T
            VeeringTriangulation("(0,~2,1)(2,~1,~0)", "RRB")
            sage: T.flip(0, RED); T
            VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB")

            sage: T.flip(1, BLUE); T
            VeeringTriangulation("(0,~2,1)(2,~1,~0)", "RBB")
            sage: T.flip(2, BLUE); T
            VeeringTriangulation("(0,~1,~2)(1,2,~0)", "RBB")

        Some examples involving linear subspaces::

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: Gx = matrix(ZZ, [s, t])
            sage: T.flip(3, 2, Gx=Gx)
            sage: T.flip(4, 2, Gx=Gx)
            sage: T.flip(5, 2, Gx=Gx)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(0), VERTICAL)
            sage: T._set_switch_conditions(T._tt_check, Gx.row(1), VERTICAL)

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(2, 3, 4, 5, 1, 2)
            sage: Gx = matrix(ZZ, [s, t])
            sage: flip_sequence = [(3, 2), (4, 1), (5, 2), (6 , 2), (5, 1), (1, 1), (5, 1)]
            sage: for e, col in flip_sequence:
            ....:     T.flip(e, col, Gx=Gx)
            ....:     T._set_switch_conditions(T._tt_check, Gx.row(0), VERTICAL)
            ....:     T._set_switch_conditions(T._tt_check, Gx.row(1), VERTICAL)
        """
        if Lx:
            raise NotImplementedError("not implemented for linear equations")
        if Gx:
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
                raise ValueError("invalid color")

        # topological flip
        e = int(e)
        assert(self.is_flippable(e))
        E = self._ep[e]

        Triangulation.flip(self, e)
        self._colouring[e] = self._colouring[E] = col

    def cylinders(self, col):
        r"""
        Return the cylinders of color ``col``.

        Each cylinder is given as a list of edges that correspond
        to the (ordered) list of edges crossed by the waist curve
        of the cylinder.

        EXAMPLES::

            sage: from veerer import *

        The torus::

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.cylinders(RED)
            [[0, 4]]
            sage: T.cylinders(BLUE)
            []

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "BRB")
            sage: T.cylinders(BLUE)
            [[0, 3]]
            sage: T.cylinders(RED)
            []

        Some examples in Q(4, -1^8)::

            sage: T = VeeringTriangulation("(0,1,2)", "BBR")
            sage: T.cylinders(BLUE)
            [[0, 1]]
            sage: T.cylinders(RED)
            []

            sage: T = VeeringTriangulation("(5,4,7)(~5,3,10)(1,~0,8)(~1,~4,11)(2,6,9)(~2,0,12)", "BBBBBBBRRRRRR")
            sage: T.cylinders(BLUE)
            [[3, 17, 4, 15, 16, 13, 6]]
            sage: T.cylinders(RED)
            []

            sage: fp = "(0,12,3)(1,2,11)(4,6,7)(5,~2,10)(8,~5,~4)(9,~1,~0)"
            sage: cols = "BBRRRBBBBRRRB"
            sage: T = VeeringTriangulation(fp, cols)
            sage: T.cylinders(BLUE)
            [[6, 7]]
            sage: T.cylinders(RED)
            [[10, 17, 11]]

            sage: fp = '(0,~5,4)(1,~3,2)(3,5,~4)(6,~1,~0)'
            sage: cols = 'RBBRBRB'
            sage: T = VeeringTriangulation(fp, cols)
            sage: T.cylinders(BLUE)
            [[6, 8, 2]]
            sage: T.cylinders(RED)
            []

            sage: T = VeeringTriangulation.from_string("RBRBBBRRBRBR_908274a531b6_ba9385764210")
            sage: T.cylinders(BLUE)
            [[5, 4, 3]]
        """
        n = self._n
        fp = self._fp
        ep = self._ep
        cols = self._colouring

        cylinders = []

        seen = [False] * n
        for a in range(n):
            if seen[a]:
                continue

            # compute the triangle (a,b,c)
            b = fp[a]
            c = fp[b]
            if seen[a] or seen[b] or seen[c] or \
               (cols[a] == col) + (cols[b] == col) + (cols[c] == col) != 2:
                seen[a] = seen[b] = seen[c] = True
                continue


            # find the two edges that can be crossed
            doors = []
            for e in (a,b,c):
                if cols[e] == col:
                    doors.append(e)
            door = door0 = doors[0]

            cc = []  # cycle to be stored
            half_turn = False
            cyl = True
            while not seen[door]:
                seen[door] = seen[ep[door]] = True
                cc.append(door)

                if door == ep[door]:
                    # folded edges... we can not continue from here, start from
                    # the other door or stop
                    if half_turn:
                        cyl = True
                        break
                    half_turn = True
                    cc = [ep[x] for x in reversed(cc)]
                    door = doors[1]
                    continue

                # pass through the door
                a = ep[door]
                b = fp[a]
                c = fp[b]
                if cols[b] == col:
                    door = b
                elif cols[c] == col:
                    door = c
                else:
                    cyl = False
                    break

            if cyl and (door == door0 or half_turn):
                cylinders.append(cc)

        return cylinders

    def is_cylindrical(self, col=None):
        r"""
        Return whether this veering triangulation is cylindrical.

        A Veering triangulation is blue cylindrical (resp red cylindrical) if
        all its triangles have two blue edges (resp red).

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.is_cylindrical()
            True
        """
        if col is None:
            return self.is_cylindrical(BLUE) or self.is_cylindrical(RED)

        n = self._n
        fp = self._fp
        cols = self._colouring
        seen = [False] * n
        for a in range(n):
            if seen[a]:
                continue
            b = fp[a]
            c = fp[b]
            if (cols[a] == col) + (cols[b] == col) + (cols[c] == col) != 2:
                return False
        return True

    def is_quadrangulable(self):
        k = 0
        n = self.num_edges()
        ep = self._ep
        for e in range(n):
            if not self.is_forward_flippable(e):
                continue

            k += 1 if ep[e] == e else 2

        return k == self.num_faces()

    def is_square_tiled(self, col=None):
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
        if col is None:
            return self.is_square_tiled(BLUE) or self.is_square_tiled(RED)

        k = 0
        n = self.num_edges()
        ep = self._ep
        for e in range(n):
            if not self.is_forward_flippable(e):
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
            81
            sage: properties_to_string(81)
            'red geometric h-balanced'
        """
        from .constants import (BLUE, RED, SQUARETILED, SQUARETILED, QUADRANGULABLE,
            GEOMETRIC, VBALANCED, HBALANCED)

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
        if self.is_balanced(VERTICAL):
            code |= VBALANCED
        if self.is_balanced(HORIZONTAL):
            code |= HBALANCED

        if code & BLUE and code & RED:
            raise RuntimeError("found a blue and red triangulations!")
        if code & SQUARETILED:
            if not code & BLUE and not code & RED:
                raise RuntimeError("square-tiled should be colored")
            if not code & GEOMETRIC:
                raise RuntimeError("square-tiled should be geometric")
            if not code & VBALANCED:
                raise RuntimeError("square-tiled should be v-balanced")
            if not code & HBALANCED:
                raise RuntimeError("square-tiled should be h-balanced")

        return code

    def flip_back(self, e, col):
        r"""
        Flip backward an edge in place

        EXAMPLES::

            sage: from veerer import *

            sage: T0 = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T = T0.copy()
            sage: T.flip(1, RED); T.flip(0, RED)
            sage: T.flip_back(0, RED); T.flip_back(1, RED)
            sage: T == T0
            True

            sage: T.flip(1, BLUE); T.flip(2, BLUE)
            sage: T.flip_back(2, BLUE); T.flip_back(1, RED)
            sage: T == T0
            True
        """
        e = int(e)
        assert(self.is_flippable(e))
        E = self._ep[e]

        Triangulation.flip_back(self, e)
        self._colouring[e] = self._colouring[E] = col

    def _set_switch_conditions(self, insert, x, slope=VERTICAL):
        r"""
        These are the linear parts of the train-track equations

        INPUT:

        - ``insert`` - a function to be called for each equation

        - ``x`` - variable factory (the variable for edge ``e`` is constructed
          via ``x[e]``)

        - ``slope`` - (default ``VERTICAL``) the slope of the train-track
          ``HORIZONTAL`` or ``VERTICAL``
        """
        if slope == VERTICAL:
            POS = BLUE
            NEG = RED
        elif slope == HORIZONTAL:
            POS = RED
            NEG = BLUE
        else:
            raise ValueError('bad slope parameter')

        for (i,j,k) in self.faces():
            # find the large edge
            # if horizontal, this is the one opposite to the RED/BLUE transition
            # if vertical this is the one opposite to the BLUE/RED transition
            i = self._norm(i)
            j = self._norm(j)
            k = self._norm(k)
            if self._colouring[i] == POS and self._colouring[j] == NEG:
                # k is large
                l,s1,s2 = k,i,j
            elif self._colouring[j] == POS and self._colouring[k] == NEG:
                # i is large
                l,s1,s2 = i,j,k
            elif self._colouring[k] == POS and self._colouring[i] == NEG:
                # j is large
                l,s1,s2 = j,k,i
            else:
                raise ValueError('can not determine the big edge on triangle (%s, %s, %s)' %
                                 (self_.edge_rep(i), self._edge_rep(j), self._edge_rep(k)))

            insert(x[self._norm(l)] == x[self._norm(s1)] + x[self._norm(s2)])

    @staticmethod
    def _tt_check(x):
        if not x:
            raise AssertionError("does not satisfy train-track contraints".format(x))

    def _set_train_track_constraints(self, insert, x, slope, low_bound, allow_degenerations):
        r"""
        Sets the equation and inequations for train tracks.

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
            AssertionError: does not satisfy train-track contraints

        Check equations with folded edges (that are "counted twice")::

            sage: T = VeeringTriangulation("(0,2,3)(~0,1,4)(~1,5,6)", [BLUE, RED, RED, BLUE, BLUE, BLUE, BLUE])
            sage: T._set_train_track_constraints(T._tt_check, [0,1,1,1,1,1,0], VERTICAL, False, False)
            sage: T._set_train_track_constraints(T._tt_check, [1,2,3,4,3,7,5], VERTICAL, False, False)
        """
        if slope == VERTICAL:
            POS = BLUE
            NEG = RED
        elif slope == HORIZONTAL:
            POS = RED
            NEG = BLUE
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
            if not low_bound or \
                (allow_degenerations and \
                 ((slope == HORIZONTAL and self.is_forward_flippable(e)) or \
                 (slope == VERTICAL and self.is_backward_flippable(e)))):
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

    def _set_balance_constraints(self, insert, x, slope):
        r"""
        """
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

    def GL2R_span(self, lx, ly):
        require_package('ppl', 'GL2R_span')

        ne = self.num_edges()

        x = [ppl.Variable(i) for i in range(ne)]
        y = [ppl.Variable(ne+i) for i in range(ne)]

        gs = ppl.Generator_System()
        gs.insert(ppl.point())

        return ppl.C_Polyhedron(gs)

    def train_track_linear_space(self, slope=VERTICAL):
        r"""
        Return the polytope determined by the switch equations (a linear subspace)

        INPUT:

        - ``slope`` - the slope for the train track (``HORIZONTAL`` or ``VERTICAL``)

        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB")
            sage: T.train_track_linear_space()
            A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 1 point, 2 lines
        """
        require_package('ppl', 'train_track_linear_space')

        cs = ppl.Constraint_System()
        ne = self.num_edges()
        x = [ppl.Variable(e) for e in range(ne)]
        self._set_switch_conditions(cs.insert, x, slope)
        return ppl.C_Polyhedron(cs)

    def train_track_polytope(self, slope=VERTICAL, low_bound=0):
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
            A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 1 point, 2 rays
            sage: P.generators()
            Generator_System {point(0/1, 0/1, 0/1), ray(1, 1, 0), ray(0, 1, 1)}

            sage: P = T.train_track_polytope(VERTICAL, low_bound=3)
            sage: P.generators()
            Generator_System {ray(1, 1, 0), ray(0, 1, 1), point(3/1, 6/1, 3/1)}
        """
        require_package('ppl', 'train_track_polytope')

        cs = ppl.Constraint_System()
        ne = self.num_edges()
        variables = [ppl.Variable(e) for e in range(ne)]
        self._set_train_track_constraints(cs.insert, variables, slope, low_bound, False)
        return ppl.C_Polyhedron(cs)

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
        require_package('ppl', 'train_track_min_solution')

        n = self.num_edges()
        M = ppl.MIP_Problem(n)

        x = [ppl.Variable(e) for e in range(n)]
        M.set_objective_function(-sum(x))
        self._set_train_track_constraints(M.add_constraint, x, slope, 1, allow_degenerations)
        return M.optimizing_point()

    def geometric_polytope(self, Lx=None, Gx=None, x_low_bound=0, y_low_bound=0, hw_bound=0):
        r"""
        Return the geometric polytope of this veering triangulation.

        The geometric polytope is the polytope of length and heights data that
        corresponds to L-infinity Delaunay triangulations.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.geometric_polytope()
            A 4-dimensional polyhedron in QQ^6 defined as the convex hull of 1 point, 7 rays
            sage: T.geometric_polytope(x_low_bound=1, y_low_bound=1, hw_bound=1)
            A 4-dimensional polyhedron in QQ^6 defined as the convex hull of 1 point, 7 rays

        TESTS::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: T = VeeringTriangulation.from_stratum(AbelianStratum(1,1))
            sage: for _ in range(100):
            ....:     T.random_forward_flip()
            ....:     test1 = T.geometric_polytope().affine_dimension() == 10
            ....:     test2 = not T.geometric_polytope(x_low_bound=1,y_low_bound=1,hw_bound=1).is_empty()
            ....:     assert test1 == test2, T.to_string()

            sage: T = VeeringTriangulation.from_stratum(QuadraticStratum(1,1,1,1))
            sage: for _ in range(100):
            ....:     T.random_forward_flip()
            ....:     test1 = T.geometric_polytope().affine_dimension() == 12
            ....:     test2 = not T.geometric_polytope(x_low_bound=1,y_low_bound=1,hw_bound=1).is_empty()
            ....:     assert test1 == test2, T.to_string()
        """
        require_package('ppl', 'geometric_polytope')

        ne = self.num_edges()
        x = [ppl.Variable(e) for e in range(ne)]
        y = [ppl.Variable(ne+e) for e in range(ne)]

        # 1. train-track conditions
        if Lx is None and Gx is None:
            P = self.train_track_polytope(VERTICAL, low_bound=x_low_bound)
            P.concatenate_assign(self.train_track_polytope(HORIZONTAL, low_bound=y_low_bound))
        else:
            if Lx is None:
                # construct linear equations from given generators
                Gy = self._complexify_generators(Gx)
                gs = ppl.Generator_System()
                gs.insert(ppl.point())
                for g in Gx:
                    g = sum(coeff * x[e] for e,coeff in enumerate(g))
                    gs.insert(ppl.line(g))
                for g in Gy:
                    g = sum(coeff * y[e] for e,coeff in enumerate(g))
                    gs.insert(ppl.line(g))
                P = ppl.C_Polyhedron(gs)

            elif Gx is None:
                # construct linear equations from given equations
                Ly = self._complexify_equations(Lx)
                cs = ppl.Constraint_System()
                for l in Lx:
                    l *= l.denominator()
                    l = sum(coeff * x[e] for e,coeff in enumerate(l))
                    cs.insert(l == 0)
                for l in Ly:
                    l *= l.denominator()
                    l = sum(coeff * y[e] for e,coeff in enumerate(l))
                    cs.insert(l == 0)
                P = ppl.C_Polyhedron(cs)

            else:
                raise ValueError("only one of Gx or Lx can be used")

            # non-negativity
            for e in range(ne):
                P.add_constraint(x[e] >= x_low_bound)
                P.add_constraint(y[e] >= y_low_bound)

        # geometric constraints
        self._set_geometric_constraints(P.add_constraint, x, y, hw_bound=hw_bound)

        return P

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
        require_package('sage', '_complexify_generators')

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
                    Gy[i,j] *= -1
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
        require_package('sage', '_complexify_equations')

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
                    Ly[i,j] *= -1
        return Ly

    def _flat_structure_from_train_track_lengths(self, VH, VV, base_ring=None):
        r"""
        Return a flat structure from two vectors ``VH`` and ``VV``
        satisfying the train track equations.
        """
        require_package('sage', '_flat_structure_from_train_track_lengths')

        from sage.modules.free_module import VectorSpace

        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ

        assert len(VH) == len(VV) == self.num_edges()
        assert all(x >=0 for x in VH)
        assert all(x >= 0 for x in VV)

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

        from .layout import FlatVeeringTriangulationLayout
        return FlatVeeringTriangulationLayout(self, vectors)

    def flat_structure_middle(self):
        r"""
        Return a flat structure with this Veering triangulation.

        Note that this triangulation must be core. The point is chosen
        by taking the interior point of the polytope obtained by
        summing each ray.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)", "RRB")
            sage: T.flat_structure_middle()
            Flat Triangulation made of 1 triangles

            sage: from surface_dynamics import *
            sage: Q = QuadraticStratum({1:4, -1:4})
            sage: CT = VeeringTriangulation.from_stratum(Q)
            sage: CT.flat_structure_middle()
            Flat Triangulation made of 16 triangles
        """
        n = self.num_edges()

        PH = self.train_track_polytope(HORIZONTAL)
        PV = self.train_track_polytope(VERTICAL)

        # pick sum of rays
        VH = [g.coefficients() for g in PH.generators() if g.is_ray()]
        VH = [sum(v[i] for v in VH) for i in range(n)]
        VV = [g.coefficients() for g in PV.generators() if g.is_ray()]
        VV = [sum(v[i] for v in VV) for i in range(n)]

        return self._flat_structure_from_train_track_lengths(VH, VV)

    def flat_structure_min(self, allow_degenerations=False):
        r"""
        Return a flat structure with this Veering triangulation.

        Note that this triangulation must be core. The point is chosen
        by taking the minimum integral point in the cone.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: from veerer import *

            sage: Q = QuadraticStratum({1:4, -1:4})
            sage: CT = VeeringTriangulation.from_stratum(Q)
            sage: CT.flat_structure_min()
            Flat Triangulation made of 16 triangles

        By allowing degenerations you can get a simpler solution but
        with some of the edges horizontal or vertical::

            sage: CT.flat_structure_min(True)
            Flat Triangulation made of 16 triangles
        """
        require_package('sage', 'flat_structure_min')

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

    def flat_structure_geometric_middle(self, L=None):
        r"""
        Return a geometric flat structure obtained by averaging the
        vertices of the geometric polytope.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.flat_structure_geometric_middle()
            Flat Triangulation made of 2 triangles
        """
        n = self.num_edges()


        P = self.geometric_polytope()
        if L is not None:
            P.intersection_assign(L)
        r = [g.coefficients() for g in P.generators() if g.is_ray()]
        VV = [sum(v[i] for v in r) for i in range(n)]
        VH = [sum(v[n+i] for v in r) for i in range(n)]

        return self._flat_structure_from_train_track_lengths(VH, VV)

    def geometric_flat_structure(self):
        raise NotImplementedError
        return self.geometric_polytope().vrep()[1:, 1:].sum(axis=0).tolist()

    def is_core(self, method='polytope'):
        r"""
        Test whether this coloured triangulation is core.

        INPUT:

        - ``method`` - a string which should either be ``'polytope'`` or ``'LP'``

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
            sage: T.flip(14,RED)
            sage: T.is_core()
            True
            sage: T.flip(14,RED)
            sage: T.is_core()
            False
        """
        # In theory LP should be much faster but in practice (in small dimensions)
        # polytope is much better
        if method == 'polytope':
            d = self.stratum().dimension()
            return self.train_track_polytope(HORIZONTAL).affine_dimension() == d and \
                   self.train_track_polytope(VERTICAL).affine_dimension() == d
        elif method == 'LP':
            require_package('ppl', 'is_core')

            n = self.num_edges()
            x = [ppl.Variable(e) for e in range(n)]
            for slope in (HORIZONTAL, VERTICAL):
                M = ppl.MIP_Problem(n)
                self._set_train_track_constraints(M.add_constraint, x, slope, True, False)

                if M.solve()['status'] == 'unfeasible':
                    return False

            return True

        else:
            raise ValueError("method must be either 'polytope' or 'LP'")

    def is_geometric(self, method=None):
        r"""
        Test whether this coloured triangulation is geometric.

        INPUT:

        - ``method`` - string - either ``'polytope'`` or ``'LP'`` (linear
          programming).

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: from veerer import *

        TESTS::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: T = VeeringTriangulation.from_stratum(AbelianStratum(1,1))
            sage: for _ in range(100):
            ....:     T.random_forward_flip()
            ....:     test1 = T.is_geometric('LP')
            ....:     test2 = T.is_geometric('polytope')
            ....:     test3 = T.is_geometric('polytope2')
            ....:     assert test1 == test2 == test3, T.to_string()

            sage: T = VeeringTriangulation.from_stratum(QuadraticStratum(1,1,1,1))
            sage: for _ in range(100):
            ....:     T.random_forward_flip()
            ....:     test1 = T.is_geometric('LP')
            ....:     test2 = T.is_geometric('polytope')
            ....:     test3 = T.is_geometric('polytope2')
            ....:     assert test1 == test2 == test3, T.to_string()
        """
        if method is None or method == 'polytope':
            d = self.stratum().dimension()
            return self.geometric_polytope().affine_dimension() == 2*d
        elif method == 'polytope2':
            return not self.geometric_polytope(x_low_bound=1,y_low_bound=1,hw_bound=1).is_empty()
        elif method == 'LP':
            require_package('ppl', 'is_geometric')

            ne = self.num_edges()
            M = ppl.MIP_Problem(2*ne)
            x = [ppl.Variable(e) for e in range(ne)]
            y = [ppl.Variable(e) for e in range(ne, 2*ne)]

            self._set_train_track_constraints(M.add_constraint, x, VERTICAL, True, False)
            self._set_train_track_constraints(M.add_constraint, y, HORIZONTAL, True, False)
            self._set_geometric_constraints(M.add_constraint, x, y, True)

            return M.solve()['status'] != 'unfeasible'

        else:
            raise ValueError("method must be 'polytope', 'polytope2' or 'LP'")

    def is_balanced(self, slope=VERTICAL):
        r"""
        Check balanceness

        EXAMPLES::

            sage: from veerer import VeeringTriangulation

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RBR")
            sage: T.is_balanced()
            True

            sage: T = VeeringTriangulation("(0,~3,2)(1,4,~2)(3,5,~4)(~5,~1,~0)", "RBBBRB")
            sage: T.is_balanced()
            False

            sage: T = VeeringTriangulation("(0,1,8)(2,~7,~1)(3,~0,~2)(4,~5,~3)(5,6,~4)(7,~8,~6)", "BRRRRBRBR")
            sage: T.is_balanced()
            False
            sage: T.rotate()
            sage: T.is_balanced()
            True
        """
        ne = self.num_edges()
        cs = ppl.Constraint_System()
        x = [ppl.Variable(e) for e in range(ne)]
        self._set_train_track_constraints(cs.insert, x, slope, False, False)
        self._set_balance_constraints(cs.insert, x, slope)
        P = ppl.C_Polyhedron(cs)

        return P.affine_dimension() == self.stratum().dimension()

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

        Flipping edge 1 in RED is fine (it remains a core triangulation)::

            sage: T = T0.copy()
            sage: T.flip(1, RED)
            sage: T.edge_has_curve(1)
            True
            sage: T.is_core()
            True

        However, flipping edge 1 in BLUE leads to a non-core triangulation
        (both edges 1 and 3 degenerate)::

            sage: T = T0.copy()
            sage: T.flip(1,BLUE)
            sage: T.edge_has_curve(1)
            False
            sage: T.edge_has_curve(3)
            False

            sage: PH = T.train_track_polytope(HORIZONTAL)
            sage: PH.affine_dimension()
            2
            sage: [g.coefficients() for g in PH.generators() if g.is_ray()]
            [(mpz(0), mpz(0), mpz(0), mpz(0), mpz(1), mpz(1)),
             (mpz(1), mpz(0), mpz(1), mpz(0), mpz(0), mpz(0))]
        """
        # TODO: we should only searching for vertex cylces; i.e. not allow more
        # than two pairs (i, ~i) to be both seen (barbell are fine but not more)

        # TODO: we want the._colouring to also take care of negative edges
        # (bad alternative: take "norm" function from flipper)
        colouring = self._colouring
        edge_rep = self._edge_rep
        ep = self._ep

        if verbose:
            print('[edge_has_curve] checking edge %s with color %s' % (edge_rep(e), colouring[e]))

        a, b, c, d = self.square_about_edge(e)
        if colouring[a] == BLUE:
            POS, NEG = BLUE, RED
            if verbose:
                print('[edge_has_curve] checking HORIZONTAL track')
        else:
            POS, NEG = RED, BLUE
            if verbose:
                print('[edge_has_curve] checking VERTICAL track')

        # check alternating condition
        assert colouring[a] != colouring[b]
        assert colouring[c] == colouring[a]
        assert colouring[d] == colouring[b]

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
        q = deque()
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

    def geometric_flips(self, Lx=None, Gx=None):
        r"""
        Return the list of geometric neighbours.

        INPUT:

        - ``Lx`` - (optional) a matrix whose rows are the linear equations
          for admissible train-track lengths (the right kernel is the space
          of admissible lengths). It is assumed that the linear equations
          contains the switch (or triangle) equations.

        - ``Gx`` - (optional) a matrix whose rows are admissible train-track
          lengths

        OUTPUT: a list of pairs (edge number, new color)

        EXAMPLES:

        L-shaped square tiled surface with 3 squares (given as a sphere with
        3 triangles). It has two geometric neighbors corresponding to simultaneous
        flipping of the diagonals 3, 4 and 5::

            sage: from veerer import *
            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1,1,1,1)
            sage: Gx = matrix(QQ, [s, t])
            sage: T.geometric_flips(Gx=Gx)
            [[(3, 2), (4, 2), (5, 2)], [(3, 1), (4, 1), (5, 1)]]

        To be compared with the geometric flips in the ambient stratum::

            sage: T.geometric_flips()
            [[(3, 2)], [(3, 1)], [(4, 2)], [(4, 1)], [(5, 2)], [(5, 1)]]

        Alternatively, instead of the generators ``Gx`` one can use the
        defining equations::

            sage: Lx = Gx.right_kernel_matrix()
            sage: T.geometric_flips(Lx=Lx)
            [[(3, 2), (4, 2), (5, 2)], [(3, 1), (4, 1), (5, 1)]]

        A more complicated example::

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(2,3,5,2,1,1)
            sage: Gx = matrix(QQ, [s, t])
            sage: T.geometric_flips(Gx=Gx)
            [[(5, 2)], [(5, 1)], [(4, 2)]]
            sage: Lx = Gx.right_kernel_matrix()
            sage: T.geometric_flips(Lx=Lx)
            [[(5, 2)], [(5, 1)], [(4, 2)]]
        """
        require_package('ppl', 'geometric_neighbours')

        ne = self.num_edges()
        x = [ppl.Variable(e) for e in range(ne)]
        y = [ppl.Variable(ne+e) for e in range(ne)]

        # construct the linear invariant subspace we are interested in
        if Lx:
            dim = ne - Lx.nrows()
        elif Gx:
            dim = Gx.nrows()
        else:
            dim = self.stratum().dimension()

        P = self.geometric_polytope(Gx=Gx, Lx=Lx)
        if P.affine_dimension() != 2 * dim:
            raise ValueError('not geometric or invalid constraints')

        # constructing the Delaunay facets
        # (this is only a subset of the facets)
        delaunay_facets = {}
        for e in self.forward_flippable_edges():
            a,b,c,d = self.square_about_edge(e)

            Q = ppl.C_Polyhedron(P)
            Q.add_constraint(x[self._norm(e)] == y[self._norm(a)] + y[self._norm(d)])
            if Q.affine_dimension() == 2*dim - 1:
                hQ = ppl_cone_to_hashable(Q)
                if hQ not in delaunay_facets:
                    delaunay_facets[hQ] = [Q, []]
                delaunay_facets[hQ][1].append(e)

        # computing colouring of the new triangulations
        neighbours = []
        for Q, edges in delaunay_facets.values():
            for cols in product([BLUE, RED], repeat=len(edges)):
                S = ppl.C_Polyhedron(Q)
                Z = list(zip(edges, cols))
                for e,col in Z:
                    a,b,c,d = self.square_about_edge(e)
                    if col == RED:
                        S.add_constraint(x[self._norm(a)] <= x[self._norm(d)])
                    else:
                        S.add_constraint(x[self._norm(a)] >= x[self._norm(d)])
                if S.affine_dimension() == 2*dim - 1:
                    neighbours.append(Z)

        return neighbours

    def _check_edge_has_curve(self):
        r"""
        Check the function ``edge_has_curve``

        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T._check_edge_has_curve()
        """
        dim = self.stratum().dimension()
        assert self.is_core()
        for slope in [HORIZONTAL, VERTICAL]:
            for e in self.mostly_sloped_edges(slope):
                for col in [BLUE, RED]:
                    T = self.copy()
                    T.flip(e, col)
                    test1 = T.edge_has_curve(e)
                    test2 = T.train_track_polytope(slope).affine_dimension() == dim
                    if test1 != test2:
                        T.edge_has_curve(e, verbose=True)
                        raise RuntimeError("failed\nT = {}\nedge = {}\ncolour = {}\nhas curve={}\nstratum dim={}\nhoriz tt dim={}\nvert tt dim={}".format(T, e, col, test1, dim,
                            T.train_track_polytope(HORIZONTAL).affine_dimension(),
                            T.train_track_polytope(VERTICAL).affine_dimension()))

    def random_forward_flip(self, repeat=1):
        r"""
        Apply a forward flip randomly among the ones that keeps the triangulation core.

        INPUT:

        - ``repeat`` - integer (default 1) - if provided make ``repeat`` flips instead of 1.
        """
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

class VeeringTriangulations(object):
    @staticmethod
    def L_shaped_surface(a1, a2, b1, b2, t1=0, t2=0):
        r"""
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
    def H11_rank1_locus():
        # return a surface in the McMullen discriminant locus in H(1,1)
        # for square discriminants it is defined over Q while for primitive ones
        # it is over a quadratic number field
        raise NotImplemented

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

