r"""
Veering triangulations of surfaces.
"""
from __future__ import print_function, absolute_import

from math import log
from collections import deque
from itertools import product

from .constants import *
from .even_permutation import *
from .misc import det2
from .triangulation import Triangulation, edge_label, norm

def dot(A, B):
    return sum(a*b for a, b in zip(A, B))

def best_rotation(X):
    return min(X[i:] + X[:i] for i in range(len(X)))

def edge_label(i):
    if i < 0:
        return "~%d" % (~i)
    else:
        return str(i)

class ColouredTriangulation(object):
    r"""
    Coloured triangulation.

    A triangulation with a label blue or red set to each color so that there
    is no monochromatic face and no monochromatic vertex.

    EXAMPLES::

        sage: from veerer import *

    Built from an explicit triangulation (in cycle or list form) and a list of colors::

        sage: ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
        [(~2, ~0, ~1), (0, 1, 2)], ['Red', 'Red', 'Blue']

    From a stratum::

        sage: from surface_dynamics import *

        sage: ColouredTriangulation.from_stratum(AbelianStratum(2))
        [(~8, 3, ~1), (~7, 1, 8), (~6, 2, 7), (~5, 0, 6), (~4, 5, ~0), (~3, 4, ~2)], ['Red', 'Red', 'Red', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue']

        sage: ColouredTriangulation.from_stratum(QuadraticStratum({1:4}))
        [(~17, 6, ~2), (~16, 4, 17), ..., (~6, 7, ~5)], ['Red', 'Red', ..., 'Blue']

    From a flipper pseudo-Anosov (TO BE DONE)!
    """
    __slots__ = ['_triangulation', '_colouring']

    def __init__(self, triangulation,  colouring, check=True):
        if isinstance(triangulation, Triangulation):
            self._triangulation = triangulation.copy()
        else:
            self._triangulation = Triangulation(triangulation)

        if isinstance(colouring, str):
            colouring = [RED if c == 'R' else BLUE for c in colouring]

        self._colouring = colouring + colouring[::-1] # A list : edge_indices --> {Red, Blue}

        if check:
            n = self._triangulation.num_edges()
            assert(all(colour in COLOURS for colour in self._colouring))
            assert(len(self._colouring) == 2 * n)
            for col in self._colouring:
                assert col == BLUE or col == RED

            cols = set([RED, BLUE])
            # no monochromatic face
            for F in self._triangulation.faces():
                assert set(self._colouring[e] for e in F) == cols

            # no monochromatic vertex
            for V in self._triangulation.vertices():
                assert set(self._colouring[e] for e in V) == cols

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
            sage: ColouredTriangulation.from_pseudo_anosov(h)
            [(~5, ~3, ~1), (~4, 1, ~0), (~2, 0, 3), (2, 5, 4)],
            ['Blue', 'Blue', 'Blue', 'Red', 'Red', 'Blue']
        """
        f = h.flat_structure()
        n = f.triangulation.zeta

        X = {i.label: e.x for i,e in f.edge_vectors.iteritems()}
        Y = {i.label: e.y for i,e in f.edge_vectors.iteritems()}

        triangles = [[x.label for x in t] for t in f.triangulation]
        colours = [RED if X[e]*Y[e] > 0 else BLUE for e in range(n)]
        return ColouredTriangulation(triangles, colours)

    @classmethod
    def from_stratum(cls, c):
        r"""
        Return a Veering triangulation from either a stratum, a component
        of stratum or a cylinder diagram.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: from veerer import *

            sage: ColouredTriangulation.from_stratum(AbelianStratum(2))
            [(~8, 3, ~1), (~7, 1, 8), (~6, 2, 7), (~5, 0, 6),
             (~4, 5, ~0), (~3, 4, ~2)],
             ['Red', 'Red', 'Red', 'Blue', 'Blue', 'Blue', 'Blue',
              'Blue', 'Blue']

            sage: Q = QuadraticStratum(9,-1)
            sage: ColouredTriangulation.from_stratum(Q.regular_component())
            [(~17, 6, ~5), (~16, ~2, 17), (~15, ~1, 16), (~14, 2, 15),
             (~13, 1, 14), (~12, 0, 13), (~11, 12, ~0), (~10, 11, 3),
             (~9, 10, 4), (~8, 9, ~3), (~7, 8, ~4), (~6, 7, 5)],
             ['Red', 'Red', 'Red', 'Red', 'Red', 'Red', 'Blue', 'Blue',
             'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue',
             'Blue', 'Blue', 'Blue']
            sage: ColouredTriangulation.from_stratum(Q.irregular_component())
            [(~17, 6, ~5), (~16, 1, 17), (~15, ~0, 16), (~14, 0, 15),
             (~13, 14, ~1), (~12, 13, 2), (~11, 12, 3), (~10, 11, 4),
             (~9, 10, 5), (~8, 9, ~2), (~7, 8, ~3), (~6, 7, ~4)],
             ['Red', 'Red', 'Red', 'Red', 'Red', 'Red', 'Blue', 'Blue', 'Blue',
              'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue',
              'Blue']

        Some examples built from cylinder diagram::

            sage: c = QuadraticCylinderDiagram('(0)-(0)')
            sage: ColouredTriangulation.from_stratum(c)
            [(~2, 1, ~0), (~1, 0, 2)], ['Red', 'Blue', 'Blue']

            sage: c = QuadraticCylinderDiagram('(0,0)-(1,1,2,2,3,3)')
            sage: ColouredTriangulation.from_stratum(c)
            [(~11, 4, ~3), (~10, ~0, 11), (~9, 0, 10), (~8, 9, 1), (~7, 8, ~1), (~6, 7, 2), (~5, 6, ~2), (~4, 5, 3)], ['Red', 'Red', 'Red', 'Red', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue']

            sage: c = CylinderDiagram('(0,6,4,5)-(3,6,5) (1,3,2)-(0,1,4,2)')
            sage: CT = ColouredTriangulation.from_stratum(c)
            sage: CT.stratum()
            H_4(6)
        """
        from . import HAS_SAGE
        if not HAS_SAGE:
            raise ValueError("this module only works within SageMath")

        try:
            from surface_dynamics.flat_surfaces.strata import Stratum, StratumComponent
            from surface_dynamics.flat_surfaces.separatrix_diagram import \
                    CylinderDiagram, QuadraticCylinderDiagram
        except ImportError:
            raise ValueError("surface_dynamics is missing")

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

        for bot,top in c.cylinders():
            # len(bot) + len(top)
            l = m + len(top) - 1
            for i in bot:
                assert isinstance(i,int)
                if seen[i]:
                    i = ~i
                else:
                    seen[i] = True
                triangles.append((i,l+1,~l))
                l += 1
            l = m + len(top) - 1
            for i in top:
                assert isinstance(i,int)
                if seen[i]:
                    i = ~i
                else:
                    seen[i] = True
                triangles.append((i,~(l-1),l))
                l -= 1
            # last one got wrong
            i, j, k = triangles.pop()
            triangles.append((i,~(k+len(bot)+len(top)-1),k))

            m += len(bot) + len(top)

        colors = [RED] * c.nseps() + [BLUE] * (2*c.nseps())
        return ColouredTriangulation(triangles, colors)

    @classmethod
    def from_string(cls, s, *args, **kwds):
        r"""
        Deserialization from the string ``s``.

        Such string is typically obtained from :meth:`to_string` or
        :meth:`iso_sig`. Not by calling ``str(triangulation)``.

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation.from_string('BBR_eacdfb')
            sage: T.to_string()
            'BBR_eacdfb'
            sage: T.iso_sig()
            'RBB_adbecf'
        """
        colours, triangles = s.split('_')
        n = len(colours)

        p = even_perm_from_base64_str(triangles, n)
        triangles = even_perm_cycles(p)[0]
        colours = [RED if colour == 'R' else BLUE for colour in colours]
        return ColouredTriangulation(triangles, colours, *args, **kwds)

    from_iso_sig = from_string

    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._triangulation == other._triangulation and self._colouring == other._colouring

    def __ne__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._triangulation != other._triangulation or self._colouring != other._colouring

    def copy(self):
        r"""
        Return a copy of this coloured triangulation

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], "RRB")
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
        T = ColouredTriangulation.__new__(ColouredTriangulation)
        T._triangulation = self._triangulation.copy()
        T._colouring = self._colouring[:]
        return T

    def __str__(self):
        n = self._triangulation.num_edges()
        return str(self._triangulation) + ', ' + str(self._colouring[:n])

    def __repr__(self):
        return str(self)

    def colour(self, e):
        e = int(e)
        return self._colouring[e]

    def angles(self):
        r"""
        Return the list of angles.

        Note that the vertices of the triangulation are labeled. The
        angles are given in the same order.

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(-3, 1, -1), (-2, 0, 2)], [RED, BLUE, BLUE])
            sage: T.angles()
            [2]

            sage: T = ColouredTriangulation([(-6, 2, -2), (-5, 1, 5), (-4, 0, 4), (-3, 3, -1)],
            ....:                           [RED, RED, BLUE, BLUE, BLUE, BLUE])
            sage: T.angles()
            [2, 2]

            sage: T = ColouredTriangulation([(-12, 4, -4), (-11, -1, 11), (-10, 0, 10),
            ....:                            (-9, 9, 1), (-8, 8, -2), (-7, 7, 2),
            ....:                            (-6, 6, -3), (-5, 5, 3)],
            ....:   [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE])
            sage: T.angles()
            [3, 1, 1, 1, 1, 1]
        """
        n = self._triangulation.num_edges()
        angles = []
        seen = [False] * (2 * n)
        vp = self._triangulation.vertex_permutation()
        for e in range(-n, n):
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

        return angles

    def is_abelian(self):
        r"""
        Return whether this coloured triangulation is Abelian.

        EXAMPLES:

            sage: from veerer import *

            sage: T = ColouredTriangulation([(-3, 1, -1), (-2, 0, 2)], [RED, BLUE, BLUE])
            sage: T.is_abelian()
            True

            sage: T = ColouredTriangulation([(-6, 2, -2), (-5, 1, 5), (-4, 0, 4), (-3, 3, -1)],
            ....:       [RED, RED, BLUE, BLUE, BLUE, BLUE])
            sage: T.is_abelian()
            True

            sage: T = ColouredTriangulation([(-12, 4, -4), (-11, -1, 11), (-10, 0, 10),
            ....:      (-9, 9, 1), (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)],
            ....:      [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE])
            sage: T.is_abelian()
            False
        """
        vp = self._triangulation.vertex_permutation()

        # Perform BFS to check that we can coherently orient the
        # imaginary part of each edge
        oris = [None] * (2 * self._triangulation.num_edges())
        oris[0] = True   # = going up
        oris[~0] = False # = going down
        q = [0, ~0]
        while q:
            e0 = q.pop()
            o = oris[e0]

            e = e0
            f = vp[e]
            while f != e0:
                if self._colouring[e] == BLUE and self._colouring[f] == RED:
                    o = not o
                oris[f] = o
                if oris[~f] is None:
                    q.append(~f)
                    oris[~f]= not o
                elif oris[~f] == oris[f]:
                    return False
                e,f = f, vp[f]

        return True

    def stratum(self):
        r"""
        Return the Abelian or quadratic stratum of this coloured triagulation.

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.stratum()
            H_1(0)
        """
        from . import HAS_SAGE
        if not HAS_SAGE:
            raise ValueError('this method is only available inside SageMath. You '
                             'can use the method .angles() to get the list of '
                             'vertex angles.')

        A = self.angles()
        if any(a%2 for a in A) or not self.is_abelian():
            from surface_dynamics import QuadraticStratum
            return QuadraticStratum([(a-2) for a in A])
        else:
            from surface_dynamics import AbelianStratum
            return AbelianStratum([(a-2)/2 for a in A])

    def stratum_dimension(self):
        dim1 = 2*self._triangulation.genus() - 2 + self._triangulation.num_vertices() + (1 if self.is_abelian() else 0)
        dim2 = self.stratum().dimension()
        assert dim1 == dim2
        return dim1

    def colours_about_edge(self, e):
        e = int(e)
        return [self._colouring[f] for f in self._triangulation.square_about_edge(e)]

    def alternating_square(self, e):
        e = int(e)
        colours = self.colours_about_edge(e)
        return all(colours[f] != colours[(f+1) % 4] for f in range(4))

    def is_flippable(self, e):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.is_flippable(0)
            True
            sage: T.is_flippable(1)
            True
            sage: T.is_flippable(2)
            False
        """
        e = int(e)
        return self._triangulation.is_flippable(e) and \
               self.alternating_square(e)

    def flippable_edges(self):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.flippable_edges()
            [0, 1]
        """
        n = self._triangulation.num_edges()
        return [e for e in range(n) if self.is_flippable(e)]

    def is_forward_flippable(self, e):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.is_forward_flippable(0)
            False
            sage: T.is_forward_flippable(1)
            True
            sage: T.is_forward_flippable(2)
            False
        """
        return self._triangulation.is_flippable(e) and \
               self.colours_about_edge(e) == [BLUE, RED, BLUE, RED]

    def forward_flippable_edges(self):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.forward_flippable_edges()
            [1]
        """
        n = self._triangulation.num_edges()
        return [e for e in range(n) if self.is_forward_flippable(e)]

    def is_backward_flippable(self, e):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.is_backward_flippable(0)
            True
            sage: T.is_backward_flippable(1)
            False
            sage: T.is_backward_flippable(2)
            False
        """
        return self._triangulation.is_flippable(e) and \
               self.colours_about_edge(e) == [RED, BLUE, RED, BLUE]

    def backward_flippable_edges(self):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.backward_flippable_edges()
            [0]
        """
        n = self._triangulation.num_edges()
        return [e for e in range(n) if self.is_backward_flippable(e)]

    def mostly_sloped_edges(self, slope):
        if slope == HORIZONTAL:
            return self.forward_flippable_edges()
        elif slope == VERTICAL:
            return self.backward_flippable_edges()
        else:
            raise ValueError

    def good_starts(self):
        r"""
        Start at a RED before a BLUE.
        """
        starts = []
        best = None

        n = self._triangulation.num_edges()
        vp = self._triangulation.vertex_permutation()
        for e in range(-n, n):
            if self._colouring[e] != RED:
                continue
            f = vp[e]
            if self._colouring[f] != BLUE:
                continue

            # build the word
            w = []
            g = f
            while True:

                # go along the v permutation and write a word in BLUE/RED
                # run through all the BLUE
                n = 0
                while self._colouring[g] == BLUE:
                    g = vp[g]
                    n += 1
                w.append(n)
                # run through all the RED
                while self._colouring[g] == RED:
                    g = vp[g]
                    n += 1
                w.append(n)

                if g == f or \
                   (best is not None and len(w) > len(best)):
                    break

            if best is None or \
               (len(w) < len(best) or (len(w) == len(best) and w < best)):
                starts = [e]
                best = w
            elif w == best:
                starts.append(e)

        return starts

    def relabel(self, p):
        r"""
        Relabel the triangulation according to the permutation ``p``.

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(-3, 1, -1), (-2, 0, 2)], [RED, BLUE, BLUE])
            sage: T.relabel([2, -2, -1, 0, 1, -3])
            [(~2, 0, ~1), (~0, 1, 2)], ['Blue', 'Blue', 'Red']
        """
        n = self._triangulation.num_edges()
        if len(p) != 2*n:
            raise ValueError('invalid relabeling permutation')
        colours = [None] * n
        for i in range(-n,n):
            j = int(p[i])
            if ~j != int(p[~i]):
                raise ValueError('invalid relabeling permutation')
            if j >= 0:
                colours[j] = self._colouring[i]

        triangles = [(p[i], p[j], p[k]) for i,j,k in self._triangulation.faces()]

        return ColouredTriangulation(triangles, colours)

    def _relabelling_from(self, start_edge):
        n = self._triangulation.num_edges()
        vp = self._triangulation.vertex_permutation()

        relabelling = [None] * (2*n)
        relabelling[start_edge] = 0
        relabelling[~start_edge] = ~0

        k = 1  # current new label

        to_process = []
        to_process.append(start_edge)
        to_process.append(~start_edge)
        while to_process:
            e0 = to_process.pop()
            e = vp[e0]
            while e != e0:
                if relabelling[e] is None:
                    relabelling[e] = k
                    relabelling[~e] = ~k
                    k += 1
                    to_process.append(~e)
                e = vp[e]

        return relabelling

    def best_relabelling(self):
        n = self._triangulation.num_edges()
        fp = self._triangulation.face_permutation(copy=False)
        best = None

        for start_edge in self.good_starts():
            relabelling = self._relabelling_from(start_edge)

            # relabelled face permutation
            fp_new = [None] * (2*n)
            for e in range(-n, n):
                #fp_new[e] = relabelling[fp[inv_relabelling[e]]]
                fp_new[relabelling[e]] = relabelling[fp[e]]

            # relabelled colouring
            colouring_new = [None] * n
            for e in range(n):
                colouring_new[norm(relabelling[e])] = self._colouring[e]

            T = (colouring_new, fp_new)
            if best is None or T < best:
                best = T
                best_relabelling = relabelling

        return best_relabelling, best

    def automorphisms(self):
        r"""
        Return the list of automorphisms of this coloured triangulation.

        The output is a list of arrays that are signed permutations given
        as lists of length 2n.

        EXAMPLES::

            sage: from veerer import *

        An example with four symmetries in H(1,1)::

            sage: p = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = 'BRBBBRRBBBBR'
            sage: T = ColouredTriangulation(p, cols)
            sage: A = T.automorphisms()
            sage: len(A)
            4
            sage: all(T.relabel(a) == T for a in A)
            True

        A very symmetric example in H(0^5)::

            sage: T = ColouredTriangulation.from_string("RBBBBRBBRBBRBBR_tDnjlwgizdfCacmbBAeyxhvukspqor")
            sage: print(len(T.automorphisms()))
            10
        """

        n = self._triangulation.num_edges()
        fp = self._triangulation.face_permutation(copy=False)

        best = None
        best_relabellings = []

        fp_new = [None] * (2*n)
        colouring_new = [None] * n
        for start_edge in self.good_starts():
            relabelling = self._relabelling_from(start_edge)

            # relabelled face permutation
            for e in range(-n, n):
                #fp_new[e] = relabelling[fp[inv_relabelling[e]]]
                fp_new[relabelling[e]] = relabelling[fp[e]]

            # relabelled colouring
            for e in range(n):
                colouring_new[norm(relabelling[e])] = self._colouring[e]

            T = (colouring_new, fp_new)
            if best is None or T == best:
                best_relabellings.append(relabelling)
                if best is None:
                    best = (colouring_new[:], fp_new[:])
            elif T < best:
                del best_relabellings[:]
                best_relabellings.append(relabelling)
                best = (colouring_new[:], fp_new[:])

        p0 = even_perm_invert(best_relabellings[0])
        return [even_perm_compose(p, p0) for p in best_relabellings]

    def to_string(self):
        r"""
        Serialization to string.

        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: T = ColouredTriangulation.from_stratum(QuadraticStratum({1:20}))
            sage: s = T.to_string()
            sage: s
            'RRR...a3aza'
            sage: TT = ColouredTriangulation.from_string(s)
            sage: T == TT
            True
        """
        n = self._triangulation.num_edges()
        colours = ''.join('R' if self._colouring[i] == RED else 'B' for i in range(n))
        perm = even_perm_base64_str(self._triangulation.face_permutation())
        return colours + '_' + perm

    def iso_sig(self):
        r"""
        Return a canonical string ("isomorphism signature").

        EXAMPLES::

            sage: from veerer import *

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = ColouredTriangulation(t, cols)
            sage: T.iso_sig()
            'RBBBBBBBRRRB_twbvcudxfhjmkigesrqponla'

        If we relabel the triangulation, the isomorphic signature does not change::

            sage: from veerer.even_permutation import signed_permutation_random
            sage: p = signed_permutation_random(12)
            sage: T.relabel(p).iso_sig()
            'RBBBBBBBRRRB_twbvcudxfhjmkigesrqponla'

        An isomorphic triangulation can be reconstructed from the isomorphic
        signature via::

            sage: s = T.iso_sig()
            sage: T2 = ColouredTriangulation.from_string(s)
            sage: T == T2
            False
            sage: T.is_isomorphic_to(T2)
            True

        TESTS::

            sage: T = ColouredTriangulation.from_string('RBBBBRBRB_aqhbprdgcfeonmjkil')
            sage: p = [2, 3, 5, 0, -5, -8, 6, -9, 1, -2, 8, -7, 7, 4, -1, -6, -4, -3]
            sage: T.iso_sig() == T.relabel(p).iso_sig()
            True
        """
        n = self._triangulation.num_edges()
        _, (colours, fp) = self.best_relabelling()

        char_colours = ''.join('R' if colours[i] == RED else 'B' for i in range(n))
        char_edges = even_perm_base64_str(fp)

        return char_colours + '_' + char_edges

    def canonical(self):
        r"""
        Return a canonically labeled isomorphic coloured triangulation.
        """
        n = self._triangulation.num_edges()
        s = self.iso_sig()
        return ColouredTriangulation.from_iso_sig(self.iso_sig())

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

    def flip(self, i, col):
        r"""
        Flip an edge inplace.

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T.flip(1, RED); T
            [(~2, 1, 0), (~1, ~0, 2)], ['Red', 'Red', 'Blue']
            sage: T.flip(0, RED); T
            [(~2, ~0, ~1), (0, 1, 2)], ['Red', 'Red', 'Blue']

            sage: T.flip(1, BLUE); T
            [(~2, 1, 0), (~1, ~0, 2)], ['Red', 'Blue', 'Blue']
            sage: T.flip(2, BLUE); T
            [(~2, 0, ~1), (~0, 1, 2)], ['Red', 'Blue', 'Blue']
        """
        i = int(i)
        assert(self.is_flippable(i))

        self._triangulation.flip(i)
        self._colouring[i] = self._colouring[~i] = col

    def back_flip(self, i, col):
        r"""
        Flip backward an edge in place

        EXAMPLES::

            sage: from veerer import *

            sage: T0 = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: T = T0.copy()
            sage: T.flip(1, RED); T.flip(0, RED)
            sage: T.back_flip(0, RED); T.back_flip(1, RED)
            sage: T == T0
            True

            sage: T.flip(1, BLUE); T.flip(2, BLUE)
            sage: T.back_flip(2, BLUE); T.back_flip(1, RED)
            sage: T == T0
            True
        """
        i = int(i)
        assert(self.is_flippable(i))

        self._triangulation.back_flip(i)
        self._colouring[i] = self._colouring[~i] = col

    def _set_train_track_constraints(self, insert, x, slope, integral, allow_degenerations):
        r"""
        Sets the equation and inequations for train tracks.

        INPUT:

        - ``insert`` - a function to be called for each equation or inequation

        - ``x`` - variable factory (the variable for edge ``e`` is constructed
          via ``x[e]``

        - ``slope`` - the slope of the train-track

        - ``integral`` - boolean - whether to set lower bounds to 0 or 1

        - ``allow_degenerations`` - boolean - when ``integral`` is also set to ``True``
          then allows certain lengths to be 0

        EXAMPLES::

            sage: from veerer import *

            sage: T  = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, False, False)
            sage: l
            [x0 >= 0, x1 >= 0, x2 >= 0, x0 == x1 + x2, x0 == x1 + x2]

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, True, False)
            sage: l
            [x0 >= 1, x1 >= 1, x2 >= 1, x0 == x1 + x2, x0 == x1 + x2]

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, True, True)
            sage: l
            [x0 >= 1, x1 >= 0, x2 >= 1, x0 == x1 + x2, x0 == x1 + x2]

        This can also be used to check that a given vector satisfies the train-track
        equations::

            sage: def check(x):
            ....:     if not x:
            ....:         raise AssertionError
            sage: T._set_train_track_constraints(check, [2,1,1], HORIZONTAL, False, False)
            sage: T._set_train_track_constraints(check, [1,1,1], HORIZONTAL, False, False)
            Traceback (most recent call last):
            ...
            AssertionError
        """
        if slope == VERTICAL:
            POS = BLUE
            NEG = RED
        elif slope == HORIZONTAL:
            POS = RED
            NEG = BLUE
        else:
            raise ValueError('bad slope parameter')

        n = self._triangulation.num_edges()

        # non-negativity
        for e in range(n):
            if not integral or \
                (allow_degenerations and \
                 ((slope == HORIZONTAL and self.is_forward_flippable(e)) or \
                 (slope == VERTICAL and self.is_backward_flippable(e)))):
                insert(x[e] >= 0)
            else:
                insert(x[e] >= 1)

        # switches relation
        for (i,j,k) in self._triangulation.faces():
            # find the large edge
            # if horizontal, this is the one opposite to the RED/BLUE transition
            # if vertical this is the one opposite to the BLUE/RED transition
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
                raise ValueError('can not determine the big edge on triangle (%s, %s, %s' %
                                 (edge_label(i), edge_label(j), edge_label(k)))

            insert(x[norm(l)] == x[norm(s1)] + x[norm(s2)])

    def train_track_polytope(self, slope=VERTICAL):
        r"""
        Return the polytope determined by the constraints.

        INPUT:

        - ``slope`` - the slope for the train track

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2),(-1,-2,-3)], [RED, RED, BLUE])
            sage: P = T.train_track_polytope(VERTICAL)
            sage: P
            A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 1 point, 2 rays
            sage: P.generators()
            Generator_System {point(0/1, 0/1, 0/1), ray(1, 1, 0), ray(0, 1, 1)}
        """
        try:
            import ppl
        except ImportError:
            raise ImportError('pplpy is not installed on your computer. Follow '
                              'the instructions at https://gitlab.com/videlec/pplpy')

        cs = ppl.Constraint_System()
        n = self._triangulation.num_edges()
        variables = [ppl.Variable(e) for e in range(n)]
        self._set_train_track_constraints(cs.insert, variables, slope, False, False)
        return ppl.C_Polyhedron(cs)

    def train_track_min_solution(self, slope=VERTICAL, allow_degenerations=False):
        r"""
        Return the minimal integral point satisfying the constraints.

        INPUT:

        - ``slope`` - the slope of the train track

        - ``allow_degenerations`` - boolean - if ``True`` then allow certain
          degenerations to occur.

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(0,1,2),(-1,-2,-3)], [RED, RED, BLUE])
            sage: T.train_track_min_solution(VERTICAL)
            [1.0, 2.0, 1.0]
            sage: T.train_track_min_solution(VERTICAL, allow_degenerations=True)
            [0.0, 1.0, 1.0]

            sage: T.train_track_min_solution(HORIZONTAL)
            [2.0, 1.0, 1.0]
            sage: T.train_track_min_solution(HORIZONTAL, allow_degenerations=True)
            [1.0, 0.0, 1.0]
        """
        from sage.numerical.mip import MixedIntegerLinearProgram

        n = self._triangulation.num_edges()
        M = MixedIntegerLinearProgram(maximization=False)
        x = M.new_variable(integer=True)
        M.set_objective(M.sum(x[e] for e in range(n))) # try to minimize length
        self._set_train_track_constraints(M.add_constraint, x, slope, True, allow_degenerations)
        M.solve()
        x = M.get_values(x)
        return [x[e] for e in range(n)]

    def geometric_polytope(self):
        raise NotImplementedError
        # 1. train-track conditions
        P = self.train_track_polytope(HORIZONTAL)
        P.concatenate_assign(self.train_track_polytope(VERTICAL))

        # 2. flip conditions
        for e in self.flippable_edges():
            raise NotImplementedError
            corner1 = self._triangulation.corner_lookup[i]
            corner2 = self._triangulation.corner_lookup[~i]

        return P

    def geometric_polytope(self, normalise=False):
        raise NotImplementedError

        # Uses PyParma.
        V = self.train_track_matrix(VERTICAL)
        H = self.train_track_matrix(HORIZONTAL)
        G = self.geometric_matrix()

        A = [[0] + [0] * i + [1] + [0] * (2*self.zeta - i - 1) for i in range(2*self.zeta)]
        for row in V:
            A.append([0] + [i for i in row] + [0] * self.zeta)
            A.append([-0] + [-i for i in row] + [-0] * self.zeta)
        for row in H:  ### 2017-12-20 - SS - Was V for some reason???
            A.append([0] + [0] * self.zeta + [i for i in row])
            A.append([-0] + [-0] * self.zeta + [-i for i in row])
        for row in G:
            A.append([0] + row)
        if normalise:
            A.append([-1] + [1] * self.zeta * 2)
            A.append([1] + [-1] * self.zeta * 2)

        return Polyhedron(ieqs=A)

    def geometric_dimension(self):
        raise NotImplementedError
        return self.geometric_polytope().dim()

    def _flat_structure_from_train_track_lengths(self, VH, VV, base_ring=None):
        r"""
        Return a flat structure from two vectors ``VH`` and ``VV``
        satisfying the train track equations.
        """
        from sage.modules.free_module import VectorSpace

        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ

        n = self._triangulation.num_edges()
        assert len(VH) == len(VV) == n
        assert all(x >=0 for x in VH)
        assert all(x >= 0 for x in VV)

        V = VectorSpace(base_ring, 2)
        vectors = [V((x, y if self._colouring[i] == RED else -y)) for \
                   i,(x,y) in enumerate(zip(VV, VH))]
        vectors.extend(vectors[::-1])

        # get correct signs for each triangle
        for i,j,k in self._triangulation:
            if det2(vectors[i], vectors[j]) < 0:
                vectors[j] = -vectors[j]
            if det2(vectors[j], vectors[k]) < 0:
                vectors[k] = -vectors[k]
            if vectors[i] + vectors[j] + vectors[k]:
                raise RuntimeError('bad vectors:\n vec[%s] = %s\n vec[%s] = %s\n vec[%s] = %s' \
                                   % (edge_label(i), vectors[i], edge_label(j), vectors[j], edge_label(k), vectors[k]))

            if det2(vectors[k], vectors[i]) < 0:
                raise RuntimeError

        from .layout import FlatTriangulation
        return FlatTriangulation(self._triangulation, vectors)

    def flat_structure_middle(self):
        r"""
        Return a flat structure with this Veering triangulation.

        Note that this triangulation must be core. The point is chosen
        by taking the interior point of the polytope obtained by
        summing each ray.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: from veerer import *

            sage: Q = QuadraticStratum({1:4, -1:4})
            sage: CT = ColouredTriangulation.from_stratum(Q)
            sage: F = CT.flat_structure_middle()
            sage: F
            Flat Triangulation made of 16 triangles
        """
        n = self._triangulation.num_edges()

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
            sage: CT = ColouredTriangulation.from_stratum(Q)
            sage: CT.flat_structure_min()
            Flat Triangulation made of 16 triangles

        By allowing degenerations you can get a simpler solution but
        with some of the edges horizontal or vertical::

            sage: CT.flat_structure_min(True)
            Flat Triangulation made of 16 triangles
        """
        VH = self.train_track_min_solution(HORIZONTAL, allow_degenerations=allow_degenerations)
        VV = self.train_track_min_solution(VERTICAL, allow_degenerations=allow_degenerations)

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

            sage: T = ColouredTriangulation(triangles, colours)
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
            from . import HAS_SAGE
            if not HAS_SAGE:
                raise ValueError('the option LP is only available in SageMath')
            from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

            n = self._triangulation.num_edges()
            for slope in (HORIZONTAL, VERTICAL):
                M = MixedIntegerLinearProgram(maximization=False)
                x = M.new_variable()
                M.set_objective(M.sum(x[e] for e in range(n))) # try to minimize length
                self._set_train_track_constraints(M.add_constraint, x, slope, True, False)

                try:
                    M.solve()
                except MIPSolverException:
                    return False

            return True

        else:
            raise ValueError("method must be either 'polytope' or 'LP'")

    def is_geometric(self, d=None):
        raise NotImplementedError
        if d is None: d = self.stratum_dimension()
        return all(self.train_track_dimension(slope) == d for slope in SLOPES) and self.geometric_dimension() == 2*d

    def edge_has_curve(self, e, verbose=False):
        r"""
        INPUT:

        - e - edge label

        OUTPUT: boolean - whether there is a curve which crosses the associated
        train-track in the correct direction.

        EXAMPLES::

            sage: from veerer import *

            sage: t = [(-6, -4, -5), (-3, -1, 3), (-2, 0, 4), (1, 5, 2)]
            sage: c = [RED, RED, BLUE, RED, BLUE, RED]
            sage: T0 = ColouredTriangulation(t, c)
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
        # TODO: we should not allow more than one pair (i, ~i) to be
        # both seen

        # TODO: we want the._colouring to also take care of negative edges
        # (bad alternative: take "norm" function from flipper)
        colouring = self._colouring

        if verbose:
            print('[edge_has_curve] checking edge %s with color %s' % (edge_label(e), colouring[e]))

        a,b,c,d = self._triangulation.square_about_edge(e)
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
            end = ~d
            if verbose:
                print('[edge_has_curve] try to find path from b=%s to ~d=%s' %
                         (edge_label(start), edge_label(end)))
        else:
            start = a
            end = ~c
            if verbose:
                print('[edge_has_curve] try to find path from a=%s to ~c=%s' %
                         (edge_label(start), edge_label(end)))

        if start == end:
            return True

        n = self._triangulation.num_edges()
        fp = self._triangulation.face_permutation()
        seen = [False] * (2 * n)
        seen[start] = True
        q = deque()
        q.append(start)
        while q:
            f = q.popleft()
            if verbose:
                print('[edge_has_curve] crossing %s' % edge_label(f))

            # NOTE: would be much nicer to have a face permutation!
            # here we set r and s so that the triangle is (r, s, ~f)
            r = fp[~f]
            s = fp[r]
            if verbose:
                print('[edge_has_curve] switch with r=%s s=%s' % (edge_label(r), edge_label(s)))
            if not seen[r] and not (colouring[r] == POS and colouring[~f] == NEG):
                if r == end:
                    if verbose:
                        print('[edge_has_curve] done at %s' % edge_label(r))
                    return True
                seen[r] = True
                q.append(r)
                if verbose:
                    print('[edge_has_curve] adding %s on top of the queue' % edge_label(r))
            if not seen[s] and not (colouring[s] == NEG and colouring[~f] == POS):
                if s == end:
                    if verbose:
                        print('[edge_has_curve] done at %s' % edge_label(s))
                    return True
                seen[s] = True
                q.append(s)
                if verbose:
                    print('[edge_has_curve] adding %s on top of the queue' % edge_label(s))

        return False

    def _check_edge_has_curve(self):
        r"""
        Check the function ``edge_has_curve``

        EXAMPLES::

            sage: from veerer import *
            sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
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

def ngon(n):
    n = int(n)
    assert(n > 4 and n % 2 == 0)
    m = (n - 2) // 2
    T = [(i, 2*m+i, ~(i+1)) for i in range(m)] + \
        [(~0, ~(2*m), ~(m+1))] + \
        [(m+i, ~(2*m+i), ~(m+i+1)) for i in range(1, m-1)] + \
        [(2*m-1, ~(3*m-1), m)]
    colouring = [RED] * (2*m) + [BLUE] * m
    return ColouredTriangulation(T, colouring)


