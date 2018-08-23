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
from .triangulation import Triangulation

def dot(A, B):
    return sum(a*b for a, b in zip(A, B))

def best_rotation(X):
    return min(X[i:] + X[:i] for i in range(len(X)))

# To be deleted.

def edge_label(e):
    raise ValueError

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
        [(~17, 6, ~3), (~16, 4, 17), ..., (~6, 7, ~2)], ['Red', 'Red', ..., 'Blue']

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
        
        n = self._triangulation.num_half_edges()
        # A list : edge_indices --> {Red, Blue}
        self._colouring = [colouring[self._triangulation._norm(i)] for i in range(n)]

        if check:
            assert(all(colour in COLOURS for colour in self._colouring))
            for col in self._colouring:
                assert col == BLUE or col == RED

            cols = set([RED, BLUE])
            # no monochromatic face
            for f in self._triangulation.faces():
                assert set(self._colouring[e] for e in f) == cols

            # no monochromatic vertex
            for v in self._triangulation.vertices():
                assert set(self._colouring[e] for e in v) == cols

    def _edge_rep(self, e):
        return self._triangulation._edge_rep(e)

    def _norm(self, e):
        return self._triangulation._norm(e)

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
        FS = h.flat_structure()
        n = FS.triangulation.zeta

        X = {i.label: e.x for i,e in f.edge_vectors.iteritems()}
        Y = {i.label: e.y for i,e in f.edge_vectors.iteritems()}

        triangles = [[x.label for x in t] for t in FS.triangulation]
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
            [(~8, 3, ~1), (~7, 1, 8), ..., (~3, 4, ~2)], ['Red', 'Red', ..., 'Blue',
              'Blue', 'Blue']

            sage: Q = QuadraticStratum(9,-1)
            sage: CTreg = ColouredTriangulation.from_stratum(Q.regular_component())
            sage: CTreg.stratum()
            Q_3(9, -1)

            sage: CTirr = ColouredTriangulation.from_stratum(Q.irregular_component())
            sage: CTirr.stratum()
            Q_3(9, -1)

        Some examples built from cylinder diagram::

            sage: c = QuadraticCylinderDiagram('(0)-(0)')
            sage: ColouredTriangulation.from_stratum(c)
            [(~2, 1, ~0), (~1, 0, 2)], ['Red', 'Blue', 'Blue']

            sage: c = QuadraticCylinderDiagram('(0,0)-(1,1,2,2,3,3)')
            sage: ColouredTriangulation.from_stratum(c)
            [(~11, 4, ~3), (~10, ~0, 11), ..., (~4, 5, 3)], ['Red', 'Red', ..., 'Blue']

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
        Return the list of angles (divided by \pi). 

        Note that the vertices of the triangulation are labeled. The
        angles are given in the same order.  

        # ??? Confused - there is no labelling yet. 

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
        n = self._triangulation.num_half_edges()
        angles = []
        seen = [False] * n
        vp = self._triangulation.vertex_permutation()
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

        return angles + [1]*self.num_folded_edges()

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
        if self._triangulation.num_folded_edges() > 0:
            return False
        
        ep = self._triangulation.edge_permutation()
        vp = self._triangulation.vertex_permutation()

        # Perform BFS to check that we can coherently orient the
        # imaginary part of each edge
        oris = [None] * self._triangulation.num_half_edges()
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
        dim1 = 2*self._triangulation.genus() - 2 \
               + self._triangulation.num_vertices() \
               + self._triangulation.num_folded_edges() \ # poles!
               + (1 if self.is_abelian() else 0)
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
        best_word = None

        n = self._triangulation.num_half_edges()
        vp = self._triangulation.vertex_permutation(copy=False)
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
        fp = self._triangulation.face_permutation(copy=False)
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

    def relabel(self, p):
        r"""
        Relabel the triangulation according to the permutation ``p``.

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation([(-3, 1, -1), (-2, 0, 2)], [RED, BLUE, BLUE])
            sage: T.relabel([2, -2, -1, 0, 1, -3])
            [(~2, 0, ~1), (~0, 1, 2)], ['Blue', 'Blue', 'Red']
        """
        n = self._triangulation.num_half_edges()
        ep = self._triangulation.edge_permutation()
        if len(p) != n:
            raise ValueError('invalid relabeling permutation')
        colours = [None] * n
        for i in range(n):
            j = int(p[i])
            if ep[j] != int(p[ep[i]]):
                raise ValueError('invalid relabeling permutation')
            if j >= 0:
                colours[j] = self._colouring[i]

        triangles = [(p[i], p[j], p[k]) for i,j,k in self._triangulation.faces()]

        return ColouredTriangulation(triangles, colours)

    def _relabelling_from(self, start_edge):
        n = self._triangulation.num_half_edges()
        ep = self._triangulation.edge_permutation()
        vp = self._triangulation.vertex_permutation()
        
        k = 0 # current available label at the front. 
        m = n - 1 # current available label at the back. 
        relabelling = [None] * n
        relabelling[start_edge] = 0
        k = k + 1
        if ep[start_edge] != start_edge:
            relabelling[ep[start_edge]] = m
            m = m - 1

        to_process = []
        to_process.append(start_edge)
        if ep[start_edge] != start_edge:
            to_process.append(ep[start_edge])
        while to_process:
            e0 = to_process.pop()
            e = vp[e0]
            while e != e0:
                if relabelling[e] is None:
                    relabelling[e] = k
                    k = k + 1
                    if ep[e] != e:
                        relabelling[ep[e]] = m
                        m = m - 1
                    to_process.append(ep[e])
                e = vp[e]

        return relabelling

    def best_relabelling(self):
        n = self._triangulation.num_half_edges()
        fp = self._triangulation.face_permutation(copy=False)
        best = None

        for start_edge in self.good_starts():
            relabelling = self._relabelling_from(start_edge)

            # relabelled face permutation
            fp_new = [None] * n
            for e in range(n):
                #fp_new[e] = relabelling[fp[inv_relabelling[e]]]
                fp_new[relabelling[e]] = relabelling[fp[e]]

            # relabelled colouring
            colouring_new = [None] * n
            for e in range(n):
#                colouring_new[self._norm(relabelling[e])] = self._colouring[e]
                colouring_new[relabelling[e]] = self._colouring[e]

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
        n = self._triangulation.num_half_edges()
        fp = self._triangulation.face_permutation(copy=False)

        best = None
        best_relabellings = []

        fp_new = [None] * n
        colouring_new = [None] * n
        for start_edge in self.good_starts():
            relabelling = self._relabelling_from(start_edge)

            # relabelled face permutation
            for e in range(n):
                #fp_new[e] = relabelling[fp[inv_relabelling[e]]]
                fp_new[relabelling[e]] = relabelling[fp[e]]

            # relabelled colouring
            for e in range(n):
                # removed norm... 
                colouring_new[relabelling[e]] = self._colouring[e]

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

    def rotate(self):
        r"""
        Apply the pi/2 rotation.

        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: T0 = ColouredTriangulation.from_stratum(QuadraticStratum(1,1,1,1))
            sage: T = T0.copy()
            sage: T.rotate()
            sage: T.conjugate()
            sage: S = T0.copy()
            sage: S.conjugate()
            sage: S.rotate()
            sage: S == T
            True
        """
        self._colouring = [RED if x == BLUE else BLUE for x in self._colouring]

    def conjugate(self):
        r"""
        Conjugate this triangulation.

        EXAMPLES::

            sage: from veerer import *

            sage: faces = "(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)"
            sage: cols = [BLUE,RED,RED,BLUE,RED,RED]
            sage: T = ColouredTriangulation(faces, cols)
            sage: T.conjugate()
            sage: T
            [(~5, ~4, ~3), (~2, ~1, ~0), (0, 2, 4), (1, 3, 5)], ['Red', 'Blue', 'Blue', 'Red', 'Blue', 'Blue']
        """
        self._colouring = [RED if x == BLUE else BLUE for x in self._colouring]
        self._triangulation.conjugate()

    # TODO: finish this!!
    def automorphism_quotient(self, aut):
        r"""
        Return the canonical triangulation of the quotient.

        (vertex fixed, edge fixed, faces fixed)

        EXAMPLES::

            sage: from veerer import *

            sage: p = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = 'BRBBBRRBBBBR'
            sage: T = ColouredTriangulation(p, cols)
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
        verts, edge_to_vert = even_perm_cycles(self._vp)
        for i,v in enumerate(verts):
            if edge_to_vert[aut[v[0]]] == i:
                # deg does not change
                nb_edges += 1
            else:
                # deg is wrapped around
                pass

        # faces
        faces, edge_to_face = even_perm_cycles(self._fp)
        for i,f in enumerate(faces):
            nb_faces += edge_to_face[aut[f[0]]] == i

        return (nb_verts, nb_edges, nb_faces)

    def to_string(self):
        r"""
        Serialization to string.

        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: T = ColouredTriangulation.from_stratum(QuadraticStratum({1:20}))
            sage: s = T.to_string()
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

    def flip(self, e, col):
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
        e = int(e)
        assert(self.is_flippable(e))
        E = self._triangulation.edge_permutation()[e]

        self._triangulation.flip(e)
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

            sage: T = ColouredTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.cylinders(RED)
            [[0, -2]]
            sage: T.cylinders(BLUE)
            []

            sage: T = ColouredTriangulation("(0,1,2)(~0,~1,~2)", "BRB")
            sage: T.cylinders(BLUE)
            [[0, -3]]
            sage: T.cylinders(RED)
            []

        An example with both blue and red cylinders in Q(1,1,1,1)::

            sage: T = ColouredTriangulation.from_string('BBRBBBRRBRBRRBBRRB_uBFyjzCdocvwqsemEbrgtxlAJkHDnGihfaIp')
            sage: T.cylinders(BLUE)
            [[3, 8, 4, -11]]
            sage: T.cylinders(RED)
            [[9, 15]]

        TESTS::

            sage: from veerer import *

            sage: T = ColouredTriangulation.from_string('RBBBRBBBBRRR_vwkesxgabijfdcuthrqpmnlo')
            sage: T.cylinders(BLUE)
            [[1, 2], [6, 5]]
            sage: T.cylinders(RED)
            []
        """
        n = self._triangulation.num_half_edges()
        fp = self._triangulation.face_permutation(copy=False)
        ep = self._triangulation.edge_permutation(copy=False)
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
                seen[a] = seen[ep[a]] = seen[b] = seen[ep[b]] = seen[c] = seen[ep[c]] = True
                continue

            # find an edge to cross
            if cols[a] == col:
                door = a
            else:
                door = b
            door0 = door

            cc = []  # cycle to be stored
            cyl = True
            while not seen[door]:
                seen[door] = seen[ep[door]] = True
                cc.append(door)

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
            if cyl and door == door0:
                cylinders.append(cc)

        return cylinders

    def is_cylindrical(self, col=None):
        r"""
        Return whether this veering triangulation is cylindrical.

        A Veering triangulation is blue cylindrical (resp red cylindrical) if
        all its triangles have two blue edges (resp red).

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.is_cylindrical()
            True
        """
        if col is None:
            return self.is_cylindrical(BLUE) or self.is_cylindrical(RED)

        n = self._triangulation.num_half_edges()
        fp = self._triangulation.face_permutation(copy=False)
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

    def back_flip(self, e, col):
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
        e = int(e)
        assert(self.is_flippable(e))
        E = self._triangulation.edge_permutation()[e]
        
        self._triangulation.back_flip(i)
        self._colouring[e] = self._colouring[E] = col

    def _set_train_track_constraints(self, insert, x, slope, low_bound, allow_degenerations):
        r"""
        Sets the equation and inequations for train tracks.

        INPUT:

        - ``insert`` - a function to be called for each equation or inequation

        - ``x`` - variable factory (the variable for edge ``e`` is constructed
          via ``x[e]``

        - ``slope`` - the slope of the train-track

        - ``low_bound`` - boolean - whether to set lower bounds to 0 or 1

        - ``allow_degenerations`` - boolean - allow to ignore the lower bound
          in appropriate situations

        EXAMPLES::

            sage: from veerer import *

            sage: T  = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, 0, False)
            sage: l
            [x0 >= 0, x1 >= 0, x2 >= 0, x0 == x1 + x2, x0 == x1 + x2]

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, 1, False)
            sage: l
            [x0 >= 1, x1 >= 1, x2 >= 1, x0 == x1 + x2, x0 == x1 + x2]

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, 3, False)
            sage: l
            [x0 >= 3, x1 >= 3, x2 >= 3, x0 == x1 + x2, x0 == x1 + x2]

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: T._set_train_track_constraints(l.append, x, HORIZONTAL, 2, True)
            sage: l
            [x0 >= 2, x1 >= 0, x2 >= 2, x0 == x1 + x2, x0 == x1 + x2]

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
            if not low_bound or \
                (allow_degenerations and \
                 ((slope == HORIZONTAL and self.is_forward_flippable(e)) or \
                 (slope == VERTICAL and self.is_backward_flippable(e)))):
                insert(x[e] >= 0)
            else:
                insert(x[e] >= low_bound)

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
                raise ValueError('can not determine the big edge on triangle (%s, %s, %s)' %
                                 (self_.edge_rep(i), self._edge_rep(j), self._edge_rep(k)))

            insert(x[self._norm(l)] == x[self._norm(s1)] + x[self._norm(s2)])

    def _set_geometric_constraints(self, insert, x, y, hw_bound=0):
        r"""
        Set the geometric constraints.

        INPUT:

        - ``insert`` - function

        - ``x``, ``y`` - variables

        - ``hw_bound`` - a nonegative number

        EXAMPLES::

            sage: from veerer import *

            sage: T  = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])

            sage: l = []
            sage: x = [SR.var("x0"), SR.var("x1"), SR.var("x2")]
            sage: y = [SR.var("y0"), SR.var("y1"), SR.var("y2")]
            sage: T._set_geometric_constraints(l.append, x, y)
            sage: l
            [x1 <= y0 + y2, y0 <= x1 + x2]
        """
        for e in self.forward_flippable_edges():
            a, b, c, d = self._triangulation.square_about_edge(e)
            insert(x[self._norm(e)] <= y[self._norm(a)] + y[self._norm(d)] - hw_bound)
        for e in self.backward_flippable_edges():
            a, b, c, d = self._triangulation.square_about_edge(e)
            insert(y[self._norm(e)] <= x[self._norm(a)] + x[self._norm(d)] - hw_bound)

    def train_track_polytope(self, slope=VERTICAL, low_bound=0):
        r"""
        Return the polytope determined by the constraints.

        INPUT:

        - ``slope`` - the slope for the train track

        - ``low_bound`` - integer - optional lower bound for the lengths
          (default to 0)

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
        self._set_train_track_constraints(cs.insert, variables, slope, low_bound, False)
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
            [1, 2, 1]
            sage: T.train_track_min_solution(VERTICAL, allow_degenerations=True)
            [0, 1, 1]

            sage: T.train_track_min_solution(HORIZONTAL)
            [2, 1, 1]
            sage: T.train_track_min_solution(HORIZONTAL, allow_degenerations=True)
            [1, 0, 1]
        """
        from sage.numerical.mip import MixedIntegerLinearProgram

        n = self._triangulation.num_edges()
        M = MixedIntegerLinearProgram(solver='PPL', maximization=False)
        x = M.new_variable()
        M.set_objective(M.sum(x[e] for e in range(n))) # try to minimize length
        self._set_train_track_constraints(M.add_constraint, x, slope, 1, allow_degenerations)
        M.solve()
        x = M.get_values(x)
        return [x[e] for e in range(n)]

    def geometric_polytope(self, x_low_bound=0, y_low_bound=0, hw_bound=0):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.geometric_polytope()
            A 4-dimensional polyhedron in QQ^6 defined as the convex hull of 1 point, 7 rays

            sage: s = 'RRRRBBBBBBBB_uvwfgkrpedchjaonmtlixbsq'
            sage: T = ColouredTriangulation.from_string(s)
            sage: T.geometric_polytope()
            A 8-dimensional polyhedron in QQ^24 defined as the convex hull of 1 point, 20 rays

        TESTS::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: T = ColouredTriangulation.from_stratum(AbelianStratum(1,1))
            sage: for _ in range(100):
            ....:     T.random_forward_flip()
            ....:     test1 = T.geometric_polytope().affine_dimension() == 10
            ....:     test2 = not T.geometric_polytope(1,1,1).is_empty()
            ....:     assert test1 == test2, T.to_string()

            sage: T = ColouredTriangulation.from_stratum(QuadraticStratum(1,1,1,1))
            sage: for _ in range(100):
            ....:     T.random_forward_flip()
            ....:     test1 = T.geometric_polytope().affine_dimension() == 12
            ....:     test2 = not T.geometric_polytope(1,1,1).is_empty()
            ....:     assert test1 == test2, T.to_string()
        """
        try:
            import ppl
        except ImportError:
            raise ImportError('pplpy is not installed on your computer. Follow '
                              'the instructions at https://gitlab.com/videlec/pplpy')

        n = self._triangulation.num_edges()

        # 1. train-track conditions
        P = self.train_track_polytope(VERTICAL, low_bound=x_low_bound)
        P.concatenate_assign(self.train_track_polytope(HORIZONTAL, low_bound=y_low_bound))

        x = [ppl.Variable(e) for e in range(n)]
        y = [ppl.Variable(n+e) for e in range(n)]

        self._set_geometric_constraints(P.add_constraint, x, y, hw_bound=hw_bound)

        return P

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
                raise RuntimeError('bad vectors for %s:\n vec[%s] = %s\n vec[%s] = %s\n vec[%s] = %s' \
                                   % (self.to_string(), self._edge_rep(i), vectors[i], self._edge_rep(j), \
                                      vectors[j], self._edge_rep(k), vectors[k]))

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

        A bad example::

            sage: T = ColouredTriangulation.from_string("RBRRBBRBR_gjbpaqechdfonmlkir")
            sage: T.flat_structure_min()
            Flat Triangulation made of 6 triangles
        """
        VH = self.train_track_min_solution(HORIZONTAL, allow_degenerations=allow_degenerations)
        VV = self.train_track_min_solution(VERTICAL, allow_degenerations=allow_degenerations)

        return self._flat_structure_from_train_track_lengths(VH, VV)

    def flat_structure_geometric_middle(self):
        r"""
        Return a geometric flat structure obtained by averaging the
        vertices of the geometric polytope.

        EXAMPLES::

            sage: from veerer import *

            sage: T = ColouredTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.flat_structure_geometric_middle()
            Flat Triangulation made of 2 triangles
        """
        n = self._triangulation.num_edges()

        P = self.geometric_polytope()
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
            try:
                import ppl
            except ImportError:
                raise ImportError('pplpy is not installed on your computer. Follow '
                                  'the instructions at https://gitlab.com/videlec/pplpy')

            n = self._triangulation.num_edges()
            x = [ppl.Variable(i) for i in range(n)]
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

            sage: s = 'RRRRBBBBBBBB_uvwfjilpedcbqkonmtsrhaxg'
            sage: T = ColouredTriangulation.from_string(s)
            sage: T.is_geometric('polytope')
            False
            sage: T.is_geometric('LP')
            False

            sage: s = 'RRRRBBBBBBBB_uvwfgkrpedchjaonmtlixbsq'
            sage: T = ColouredTriangulation.from_string(s)
            sage: T.is_geometric('polytope')
            True
            sage: T.is_geometric('LP')
            True

        TESTS::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: T = ColouredTriangulation.from_stratum(AbelianStratum(1,1))
            sage: for _ in range(100):
            ....:     T.random_forward_flip()
            ....:     test1 = T.is_geometric('LP')
            ....:     test2 = T.is_geometric('polytope')
            ....:     test3 = T.is_geometric('polytope2')
            ....:     assert test1 == test2 == test3, T.to_string()

            sage: T = ColouredTriangulation.from_stratum(QuadraticStratum(1,1,1,1))
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
            return not self.geometric_polytope(1,1,1).is_empty()
        elif method == 'LP':
            raise RuntimeError
            try:
                import ppl
            except ImportError:
                raise ImportError('pplpy is not installed on your computer. Follow '
                                  'the instructions at https://gitlab.com/videlec/pplpy')

            n = self._triangulation.num_edges()
            M = ppl.MIP_Problem(2*n)
            x = [ppl.Variable(i) for i in range(n)]
            y = [ppl.Variable(n+i) for i in range(n)]

            self._set_train_track_constraints(M.add_constraint, x, VERTICAL, True, False)
            self._set_train_track_constraints(M.add_constraint, y, HORIZONTAL, True, False)
            self._set_geometric_constraints(M.add_constraint, x, y, True)

            return M.solve()['status'] != 'unfeasible'

        else:
            raise ValueError("method must be 'polytope' or 'LP'")

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
        edge_rep = self._edge_rep
        ep = self._triangulation.edge_permutation()
        
        if verbose:
            print('[edge_has_curve] checking edge %s with color %s' % (edge_rep(e), colouring[e]))

        a, b, c, d = self._triangulation.square_about_edge(e)
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

        n = self._triangulation.num_half_edges()
        fp = self._triangulation.face_permutation()
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

    def random_forward_flip(self, repeat=1):
        r"""
        Apply a forward flip randomly among the ones that keeps the triangulation core.

        INPUT:

        - ``repeat`` - integer (default 1) - if provided make ``repeat`` flips instead of 1.
        """
        from . import HAS_SAGE

        if HAS_SAGE:
            from sage.misc.prandom import choice, shuffle
        else:
            from random import choice, shuffle

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
                    self.back_flip(e, old_col)

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


