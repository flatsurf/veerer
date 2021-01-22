r"""
Veering triangulations endowed with a flat structure.
"""

from sage.structure.sequence import Sequence
from sage.categories.fields import Fields
from sage.rings.all import ZZ, QQ, AA, RDF, NumberField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector

from .constants import BLUE, RED, PURPLE, GREEN
from .permutation import perm_cycles
from .triangulation import Triangulation
from .veering_triangulation import VeeringTriangulation
from .misc import flipper_edge, flipper_edge_perm, flipper_nf_to_sage, flipper_nf_element_to_sage, det2, flipper_face_edge_perms

_Fields = Fields()

def vec_slope(v):
    r"""
    Return the slope of a 2d vector ``v``.

    EXAMPLES::

        sage: from veerer.constants import RED, BLUE, PURPLE, GREEN
        sage: from veerer.flat_structure import vec_slope

        sage: vec_slope((1,0)) == PURPLE
        True
        sage: vec_slope((0,1)) == GREEN
        True
        sage: vec_slope((1,1)) == vec_slope((-1,-1)) == RED
        True
        sage: vec_slope((1,-1)) == vec_slope((1,-1)) == BLUE
        True
    """
    if v[0].is_zero():
        return GREEN
    elif v[1].is_zero():
        return PURPLE
    elif v[0] * v[1] > 0:
        return RED
    else:
        return BLUE

class FlatVeeringTriangulation(Triangulation):
    r"""
    A triangulation with flat structures with veering compatible slopes.

    The vectors are kept coherently within triangles (ie a+b+c = 0). A pair of
    edges (e, E) can either be glued via translation or point symmetry. If the
    surface is a translation surface, these are only translations.

    EXAMPLES::

        sage: from veerer import FlatVeeringTriangulation

        sage: vecs = [(1, 2), (-2, -1), (1, -1), (1, -1), (-2, -1), (1, 2)]
        sage: FlatVeeringTriangulation("(0,1,2)(~0,~1,~2)", vecs)
        FlatVeeringTriangulation(Triangulation("(0,1,2)(~2,~0,~1)"), [(1, 2), (-2, -1), (1, -1), (-1, 1), (2, 1), (-1, -2)])

    For translation structures (or Abelian differential) the holonomies might be modified
    by a sign so that ``holonomies[e] = - holonomies[~e]``::

        sage: vecs = [(1, 2), (-2, -1), (1, -1), (-1, 1), (2, 1), (-1, -2)]
        sage: vecs = [vector(ZZ, v) for v in vecs]
        sage: fp = "(0,1,2)(~0,~1,~2)"
        sage: FlatVeeringTriangulation(fp, vecs)
        FlatVeeringTriangulation(Triangulation("(0,1,2)(~2,~0,~1)"), [(1, 2), (-2, -1), (1, -1), (-1, 1), (2, 1), (-1, -2)])
        sage: vecs = [vecs[0], vecs[1], vecs[2], -vecs[3], -vecs[4], -vecs[5]]
        sage: FlatVeeringTriangulation(fp, vecs)
        FlatVeeringTriangulation(Triangulation("(0,1,2)(~2,~0,~1)"), [(1, 2), (-2, -1), (1, -1), (-1, 1), (2, 1), (-1, -2)])
        sage: vecs = [-vecs[0], -vecs[1], -vecs[2], vecs[3], vecs[4], vecs[5]]
        sage: FlatVeeringTriangulation(fp, vecs)
        FlatVeeringTriangulation(Triangulation("(0,1,2)(~2,~0,~1)"), [(1, 2), (-2, -1), (1, -1), (-1, 1), (2, 1), (-1, -2)])
        sage: vecs = [-vecs[0], -vecs[1], -vecs[2], vecs[3], vecs[4], vecs[5]]
        sage: FlatVeeringTriangulation(fp, vecs)
        FlatVeeringTriangulation(Triangulation("(0,1,2)(~2,~0,~1)"), [(1, 2), (-2, -1), (1, -1), (-1, 1), (2, 1), (-1, -2)])
    """
    def __init__(self, triangulation, holonomies=None, base_ring=None, check=True):
        Triangulation.__init__(self, triangulation, check=False)
        if isinstance(triangulation, FlatVeeringTriangulation):
            self._holonomies = triangulation._holonomies[:]
            self._base_ring = triangulation._K
            self._V = triangulation._V
            self._K = triangulation._K
            self._translation = triangulation._translation
            return 

        if base_ring is None:
            S = Sequence([vector(v) for v in holonomies])
            self._V = S.universe()
            self._K = self._V.base_ring()
            holonomies = list(S)
        else:
            self._K = base_ring
            self._V = VectorSpace(self._K, 2)
            holonomies = [self._V(v) for v in holonomies]

        if self._K not in _Fields:
            self._K = self._K.fraction_field()
            self._V = self._V.change_ring(self._K)
            holonomies = [v.change_ring(self._K) for v in holonomies]

        n = self._n
        m = self.num_edges()
        ep = self._ep

        if len(holonomies) == m:
            for e in range(m):
                E = ep[e]
                if e != E and E < m:
                    raise ValueError("edge perm not in standard form")
            holonomies.extend([-holonomies[ep[e]] for e in range(m,n)])
        if len(holonomies) != n:
            raise ValueError('wrong number of vectors')
        self._holonomies = holonomies

        cols = [vec_slope(self._holonomies[e]) for e in range(n)]
        if isinstance(triangulation, VeeringTriangulation):
            # check that colours are compatible
            for e in range(triangulation.num_edges()):
                tcol = triangulation.edge_colour(e)
                scol = self.edge_colour(e)
                if scol == PURPLE or scol == GREEN:
                    continue
                if tcol != scol:
                    raise ValueError("incompatible colours")

        V = self.to_veering_triangulation()
        ans, cert = V.is_abelian(certificate=True)
        self._translation = False
        if ans:
            # translation surface (Abelian differential)
            # fix holonomies so that
            # holonomies[ep[e]] = - holonomies[e]
            self._translation = True
            ep = self._ep
            for e, right in enumerate(cert):
                E = ep[e]
                if right != (self._holonomies[e][0] > 0):
                    self._holonomies[e] = -self._holonomies[e]
                if right != (self._holonomies[E][0] < 0):
                    self._holonomies[E] = -self._holonomies[E]

        if check:
            self._check()

    @classmethod
    def from_coloured_triangulation(cls, T):
        r"""
        Construct a flat triangulation associated to a given coloured triangulation.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
            sage: FlatVeeringTriangulation.from_coloured_triangulation(T)
            FlatVeeringTriangulation(Triangulation("(0,1,2)(~2,~0,~1)"), [(1, 2), (-2, -1), (1, -1), (-1, 1), (2, 1), (-1, -2)])
        """
        return T.flat_structure_min()

    def triangle_upside_down(self, e):
        r"""
        Apply a 180 degree rotation to the triangle containing the edge ``e``.
        """
        f = self._fp[e]
        g = self._fp[f]
        self._holonomies[e] = -self._holonomies[e]
        self._holonomies[f] = -self._holonomies[f]
        self._holonomies[g] = -self._holonomies[g]

        self._check()

    def colours_about_edge(self, e):
        e = int(e)
        return [self.edge_colour(f) for f in self.square_about_edge(e)]

    def alternating_square(self, e):
        r"""
        Return whether there is an alternating square around the edge ``e``.
        """
        e = int(e)
        colours = self.colours_about_edge(e)
        if any(colours[f] == GREEN or colours[f] == PURPLE for f in range(4)):
            return False
        return all(colours[f] != colours[(f+1) % 4] for f in range(4))

    def is_flippable(self):
        return Triangulation.is_flippable(self, e) and self.alternating_square()

    def is_forward_flippable(self, e):
        return Triangulation.is_flippable(self, e) and self.colours_about_edge(e) == [BLUE, RED, BLUE, RED]

    def forward_flippable_edges(self):
        ep = self._ep
        n = self._n
        return [e for e in range(n) if e <= ep[e] and self.is_forward_flippable(e)]

    def is_backward_flippable(self, e):
        return Triangulation.is_flippable(self, e) and self.colours_about_edge(e) == [RED, BLUE, RED, BLUE]

    def backward_flippable_edges(self):
        ep = self._ep
        n = self._n
        return [e for e in range(n) if e <= ep[e] and self.is_backward_flippable(e)]

    def swap(self, e):
        r"""
        Swap the orientation of the edge ``e``.
        """
        E = self._ep[e]
        if e != E:
            Triangulation.swap(e)
            self._holonomies[e], self._holonomies[E] = self._holonomies[E], self._holonomies[e]

    def boshernitzan_criterion(self):
        r"""
        """
        if not self._translation:
            raise NotImplementedError

        n = self.num_edges()
        X = []
        for h in self._holonomies[:n]:
            if h[1] == 0:
                continue
            elif h[1] < 0:
                X.append(-h[0])
            else:
                X.append(h[0])
        X = [x.vector() for x in X]
        m = len(X[0])
        X = [[0] + [X[i][j] for i in range(len(X))] for j in range(m)]

        # non-negativity
        ieqs = []
        for i in range(n):
            e = [0] * n
            e[i] = 1
            ieqs.append([0] + e)

        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(eqns=X, ieqs=ieqs)

    @classmethod
    def from_flipper_pseudo_anosov(cls, h):
        r"""
        Construct the flat structure from a pseudo-Anosov homeomorphism.

        EXAMPLES::

            sage: from flipper import *  # optional - flipper
            sage: from veerer import *

            sage: T = flipper.load('SB_4')  # optional - flipper
            sage: h = T.mapping_class('s_0S_1s_2S_3s_1S_2').canonical()     # optional - flipper
            sage: F = FlatVeeringTriangulation.from_flipper_pseudo_anosov(h)  # optional - flipper
            sage: F                                                         # optional - flipper
            FlatVeeringTriangulation(Triangulation("(0,4,~1)(1,5,3)(2,~0,~3)(~5,~4,~2)"), [(-2, -1), (-2*a + 4, a - 1), (-2*a + 4, a - 1), (2*a - 2, -a + 2), (2*a - 2, -a + 2), (-2, -1), (-2, -1), (2*a - 2, -a + 2), (2*a - 2, -a + 2), (-2*a + 4, a - 1), (-2*a + 4, a - 1), (-2, -1)])

        The flipper and veerer triangulations carry the same edge labels::

            sage: F.to_veering_triangulation()         # optional - flipper
            VeeringTriangulation("(0,4,~1)(1,5,3)(2,~0,~3)(~5,~4,~2)", "RBBBBR")
            sage: h.source_triangulation   # optional - flipper
            [(~5, ~4, ~2), (~3, 2, ~0), (~1, 0, 4), (1, 5, 3)]
        """
        Fh = h.flat_structure()
        Th = Fh.triangulation
        n = 3 * Th.num_triangles # number of half edges

        fp, ep = flipper_face_edge_perms(Th)
        T = Triangulation.from_face_edge_perms(fp, ep)

        # extract flat structure
        x = next(iter(Fh.edge_vectors.values())).x
        K = flipper_nf_to_sage(x.field)
        V = VectorSpace(K, 2)
        # translate into Sage number field
        Vec = {flipper_edge(Th,e): V((flipper_nf_element_to_sage(v.x, K),
                               flipper_nf_element_to_sage(v.y, K)))
                  for e,v in Fh.edge_vectors.items()}
        vectors = [Vec[i] for i in range(n)]

        return FlatVeeringTriangulation(T, vectors, K)

    def to_veering_triangulation(self):
        return VeeringTriangulation(self, [self.edge_colour(e) for e in range(self._n)])

    def to_pyflatsurf(self):
        if not self._translation:
            raise ValueError("pyflatsurf only works with translation surfaces")
        from pyflatsurf.factory import make_surface
        from pyflatsurf.sage_conversion import make_vectors
        verts = [(i+1 if i >=0 else i for i in c) for c in perm_cycles(self._vp, True, self._n)]
        vectors = makeVectors(self._holonomies[e] for e in range(self.num_edges()))
        return make_surface(verts, vectors)

    def __str__(self):
        return "FlatVeeringTriangulation({}, {})".format(Triangulation.__str__(self),
                  self._holonomies)

    def __repr__(self):
        return str(self)

    def _check(self):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,3)", "BRRR")
            sage: assert T.is_core()
            sage: F = T.flat_structure_min()
            sage: F._check()
        """
        Triangulation._check(self)

        n = self.num_half_edges()
        ep = self.edge_permutation(copy=False)
        vectors = self._holonomies
        if len(vectors) != n:
            raise ValueError("invalid list of vectors")

        for a in range(n):
            A = ep[a]
            u = vectors[a]
            v = vectors[A]
            if u != v and u != -v:
                raise ValueError('ep[%s] = %s but vec[%s] = %s' % (a, u, A, v))

        for a,b,c in self.faces():
            va = vectors[a]
            vb = vectors[b]
            vc = vectors[c]
            if va + vb + vc:
                raise ValueError('vec[%s] = %s, vec[%s] = %s and vec[%s] = %s do not sum to zero' % (a, va, b, vb, c, vc))

            if det2(va, vb) <= 0 or det2(vb, vc) <= 0 or det2(vc, va) <= 0:
                raise ValueError('(%s, %s, %s) is a clockwise triangle' %
                        (a, b, c))

        if self._translation:
            for e in range(self.num_edges()):
                E = ep[e]
                assert self._holonomies[e] == -self._holonomies[E]

    def copy(self):
        res = FlatVeeringTriangulation.__new__(FlatVeeringTriangulation)
        res._n = self._n
        res._vp = self._vp[:]
        res._ep = self._ep[:]
        res._fp = self._fp[:]
        res._V = self._V
        res._K = self._K
        res._holonomies = [v.__copy__() for v in self._holonomies]
        res._translation = self._translation
        return res

    def edge_colour(self, e):
        return vec_slope(self._holonomies[e])

    def layout(self):
        r"""
        Return a layout object that can be used to produce various plots of the surface.
        """
        from .layout import FlatVeeringTriangulationLayout
        return FlatVeeringTriangulationLayout(self)

    def plot(self, *args, **kwds):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,3)", "BRRR")
            sage: F = T.flat_structure_min()
            sage: F.plot()
            Graphics object consisting of 15 graphics primitives

            sage: F.plot(horizontal_train_track=True)
            Graphics object consisting of 19 graphics primitives
            sage: F.plot(vertical_train_track=True)
            Graphics object consisting of 19 graphics primitives
            sage: F.plot(horizontal_train_track=True, vertical_train_track=True)
            Graphics object consisting of 23 graphics primitives
        """
        return self.layout().plot(*args, **kwds)

    def flip(self, e):
        r"""
        Flip the edge ``e``.

        EXAMPLES::

            sage: from veerer import Triangulation, FlatVeeringTriangulation
            sage: T = Triangulation("(0,1,2)(3,4,~0)(5,6,~1)")
            sage: fl = FlatVeeringTriangulation(T, [(-47, 51), (-27, -67), (74, 16), (22, 79), (-69, -28), (-61, -31), (34, -36)])
            sage: fl.flip(2)
            sage: fl.flip(4)
            sage: fl.flip(0)
            sage: fl
            FlatVeeringTriangulation(Triangulation("(0,1,4)(2,~0,3)(5,6,~1)"), [(2, 197), (-27, -67), (-20, 118), (22, 79), (25, -130), (-61, -31), (34, -36), (27, 67), (-2, -197)])
        """
        if not self.is_forward_flippable(e):
            raise ValueError

        E = self._ep[e]
        if e == E:
            # folded edge
            a = self._fp[e]
            b = self._fp[a]
            self._holonomies[a] = -self._holonomies[a]
            self._holonomies[e] = -(self._holonomies[a] + self._holonomies[b])
            Triangulation.flip(self, e)

        else:
            if self._holonomies[e] == self._holonomies[E]:
                # Make it so that e is well oriented
                self.triangle_upside_down(e)
            assert self._holonomies[e] == -self._holonomies[E]

            a = self._fp[e]
            b = self._fp[a]
            c = self._fp[E]
            d = self._fp[c]

            assert self._holonomies[d] + self._holonomies[a] + self._holonomies[b] + self._holonomies[c] == 0

            Triangulation.flip(self, e)
            self._holonomies[e] = self._holonomies[d] + self._holonomies[a]
            self._holonomies[E] = -self._holonomies[e]

        self._check()
