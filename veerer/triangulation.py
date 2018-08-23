r"""
Triangulation of surfaces.
"""

from array import array

from .even_permutation import *

# These are broken forever - we have to look at ep to do this, so give
# up.  I haven't deleted them yet, but I will after I copy them to
# live inside of Triangulation (where they have access to ep).

def edge_label(e):
    raise ValueError

def norm(e):
    raise ValueError
class Triangulation(object):
    r"""
    A triangulation of an oriented surface

    with attributes:

    * _n  number of half-edges (an int)
    * _fp  face permutation (an array)
    * _ep  edge permutation (an array)
    * _vp  vertex permutation (an array)

    Each of fp, ep, vp is a permutation of oriented half-edges (aka
    "darts"), given as a function.  Our conventions for the
    permutations are set out in the following figure:

          ~b           
        ----->         
    w-----------*-----------v
     \             <-----  /  
      \ \             b   / /
       \ \ c             / /~a
        \ \     F       / /
         \ v           / v 
          *         ^ *  
         ^ \       / /
          \ \   a / /
         ~c\ \   / /
            \ \   /
               \ /
                u

    Here the face permutation sends a to b, b to c, and c to a.  The
    vertex permutation send a to ~c, b to ~a, and c to ~b.  The edge
    permutation interchanges e and ~e for every edge; the edge e is
    folded if and only if e = ~e.  Thus folded edges are fixed by the
    edge permutation.

    EXAMPLES::

        sage: from veerer import *

        sage: T = Triangulation("(~2, 1, ~0)(~1, 0, 2)")
        sage: T.genus()
        1
        sage: T.num_faces()
        2
        sage: T.num_vertices()
        1
        sage: T.flip(0); T
        [(~2, 0, ~1), (~0, 1, 2)]
        sage: T.flip(0); T
        [(~2, 1, 0), (~1, ~0, 2)]
        sage: T.flip(0); T
        [(~2, ~0, ~1), (0, 1, 2)]
        sage: T.flip(0); T
        [(~2, 1, ~0), (~1, 0, 2)]

    """
    __slots__ = ['_n', '_fp', '_ep', '_vp']

    def __init__(self, triangles, check=True):
        if isinstance(triangles, Triangulation):
            self._fp = triangles.face_permutation(copy=True)
            self._ep = triangles.edge_permutation(copy=True)
        else:
            self._fp, self._ep = face_perm_init(triangles)

        fp = self._fp
        ep = self._ep
        n = self._n = len(fp)
        
        # TODO: face labels are disabled for now. We should actually
        # make it a *partition* of the faces (two faces are allowed to
        # have the same label) The trivial partition {0,1,2,...,F-1}
        # would correspond to no labels and the atomic
        # {{0},{1},...,{F-1}} would correspond to everybody being
        # labelled
        # fl = self._fl = [None] * F # face labels

        # TODO: edge labels are disabled for now. 
        # el = self._fl = [None] * E # edge labels

        v_base = [-1] * n 
        vp = self._vp = array('l', v_base)  

        for i in range(n):
            vp[fp[ep[i]]] = i

        # TODO: vertex labels are disabled for now
        # vl = self._vl = perm_cycles(vp)[1]....

        if check:
            self._check()

    def _check(self):
        n = self._n

        if not perm_check(self._fp, n):
            raise RuntimeError('fp is not a permtation')
        if not perm_check(self._ep, n):
            raise RuntimeError('ep is not permutation')
        if not perm_check(self._vp, n):
            raise RuntimeError('vp is not a permutation')

        # The face, edge, and vertex permutations fp, ep, and vp must
        # satisfy the relations fp^3 = ep^2 = fp.ep.vp = Id.  Note
        # that this is a representation of PSL(2, \ZZ) \isom S^2(2, 3,
        # \infty) into the symmetric group on the half edges.

        for i in range(n):
            if self._fp[self._fp[self._fp[i]]] != i:
                raise RuntimeError('broken face permutation')
            if self._ep[self._ep[i]] != i:
                raise RuntimeError('broken edge permutation') 
             if self._fp[self._ep[self._vp[i]]] != i:
                raise RuntimeError('fev relation not satisfied')

    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n == other._n and self._fp == other._fp and self._ep == other._ep

    def __ne__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n == other._n and self._fp != other._fp and self._ep != other._ep

    def copy(self):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation([[0,1,2],[-1,-2,-3]])
            sage: U = T.copy()
            sage: T == U
            True
            sage: T.flip(0)
            sage: T == U
            False
        """
        T = Triangulation.__new__(Triangulation)
        T._n = self._n
        T._fp = self._fp[:]
        T._ep = self._ep[:]
        T._vp = self._vp[:]
        return T

    def face_permutation(self, copy=True):
        if copy:
            return self._fp[:]
        else:
            return self._fp

    def edge_permutation(self, copy=True):
        if copy:
            return self._ep[:]
        else:
            return self._ep

    def vertex_permutation(self, copy=True):
        if copy:
            return self._vp[:]
        else:
            return self._vp

    def num_half_edges(self):
        return self._n

    def folded_edges(self):
        n = self._n
        ep = self._ep
        return [i for i in range(n) if ep[i]==i]
        
    def num_folded_edges(self):
        return len(self.folded_edges())

    def num_edges(self):
        return (self._n + self.num_folded_edges()) / 2

    def _edge_rep(self, e):
        f = self._ep[e]
        if f < e: return '~%d' % f
        else: return str(e)

    def _norm(self, e):
        f = self._ep[e]
        return f if f < e else e
    
    def num_faces(self):
        return self._n / 3

    def faces(self):
        r"""
        Return the list of faces as triples of half-edges

        EXAMPLES::

            sage: from veerer import *
            sage: t = Triangulation([1, 2, 0, -1, -3, -2]) # Fix
            sage: t.faces()
            [(-3, -1, -2), (0, 1, 2)]

        """
        return perm_cycles(self._fp)[0]

    def edges(self):
        r"""
        Return the list of faces as tuples of half-edges
        """
        l = perm_cycles(self._ep)[0]
        return l

    def vertices(self):
        r"""
        Return the list of faces as tuples of half-edges
        """
        l = perm_cycles(self._vp)[0]
        return l

    def num_vertices(self):
        return len(self.vertices())

    def euler_characteristic(self):
        return self.num_faces() - self.num_edges() + (self.num_vertices() + self.num_folded_edges())

    def genus(self):
        r"""
        Return the genus of this triangulation.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,~0)(2,~1,~2)")
            sage: T.genus()
            0

            sage: T = Triangulation([(-3,1,-1), (-2,0,2)])
            sage: T.genus()
            1
        """
        # 2 - 2g = \chi so 
        return (2 - self.euler_characteristic()) / 2

    def __repr__(self):
        faces = self.faces()        
        return '[' + ', '.join('(' + ', '.join(edge_label(e) for e in f) + ')' for f in l) + ']'

    def is_flippable(self, e):
        e = int(e)
        E = self._ep[e]
        a = self._fp[e]
        b = self._fp[a]
        return a != E and b != E

    def square_about_edge(self, e):
        # x<----------x
        # |     a    ^^
        # |         / |
        # |        /  |
        # |       /   |
        # |b    e/   d|
        # |     /     |
        # |    /      |
        # |   /       |
        # |  /        |
        # | /         |
        # v/    c     |
        # x---------->x

        e = int(e)
        E = self._ep[e]

        a = self._fp[e]
        b = self._fp[a]
        c = self._fp[E]
        d = self._fp[c]

        return a,b,c,d

    def flip(self, e):
        r"""
        Flip the edge ``e``.
        """
        # v<----------u     v<----------u
        # |     a    ^^     |^    a     ^
        # |         / |     | \         |
        # |  F     /  |     |  \     G  |
        # |       /   |     |   \       |
        # |b    e/   d| --> |b   \e    d|
        # |     /     |     |     \     |
        # |    /      |     |      \    |
        # |   /       |     |       \   |
        # |  /     G  |     | F      \  |
        # | /         |     |         \ |
        # v/    c     |     v     c    \|
        # w---------->x     w---------->x

        e = int(e)
        E = self._ep[e]

        a = self._fp[e]
        b = self._fp[a]
        if a == E or b == E:
            raise ValueError('edge %s is not flippable' % edge_label(e))
        c = self._fp[E]
        d = self._fp[c]

        A = self._ep[a]
        B = self._ep[b]
        C = self._ep[c]
        D = self._ep[d]
        
        # Disabled for now
        # F = self._fl[e]
        # G = self._fl[E]

        # v = self._vl[b]
        # x = self._vl[d]

        # fix face perm and cycles
        self._fp[e] = b
        self._fp[b] = c
        self._fp[c] = e
        self._fp[a] = E
        self._fp[E] = d
        self._fp[d] = a

        # Face labels
        # self._fl[a] = G
        # self._fl[c] = F

        # fix vertex perm
        self._vp[a] = D
        self._vp[b] = E
        self._vp[E] = A
        self._vp[c] = B
        self._vp[d] = e
        self._vp[e] = C

        # Vertex labels
        # self._vl[e] = x
        # self._vl[~e] = v

    def back_flip(self, e):
        r"""
        Flip back the edge ``e``.

        EXAMPLES::

            sage: from veerer import *
            sage: T0 = Triangulation([(0,1,2),(-1,-2,-3)])
            sage: T = T0.copy()
            sage: T.flip(0)
            sage: T.back_flip(0)
            sage: T == T0
            True

            sage: T.flip(1); T.flip(2)
            sage: T.back_flip(2); T.back_flip(1)
            sage: T == T0
            True
        """
        # Use the following for reference:
        # v<----------u     v<----------u
        # |     a    ^^     |\    a     ^
        # |         / |     | \         |
        # |  F     /  |     |  \     F  |
        # |       /   |     |   \       |
        # |b    e/   d| --> |b   \e    d|
        # |     /     |     |     \     |
        # |    /      |     |      \    |
        # |   /       |     |       \   |
        # |  /     G  |     | G      \  |
        # | /         |     |         \ |
        # v/    c     |     v     c    v|
        # w---------->x     w---------->x

        e = int(e)
        E = self._ep[e]

        a = self._fp[e]
        b = self._fp[a]
        if a == E or b == E:
            raise ValueError('edge %s is not flippable' % edge_label(e))
        c = self._fp[E]
        d = self._fp[c]

        A = self._ep[a]
        B = self._ep[b]
        C = self._ep[c]
        D = self._ep[d]
        
        # Disabled for now
        # F = self._fl[e]
        # G = self._fl[E]

        # v = self._vl[b]
        # x = self._vl[d]

        # fix face perm and cycles
        self._fp[e] = d
        self._fp[d] = a
        self._fp[a] = e
        self._fp[b] = c
        self._fp[c] = E
        self._fp[E] = b

        # Face labels
        # Disabled for now
        # self._fl[a] = G
        # self._fl[c] = F

        # fix vertex perm
        self._vp[a] = D
        self._vp[b] = e
        self._vp[e] = A
        self._vp[c] = B
        self._vp[d] = E
        self._vp[E] = C

        # Vertex labels
        # Disabled for now
        # self._vl[e] = x
        # self._vl[E] = v

    def to_string(self):
        r"""
        Serialize this triangulation as a string.

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation([(0,1,2),(-1,-2,-3)])
            sage: T.to_string()
            'efdcab'
        """
        return even_perm_base64_str(self._fp)

    @staticmethod
    def from_string(s, n):
        r"""
        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: Triangulation.from_string(T.to_string(), T.num_edges()) == T
            True
        """
        fp = even_perm_from_base64_str(s, n)
        return Triangulation(fp)

    def conjugate(self):
        r"""
        Conjugate this triangulation.

        The face perm is replaced by its inverse and conjugated by the
        ~ operation and the vertex permutation is replaced by its
        inverse.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)")
            sage: T.conjugate()
            sage: T
            [(~5, ~4, ~3), (~2, ~1, ~0), (0, 2, 4), (1, 3, 5)]
        """
        self._fp = even_perm_conjugate(even_perm_invert(self._fp))
        self._vp = even_perm_invert(self._vp)
