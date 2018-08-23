r"""
Triangulation of surfaces.
"""

from array import array

from .even_permutation import *

# These are broken forever - we have to look at ep to do this, so give
# up.  I haven't deleted them yet, but I will after I copy them to
# live inside of Triangulation (where they have access to ep).

#def edge_label(e):
#    if e < 0: return '~%d' % (~e)
#    else: return str(e)

#def norm(e):
#    return e if e >= 0 else ~e

class Triangulation(object):
    r"""
    A triangulation 

    with attributes:

    * _n  number of half-edges (an int)
    * _nf  number of folded edges (an int) # Perhaps we don't need this. 
    * _vp  vertex permutation (an array)
    * _ep  edge permutation (an array)
    * _fp  face permutation (an array)

    Each of vp, ep, fp is a permutation of half-edges, given as a
    function.  We allow folded edges; these are fixed by ep.  Thus the
    number of half-edges may be odd.

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
    __slots__ = ['_n', '_nf, ''_vp', '_ep', '_fp']

    def __init__(self, triangles, check=True):
        if isinstance(triangles, Triangulation):
            self._fp = triangles.face_permutation(copy=True)
            self._ep = triangles.edge_permutation(copy=True)
        else:
            self._fp, self._ep = perm_init(triangles)

        fp = self._fp
        ep = self._ep
        n = self._n = len(fp)
        nf = self._nf = len([i for i in range(n) if ep[i]==i])
        
        # TODO: face labels are disabled for now. We should actually
        # make it a *partition* of the faces (two faces are allowed to
        # have the same label) The trivial partition {0,1,2,...,F-1}
        # would correspond to no labels and the atomic
        # {{0},{1},...,{F-1}} would correspond to everybody being
        # labelled
        # fl = self._fl = [None] * F # face labels

        # TODO: edge labels are disabled for now. 
        # el = self._fl = [None] * E # edge labels

        v_base = [None] * n # Replaced 0 by None, as 0 is a possible edge label. 
        vp = self._vp = array('l', v_base)  

        for i in range(n):
            vp[fp[ep[i]]] = i

        # TODO: vertex labels are disabled for now
        # vl = self._vl = perm_cycles(vp)[1]....

        if check:
            self._check()

    def _check(self):
        n = self._n

        if not isinstance(self._fp, array) or len(self._fp) != n:
            raise RuntimeError('broken face perm')
        if not isinstance(self._ep, array) or len(self._ep) != n:
            raise RuntimeError('broken edge perm')
        if not isinstance(self._vp, array) or len(self._vp) != n:
            raise RuntimeError('broken vertex perm')
        perm_check(self._fp)
        perm_check(self._ep)
        perm_check(self._vp)

        # The face, edge, and vertex permutations fp, ep, and vp must
        # satisfy the relations fp^3 = ep^2 = fp.ep.vp = Id.  Note
        # that this is a representation of PSL(2, \ZZ) \isom S^2(2, 3,
        # \infty) into the symmetric group (on the half edges).

        for i in range(n):
            if self._fp[self._fp[self._fp[i]]] != i:
                raise RuntimeError('broken face permutation') # used to be ValueError
            if self._ep[self._ep[i]] != i:
                raise RuntimeError('broken edge permutation') 
            if self._fp[self._ep[self._vp[i]]] != i:
                raise RuntimeError('fev condition not satisfied')

    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n == other._n and self._fp == other._fp and self._ep == other._ep

    def __ne__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n == other._n and self._fp != other._fp and self._ep != other._ep
    
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

    def faces(self):
        r"""
        Return the list faces as triple of signed integers.

        EXAMPLES::

            sage: from veerer import *
            sage: t = Triangulation([1, 2, 0, -1, -3, -2])
            sage: t.faces()
            [(-3, -1, -2), (0, 1, 2)]
        """
        return list(self)

    def __iter__(self):
        n = self._n
        seen = [False] * n
        for e in range(n):
            if seen[e]:
                continue
            r = self._fp[e]
            s = self._fp[r]
            seen[e] = seen[r] = seen[s] = True
            yield (e, r, s)

    def __len__(self):
        return (self._n + self._nf) / 3

    def vertices(self):
        r"""
        Return the vertices ordered by their labels.
        """
        l = perm_cycles(self._vp)[0]
        return l

    def __repr__(self):
        l = self.faces()
        return '[' + ', '.join('(' + ', '.join(edge_label(e) for e in f) + ')' for f in l) + ']'

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
        T._vp = self._vp[:]
        T._fp = self._fp[:]
        return T

    def num_vertices(self):
        return len(self.vertices())

    def num_edges(self):
        return self._n

    def num_faces(self):
        return len(self.faces())

    def euler_characteristic(self):
        return self.num_faces() - self.num_edges()

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
        return 1 - (self.num_faces() - self.num_edges() + self.num_vertices()) / 2

    def is_flippable(self, e):
        e = int(e)
        return self._fp[e] != ~e and self._fp[~e] != e

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

        a = self._fp[e]
        b = self._fp[a]
        c = self._fp[~e]
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

        a = self._fp[e]
        b = self._fp[a]
        if a == ~e or b == ~e:
            raise ValueError('edge %s is not flippable' % edge_label(e))
        c = self._fp[~e]
        d = self._fp[c]

        # Disabled for now
        # F = self._fl[e]
        # G = self._fl[~e]

        # v = self._vl[b]
        # x = self._vl[d]

        # fix face perm and cycles
        self._fp[e] = b
        self._fp[b] = c
        self._fp[c] = e
        self._fp[a] = ~e
        self._fp[~e] = d
        self._fp[d] = a

        # Face labels
        # self._fl[a] = G
        # self._fl[c] = F

        # fix vertex perm
        self._vp[a] = ~d
        self._vp[b] = ~e
        self._vp[~e] = ~a
        self._vp[c] = ~b
        self._vp[d] = e
        self._vp[e] = ~c

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

        a = self._fp[e]
        b = self._fp[a]
        if a == ~e or b == ~e:
            raise ValueError('edge %s is not flippable' % edge_label(e))
        c = self._fp[~e]
        d = self._fp[c]

        # Disabled for now
        # F = self._fl[e]
        # G = self._fl[~e]

        # v = self._vl[b]
        # x = self._vl[d]

        # fix face perm and cycles
        self._fp[e] = d
        self._fp[d] = a
        self._fp[a] = e
        self._fp[b] = c
        self._fp[c] = ~e
        self._fp[~e] = b

        # Face labels
        # Disabled for now
        # self._fl[a] = G
        # self._fl[c] = F

        # fix vertex perm
        self._vp[a] = ~d
        self._vp[b] = e
        self._vp[e] = ~a
        self._vp[c] = ~b
        self._vp[d] = ~e
        self._vp[~e] = ~c

        # Vertex labels
        # Disabled for now
        # self._vl[e] = x
        # self._vl[~e] = v

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
