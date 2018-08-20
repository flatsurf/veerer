r"""
Triangulation of surfaces.
"""

from array import array

from .even_permutation import *

def edge_label(e):
    if e < 0: return '~%d' % (~e)
    else: return str(e)

def norm(e):
    return e if e >= 0 else ~e

class Triangulation(object):
    r"""
    A triangulation

    attributes

    * _n  number of edges (an int)
    * _vp  vertex permutations (an array)
    * _fp  face permutation (an array)

    EXAMPLES::

        sage: from veerer import *

        sage: preparser(False)

        sage: T = Triangulation([(-3, 1, -1), (-2, 0, 2)])
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
    __slots__ = ['_n', '_vp', '_fp']

    def __init__(self, triangles, check=True):
        if isinstance(triangles, Triangulation):
            triangles = triangles.faces()

        n = self._n = (3 * len(triangles)) / 2  # number of edges

        base = [0] * (2 * n)
        fp = self._fp = array('l', base)  # face permutation
        # fl = self._fl = [None] * (2 * n)  # face labels
        vp = self._vp = array('l', base)  # vertex permutation

        for i,t in enumerate(triangles):
            for j in range(3):
                # self._fl[t[j]] = i
                fp[t[j]] = t[(j+1) % 3]

        for i in range(-n, n):
            vp[fp[i]] = ~i
        # self._vl = even_perm_cycles(vp)[1]

        if check:
            self._check()

    def _check(self):
        n = self._n
        
        if not isinstance(self._fp, array) or len(self._fp) != 2*n:
            raise RuntimeError('wrong face perm')
        if not isinstance(self._vp, array) or len(self._vp) != 2*n:
            raise RuntimeError('wrong vertex perm')
        even_perm_check(self._fp)
        even_perm_check(self._vp)
        for i in range(-n, n):
            if ~self._vp[self._fp[i]] != i:
                raise RuntimeError('vef condition not satisfied')
            if self._fp[~self._vp[i]] != i:
                raise RuntimeError('vef condition not satisfied')

            # if self._vl[i] != self._vl[self._vp[i]]:
            #    raise RuntimeError("vl bad")
            # if self._fl[i] != self._fl[self._fp[i]]:
            #    raise RuntimeError("fl bad")

    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n == other._n and self._fp == other._fp

    def __ne__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n != other._n or self._fp != other._fp
    
    def vertex_permutation(self):
        return self._vp

    def face_permutation(self):
        return self._fp

    def faces(self):
        r"""
        Return the faces ordered by their labels.
        """
        return list(self)

    def __iter__(self):
        n = self._n
        seen = [False] * (2 * n)
        for e in range(-n, n):
            if seen[e]:
                continue
            r = self._fp[e]
            s = self._fp[r]
            seen[e] = seen[r] = seen[s] = True
            yield (e, r, s)

    def __len__(self):
        return (2 * self._n) / 3

    def vertices(self):
        r"""
        Return the vertices ordered by their labels.
        """
        l = even_perm_cycles(self._vp)[0]
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

            sage: T = Triangulation([(0,1,-1),(2,-2,-3)])
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
        # Use the following for reference:
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

        # TODO: below is wrong
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
        # TODO: below is wrong
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
        # TODO: below is wrong
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

            sage: T = Triangulation([(-12, 4, -4), (-11, -1, 11), (-10, 0, 10),
            ....:                    (-9, 9, 1), (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)])
            sage: Triangulation.from_string(T.to_string(), T.num_edges()) == T
            True
        """
        fp = even_perm_from_base64_str(s, n)
        faces = even_perm_cycles(fp)[0]
        return Triangulation(faces)

