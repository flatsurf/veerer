r"""
Covering of triangulations
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

from array import array

from .permutation import *
from .triangulation import Triangulation


class TriangulationCover(object):
    r"""
    A triangulation cover

    We consider ramified cover with arbitrary ramifications at the vertices and
    at most double ramification at the middle of folded edges. A covering of
    degree d of a triangulation with n half edges is determined by a tuple of n
    permutations in S_d.

    As this is also a triangulation, we choose a canonical labeling of its
    edges as follows: i + n*j for the i-th edge in the j-th copy of the
    original triangulation. This is used in homology computations.

    This object has attributes:

    * _t: underlying triangulation of type ``Triangulation``
    * _d: degree of the cover of type ``int``
    * _c: covering cocycle, a list of permutations, one per half-edge

    EXAMPLES:

    The torus as a double of a triangle::

        sage: from veerer import Triangulation  # random output due to deprecation warnings from realalg
        sage: T = Triangulation("(0,1,2)")
        sage: C = T.cover([[1,0], [1,0], [1,0]])
        sage: C.euler_characteristic()
        0
        sage: C.as_triangulation().euler_characteristic()
        0

    The quadratic L-shape surface::

        sage: from veerer import Triangulation
        sage: T = Triangulation("(0,1,2)")
        sage: C = T.cover([[1,0,2], [2,1,0], [0,1,2]])
        sage: C.euler_characteristic()
        2
        sage: C.as_triangulation().euler_characteristic()
        2
    """
    __slots__ = ['_t', '_d', '_c']

    def __init__(self, triangulation, cover, mutable=False, check=True):
        if not isinstance(triangulation, Triangulation):
            self._t = Triangulation(triangulation, mutable=mutable, check=False)
        else:
            self._t = triangulation.copy(mutable=mutable)

        if not isinstance(cover, (tuple, list)) or \
           (len(cover) != self._t.num_half_edges() and len(cover) != self._t.num_edges()):
            raise ValueError("the argument 'cover' must be a list with as many elements as the number of edges (or half edges) of the triangulation")

        d = None
        cover = list(cover)
        for i,p in enumerate(cover):
            if p is not None:
                cover[i] = p = perm_init(p)
                if d is None:
                    d = len(p)
                elif len(p) != d:
                    raise ValueError("inconsistent cover degree")
        if d is None:
            raise ValueError("can not determine cover degree")
        for i,p in enumerate(cover):
            if p is None:
                cover[i] = perm_init(range(d))

        ep = self._t.edge_permutation(copy=False)
        n = self._t.num_half_edges()
        ne = self._t.num_edges()
        nf = self._t.num_folded_edges()
        if len(cover) == ne:
            for i in range(ne, n):
                if ep[i] >= ne:
                    raise ValueError("non standard edge labeling, you need to provide the full list of permutations for the cover")
                cover.append(perm_invert(cover[ep[i]]))

        self._d = d
        self._c = cover

        if check:
            self._check(ValueError)

    def _check(self, error=RuntimeError):
        n = self._t.num_half_edges()
        d = self._d

        if len(self._c) != n:
            raise error("wrong length for _c attribute")

        for p in self._c:
            if not perm_check(p, d):
                raise error("invalid %d-th permutation %s in covering data" % (d, p))

        ep = self._t.edge_permutation(copy=False)
        for i,p in enumerate(self._c):
            q = self._c[ep[i]]
            if perm_invert(p) != q:
                raise error("inconsistent covering data on edge (%d,%d)" %(i, ep[i]))

        # NOTE: should we check connectedness?

    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._t == other._t and self._c == other._c

    def __ne__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._t != other._t or self._c != other._c

    def copy(self, mutable=None):
        T = TriangulationCover.__new__(TriangulationCover)
        T._t = self._t.copy(mutable=mutable)
        T._d = self._d
        T._c = tuple(p[:] for p in self._c)
        return T

    def __str__(self):
        t = self._t
        ep = t.edge_permutation(copy=False)
        n = t.num_half_edges()
        faces = t.faces()
        fstr = ''.join('(' + ','.join(t._edge_rep(e) for e in f) + ')' for f in faces)
        cstr = ',\n  '.join('[' + ','.join(str(i) for i in self._c[e]) + ']' for e in range(n) if ep[e] >= e)
        return 'TriangulationCover("%s",\n [%s])' % (fstr, cstr)

    def __repr__(self):
        return str(self)

    def base(self, copy=True):
        if copy:
            return self._t.copy()
        else:
            return self._t

    ############
    # topology #
    ############

    def degree(self):
        return self._d

    def vertex_permutation(self):
        n = self._t.num_half_edges()
        d = self._d
        vp = self._t.vertex_permutation(copy=False)
        ep = self._t.edge_permutation(copy=False)
        cvp = array('i', [-1] * (d * n))

        for a in range(n):
            b = vp[a]
            p = self._c[ep[b]]
            for i in range(d):
                ai = a + n * i
                bj = b + n * p[i]
                cvp[ai] = bj

        return cvp

    def vertices(self):
        return perm_cycles(self.vertex_permutation())

    def edge_permutation(self):
        n = self._t.num_half_edges()
        d = self._d
        ep = self._t.edge_permutation(copy=False)
        c = self._c
        cep = array('i', [-1] * (d * n))

        for a in range(n):
            b = ep[a]
            p = c[a]
            for j in range(d):
                cep[a + n * j] = b + n * p[j]

        return cep

    def edges(self):
        return perm_cycles(self.edge_permutation())

    def face_permutation(self):
        n = self._t.num_half_edges()
        d = self._d
        fp = self._t.face_permutation(copy=False)
        cfp = array('i', [-1] * (d * n))

        for a in range(n):
            b = fp[a]
            for i in range(d):
                cfp[a + n * i] = b + n * i

        return cfp

    def faces(self):
        return perm_cycles(self.face_permutation())

    def as_triangulation(self, mutable=None):
        r"""
        Return this cover as a triangulation.

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,~1)(~0,2,3)(~2,4,~4)")
            sage: C = T.cover([[0,1,3,2],[1,0,3,2],[1,2,3,0],[3,2,1,0],[0,1,3,2]])
            sage: C.as_triangulation()
            Triangulation("(0,1,7)(2,3,~0)(4,~4,6)...(~6,~3,~18)(~23,~22,~20)")
        """
        if mutable is None:
            mutable = self._t._mutable
        return Triangulation.from_face_edge_perms(
                fp=self.face_permutation(),
                ep=self.edge_permutation(),
                vp=self.vertex_permutation(),
                mutable=mutable)

    def num_folded_edges(self):
        n = self._t.num_half_edges()
        d = self._d
        ep = self._t.edge_permutation(copy=False)
        c = self._c
        nfe = 0
        for i in range(n):
            if ep[i] != i:
                continue
            for j in range(d):
                nfe += c[i][j] == j
        return nfe

    def num_edges(self):
        return ((self._d * self._t.num_half_edges()) + self.num_folded_edges()) // 2

    def num_faces(self):
        return self._t.num_faces() * self._d

    def num_vertices(self):
        return len(self.vertices())

    def euler_characteristic(self):
        return self.num_faces() - self.num_edges() + (self.num_vertices() + self.num_folded_edges())

    def genus(self):
        # 2 - 2g = \chi so
        return (2 - self.euler_characteristic()) // 2

    ########
    # flip #
    ########

    def flip(self, e):
        r"""
        Flip the edge ``e`` in the underlying triangulation.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)")
            sage: C = T.cover([[1,0], [1,0], [1,0]], mutable=True)
            sage: C.flip(0)
            sage: C
            TriangulationCover("(0,2,1)",
             [[1,0],
              [1,0],
              [1,0]])

            sage: T = Triangulation("(0,1,2)")
            sage: C = T.cover([[1,0,2], [2,1,0], [0,1,2]], mutable=True)
            sage: C.flip(0); C._check()
            sage: C.flip(1); C._check()
            sage: C.flip(2); C._check()

        TESTS::

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: n = T.num_half_edges()
            sage: d = 5
            sage: for e in range(n):
            ....:     if T.is_flippable(e):
            ....:         C = T.cover([[2,4,1,3,0],[4,1,3,2,0],[3,0,1,2,4]], mutable=True)
            ....:         TT = C.as_triangulation()
            ....:         C.flip(e)
            ....:         C._check()
            ....:         for i in range(d):
            ....:             assert TT.is_flippable(e + n * i)
            ....:             TT.flip(e + n * i)
            ....:         assert C.as_triangulation().is_isomorphic_to(TT)
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
        if not self._t._mutable:
            raise ValueError("immutable veering triangulation cover; use a mutable copy instead")

        e = int(e)
        ep = self._t.edge_permutation(copy=False)
        fp = self._t.face_permutation(copy=False)
        E = ep[e]
        a = fp[e]
        A = ep[a]
        c = fp[E]
        C = ep[c]

        self._t.flip(e)

        # TODO: do the compositions inplace
        self._c[a] = perm_compose(self._c[E], self._c[a])
        self._c[A] = perm_compose(self._c[A], self._c[e])
        self._c[c] = perm_compose(self._c[e], self._c[c])
        self._c[C] = perm_compose(self._c[C], self._c[E])

    ################################
    # relabeling and automorphisms #
    ################################

    # In a triangle (a, b, c) we are allowed to multiply the three
    # permutations (theta_a, theta_b, theta_c) that define the cover
    # by a common permutation g, (g theta_a, g theta_b, g theta_c)

    # To make a canonical numbering, we choose a spanning tree of the dual and
    # set the permutations on this spanning tree to be identity


    ###################
    # homology action #
    ###################

    def _check_homology_matrix(self, m):
        pass
