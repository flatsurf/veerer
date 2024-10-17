# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2024 Vincent Delecroix
#                     2024 Kai Fu
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


import numbers
from array import array

from sage.rings.all import ZZ

from .permutation import perm_cycles, perm_cycles_to_string, str_to_cycles_and_data
from .triangulation import face_edge_perms_init, boundary_init
from .constellation import Constellation


class StrebelGraph(Constellation):
    r"""
    Strebel graph.

    A Strebel graph encodes a specific saddle connection configuration
    associated to a meromorphic Abelian or quadratic differential on a Riemann
    surface.

    EXAMPLES::

        sage: from veerer import StrebelGraph

        sage: StrebelGraph("(0,~0:1)")
        StrebelGraph("(0,~0:1)")
    """
    def __init__(self, faces, mutable=False, check=True):
        if isinstance(faces, StrebelGraph):
            fp = faces.face_permutation(copy=True)
            ep = faces.edge_permutation(copy=True)
            bdry = faces.boundary_vector(copy=True)
        else:
            faces, boundary = str_to_cycles_and_data(faces)
            fp, ep = face_edge_perms_init(faces)
            bdry = boundary_init(fp, ep, boundary)

        Constellation.__init__(self, len(fp), None, ep, fp, (bdry,), mutable, check)

    def _set_data_pointers(self):
        self._bdry = self._data[0]

    def __str__(self):
        r"""
        EXAMPLES::

            sage: from veerer import StrebelGraph
            sage: T = StrebelGraph("(0,1,2)(~0,~1:1,~2:2)")
            sage: str(T)
            'StrebelGraph("(0,1,2)(~2:2,~0,~1:1)")'
        """
        bdry_cycles = perm_cycles_to_string(perm_cycles(self._fp, n=self._n), involution=self._ep, data=self._bdry)
        return 'StrebelGraph("%s")' % bdry_cycles

    def __repr__(self):
        return str(self)

    def flip(self, e, check=True):
        r"""
        Flip a half-edge of the Strebel graph.

        All half-plane excesses are required to be zero within this method.

        EXAMPLES::

            sage: from veerer import StrebelGraph

            sage: G = StrebelGraph("(0,1,2,3)(~0,4,5,6)(~1,~6,~5,~4,~3,~2)", mutable=True)
            sage: G.flip(0)
            sage: G1 = StrebelGraph("(0,2,3,4)(~0,5,6,1)(~1,~6,~5,~4,~3,~2)")
            sage: G == G1
            True
        """
        # v<----------u     v<----------u
        # |     a    ^^     |^    a     ^
        # |b        / |     |b\         |
        # v  F     / k|     v  \     G k|
        # |       /   ^     |   \       ^
        # |     e/    | --> |    \      |
        # |     /     |     |     \     |
        # v    /      |     v     e\    |
        # |h  /       ^     |h      \   ^
        # |  /     G  |     | F      \  |
        # | /        d|     |         \d|
        # v/    c     |     v     c    \|
        # w---------->x     w---------->x

        n = self._n

        for i in range(n):
            if self._bdry[i] != 0:
                raise ValueError('graph with boundary faces; we do not consider boundary face in this flipping')

        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        if check:
            e = self._check_half_edge(e)

        E = self._ep[e]

        a = self._fp[e]
        b = self._fp[a]
        if a == E or b == E:
            raise ValueError('edge %s is not flippable' % self._norm(e))
        c = self._fp[E]
        d = self._fp[c]

        A = self._ep[a]
        B = self._ep[b]
        C = self._ep[c]
        D = self._ep[d]

        H = self._vp[e]
        h = self._ep[H]
        K = self._vp[E]
        k = self._ep[K]

        # fix face perm and cycles
        self._fp[e] = b
        self._fp[h] = c
        self._fp[c] = e
        self._fp[a] = E
        self._fp[E] = d
        self._fp[k] = a

        # fix vertex perm
        self._vp[a] = K
        self._vp[b] = E
        self._vp[E] = A
        self._vp[c] = H
        self._vp[d] = e
        self._vp[e] = C

    def halfdualflip(self, e, check = True):
        r"""
        flip the half-edge ``e`` of a strebel graph.

        EXAMPLES::

            sage: from veerer import StrebelGraph

        Flip a half-edge with zero half-plane excess::

            sage: G = StrebelGraph("(0,1,2)(~0)(~1:1,~2:1)", mutable=True)
            sage: G.halfdualflip(0)
            sage: G1 = StrebelGraph("(0,~2,~1:1)(1,2:1)(~0)")
            sage: G == G1
            True

        Flip a half-edge with non-zero half-plane excess::

            sage: G = StrebelGraph("(0,1,2)(~0:0)(~1:1,~2:1)", mutable=True)
            sage: G.halfdualflip(4)
            sage: G1 = StrebelGraph("(0,1,2)(~1:0,~2:2)(~0)")
            sage: G == G1
            True

        The graph remains the same if ``e`` is in a loop edge::

            sage: G = StrebelGraph("(0,1,2)(~0:0)(~1:1,~2:1)", mutable=True)
            sage: G1 = G.copy()
            sage: G.halfdualflip(5)
            sage: G == G1
            True

        The graph remains the same if ``e`` is an end of the graph with half-plane excess::

            sage: G = StrebelGraph("(0,~1,1,~2,2,~0)", mutable=True)
            sage: G1 = G.copy()
            sage: G.halfdualflip(0)
            sage: G == G1
            True
        """
        #
        #       |A                         \A
        #       |                           \
        #       v                            v
        #       |             -->             \
        #       |                              \
        #       ^           ^                   ^   ^
        #      a|          c|                  a \ c|
        #       |           |                     \ |
        #       |  e        |              e       \|
        # <-----w-----------x     <-----w-----------x
        #    b         E            b           E
        #  or
        #                ^                          ^
        #                |                          |
        #               c|   -->                   c|
        #  angle0  angle1|        angle0-1  angle1+1|
        #  w-------------x        w-----------------x


        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        if check:
            e = self._check_half_edge(e)

        a = self._vp[e]
        b = self._vp[a]
        c = self._fp[e]
        E = self._ep[e]
        A = self._ep[a]
        B = self._ep[b]

        angle0 = self._bdry[e]
        angle1 = self._bdry[c]

        if angle0 == 0:
            # we do not flip anything when a is the same as e or E(simple pole)

            if (a == e) or (a == E):
                pass

            else:
                self._fp[B] = e
                self._fp[e] = a
                self._fp[A] = c

                self._vp[e] = b
                self._vp[c] = a
                self._vp[a] = E
        else:
            self._bdry[e] = angle0 - 1
            self._bdry[c] = angle1 + 1

    def dualflip(self,e,check = True):
        r"""
        flip the edge of the Strebel graph containg the half-edge e

        EXAMPLES::

            sage: from veerer import *

            sage: G = StrebelGraph("(0,1,2)(~0:0)(~1:1,~2:1)", mutable=True)
            sage: G.dualflip(1)
            sage: G1 = StrebelGraph("(0,2)(~0,1)(~1,~2:2)")
            sage: G == G1
            True

        We obtain the same Strebel graphs no matter we input ``e``
        or the edge permutation of ``e``::

            sage: G = StrebelGraph("(0,1,2)(~0:0)(~1:1,~2:1)", mutable=True)
            sage: G.dualflip(4)
            sage: G == G
            True

        Flip the edge where ``e`` is an end of the Strebel graph::

            sage: G = StrebelGraph("(0,~1,1,~2,2,~0)", mutable=True)
            sage: G.dualflip(0)
            sage: G1 = StrebelGraph("(0,~2,2,~1,1,~0)")
            sage: G == G1
            True
        """
        #
        # |                                  \
        # |                                   \
        # |  e                      e          \
        # v------------->w  ===> v------------->w
        #             E  |        \         E
        #                |         \
        #                |          \
        #
        #and
        #                          ------------
        #                         |            |
        #    e                    |            v
        # v -------> w      ===>  | v -------> w
        #            ^            |
        #            |            |
        #            |            \_____________
        #            |


        E = self._ep[e]
        a = self._fp[E]

        if a == e:
            self.halfdualflip(E)
            self.halfdualflip(e)
        else:
            self.halfdualflip(e)
            self.halfdualflip(E)

    def relhalfdualflip(self, e, S):

        #make S invariant under the edge permutation
        for b in S:
            B = self._ep[b]
            if B not in S:
                S.append(B)

        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        a = self._vp[e]

        if (a in S) and (self._bdry[e] == 0):
            pass
        else:
            self.halfdualflip(e)

    def fp_ordered(self, S):
        r"""
        We define a partial order on S andoutput a list of list of
        the half-edges in ``S``:

        - in each list, the half-edges are listed according to the
          partial order,
        - the half-edges from different lists are non-comparable

        Note: if a half-edge e is contained in S, then the method add the
        edge permutation of e in S automatically.

        EXAMPLES::

            sage: from veerer import *

            sage: G = StrebelGraph("(0,1,2)(~0:0)(~1:1,~2:1)", mutable=True)
            sage: S = [2,1]
            sage: G.fp_ordered(S)
            [[4], [3], [1, 2]]

        The following is an example when ``S`` contains an end of the graph::

            sage: G = StrebelGraph("(0,~1,1,~2,2,~0)", mutable=True)
            sage: S = [5]
            sage: G.fp_ordered(S)
            [[5, 0]]
        """

        S1 = S.copy()

        #make S invariant under the edge permutation
        for e in S1:
            E = self._ep[e]
            if E not in S1:
                S1.append(E)

        orderedlist = []


        while len(S1) > 0:

            e = S1.pop()
            l = [e]

            a0 = e
            angle0 = self._bdry[e]

            #extend half-edge e backwards
            while angle0 == 0:

                E0 = self._vp[a0]
                a0 = self._ep[E0]

                if a0 in S1:

                    S1.remove(a0)
                    l = [a0] + l

                    angle0 = self._bdry[a0]

                else: #if a0 = e0, then e0 is a loop surrounding a simple pole. in this case, l = [e0]
                    break

            a1 = self._fp[e]
            angle1 = self._bdry[a1]

            #extend half-edge e forwards
            while angle1 == 0:

                if a1 in S1:

                    S1.remove(a1)
                    l.append(a1)

                    a1 = self._fp[a1]
                    angle1 = self._bdry[a1]

                else:
                    break

            orderedlist.append(l)

        return orderedlist


    def multidualflip(self, S):
        r"""
        Flip several half-edges of the Strebel graph.

        Note: if a half-edge e is contained in S, then the method add the
        edge permutation of e in S automatically.

        EXAMPLES::

            sage: from veerer import *

            sage: G = StrebelGraph("(0,1,2)(~0:0)(~1:1,~2:1)", mutable=True)
            sage: S = [1,2]
            sage: G.multidualflip(S)
            sage: G1 = StrebelGraph("(~0,1,2)(0)(~1:1,~2:1)")
            sage: G == G1
            True

        The following is an example when ``S`` contains an end of the graph::

            sage: G = StrebelGraph("(0,~1,1,~2,2,~0)", mutable=True)
            sage: S = [5]
            sage: G.multidualflip(S)
            sage: G1 = StrebelGraph("(0,~2,2,~1,1,~0)")
            sage: G == G1
            True
        """

        ll = self.fp_ordered(S)

        for l in ll:
            for e in l:
                self.relhalfdualflip(e, S)

    def is_abelian(self, certificate=False):
        r"""
        return whether a strebel graph is abelian

        EXAMPLES::

        """

        ep = self._ep
        vp = self._vp

        # The code attempt to give a coherent holonomy with signs for each half
        # edge. For that purpose, it is enough to store a boolean for each edge:
        #   True: +
        #   False: -
        # In other words, the half edge is stored with True if its x-coordinate
        # is positive.
        #
        # The propagation of half-edge orientations is done by walking around
        # vertices

        oris = [None] * self._n  # list of orientations decided so far
        oris[0] = True
        oris[ep[0]] = False
        q = [0, ep[0]] # queue of half-edges to be treated

        while q:
            e = q.pop()
            o = oris[e]
            assert o is not None
            f = vp[e]
            while True:
                # propagate the orientation from e
                if self._bdry[e] % 2 == 0: #if even, we flip the edge
                    o = not o

                if oris[f] is None:
                    assert oris[ep[f]] is None
                    oris[f] = o
                    oris[ep[f]] = not o
                    q.append(ep[f])
                elif oris[f] != o:
                    return (False, None) if certificate else False
                else:
                    break

                e, f = f, vp[f]

        return (True, oris) if certificate else True

    def stratum(self):
        r"""
        Return the stratum of Abelian or quadratic differentials of this Strebel graph.

        EXAMPLES::

            sage: from veerer import StrebelGraph
        """
        # degrees of holomorphic part
        hol = [len(v) - 2 + sum(self._bdry[i] for i in v) for v in self.vertices()]

        # degrees of meromorphic part
        mer = [-2 - sum(self._bdry[i] for i in f) for f in self.faces()]

        # Determine whether it is Abelian (k=1) or quadratic (k=2) stratum
        if any(x % 2 for x in hol) or any(x % 2 for x in mer) or not self.is_abelian():
            k  = 2
        else:
            hol = [x // 2 for x in hol]
            mer = [x // 2 for x in mer]
            k = 1

        from surface_dynamics import Stratum
        return Stratum(hol + mer, k)

    # TODO: deprecate
    @staticmethod
    def from_face_edge_perms(vp, ep, fp=None, boundary=None, mutable=False, check=True):
        r"""
        INPUT:

        - ``fp``, ``ep``, ``vp`` -- the face, edge and vertex permutation

        - ``check`` - boolean (default: ``True``) - if set to ``False`` no
          check are performed

        EXAMPLES::

            sage: from veerer import *
            sage: from array import array
            sage: vp = array('i', [1, 3, 0, 2])
            sage: ep = array('i', [3, 2, 1, 0])
            sage: StrebelGraph.from_face_edge_perms(vp, ep, boundary = array('i', [1, 0, 0, 1]))
            doctest:warning
            ...
            UserWarning: the StrebelGraph.from_face_edge_perms is deprecated; use the classmethod from_permutations instead
            StrebelGraph("(0:1,1,~0:1,~1)")
        """
        import warnings
        warnings.warn('the StrebelGraph.from_face_edge_perms is deprecated; use the classmethod from_permutations instead')

        n = len(vp)

        if fp is None:
            fp = array('i', [-1] * n)
            for i in range(n):
                fp[ep[vp[i]]] = i

        if boundary is None:
            bdry = array('i', [0] * n)
        else:
            bdry = array('i', boundary)

        return StrebelGraph.from_permutations(vp, ep, fp, (bdry,), mutable, check)
