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

from array import array
import itertools
import numbers

from sage.rings.all import ZZ

from .permutation import perm_check, perm_cycles, perm_cycles_to_string, str_to_cycles_and_data
from .triangulation import face_edge_perms_init, boundary_init
from .constellation import Constellation
from .constants import *


def one_edge_completion(T):
    r"""
    Return two colored triangulations (represented by a list of attributes) that adds one edge to T
    with different colors.

    T is a list: [vertex permutation, edge permutation, face permutation, colouring, boundary, bdryedge],
    where bdryedge is a list consisting of ``0`` or ``1`` and bdryedge[e] = 1 if and only
    if e is a boundary edge of T

    EXAMPLES::

        sage: from array import array
        sage: from veerer.strebel_graph import one_edge_completion
        sage: vp = array('i', [3, 2, 0, 1])
        sage: ep = array('i', [3, 2, 1, 0])
        sage: fp = array('i', [0, 1, 3, 2])
        sage: c = [2,1,1,2]
        sage: bdry = array('i', [3, 3, 0, 3])
        sage: T = [vp,ep,fp,c,bdry,[1]*4]
        sage: one_edge_completion(T)
        [[array('i', [5, 4, 1, 2, 0, 3]),
        array('i', [5, 4, 3, 2, 1, 0]),
        array('i', [0, 1, 5, 3, 2, 4]),
        array('i', [2, 1, 1, 1, 1, 2]),
        array('i', [3, 3, 0, 3, 0, 0]),
        [1, 1, 0, 1, 0, 0]],
        [array('i', [5, 4, 1, 2, 0, 3]),
        array('i', [5, 4, 3, 2, 1, 0]),
        array('i', [0, 1, 5, 3, 2, 4]),
        array('i', [2, 1, 2, 2, 1, 2]),
        array('i', [3, 3, 0, 3, 0, 0]),
        [1, 1, 0, 1, 0, 0]]]
    """
    vp, ep, fp, colouring, boundary, bdryedge = T
    n = len(bdryedge) #current number of the half-edges
    m = n // 2

    found = False
    for e in range(n):
        if bdryedge[e] and boundary[e] == 0:
            found = True
            break
    if not found:
        raise ValueError('already complete')

    #new edge permutation
    newep = array('i', list(reversed(range(n + 2))))

    #change the labels of e, vp, fp, colouring and boundary due to the added edge
    if e > m - 1:
        e = e + 2

    newfp = array('i', fp)
    newvp = array('i', vp)
    newcolouring = array('i', colouring)
    newboundary = boundary[:]
    newbdryedge = bdryedge[:]

    newfp.insert(m, -1)
    newfp.insert(m+1, -1)
    newvp.insert(m, -1)
    newvp.insert(m+1, -1)
    newcolouring.insert(m, -1)
    newcolouring.insert(m+1, -1)
    newboundary.insert(m, -1)
    newboundary.insert(m+1, -1)
    newbdryedge.insert(m, 0) # the half-edge m is always internal
    newbdryedge.insert(m+1, 1) #the half-edge m + 1 is always boundary

    for ee in range(n + 2):
        if (ee != m) and (ee != m + 1):
            ve = newvp[ee]
            fe = newfp[ee]
            if fe > m - 1:
                newfp[ee] = fe + 2
            if ve > m - 1:
                newvp[ee] = ve + 2

    #build the new triangle face
    a = newvp[e]
    b = newfp[e]
    c = newvp[newep[a]]
    newfp[e] = m
    newfp[m] = newep[a]

    #build the vertex permutation for half-edges of the triangle face
    newvp[newep[a]] = m + 1
    newvp[m] = newep[e]

    #build the new boundary face and the vertex permutation for half-edges in the new boundary face
    if c == newep[e]:
        assert b == newep[a]
        newfp[m + 1] = m + 1
        newvp[m + 1] = m
        newvp[b] = m + 1
    else:
        newfp[newep[c]] = m + 1
        newfp[m + 1] = b
        newvp[m + 1] = c
        newvp[b] = m

    #modify the bdryedge
    assert newbdryedge[e] == 1
    newbdryedge[e] = 0 #e becomes internal

    assert newbdryedge[newep[a]] == 1
    newbdryedge[newep[a]] = 0 #newep[a] becomes internal

    #build new colouring
    #claim: for e1=vp[e], e and e1 must have the same colour
    newcolouring_1 = array('i', newcolouring)
    newcolouring_2 = array('i', newcolouring)

    newcolouring_1[m] = newcolouring_1[m + 1] = 1 #RED
    newcolouring_2[m] = newcolouring_2[m + 1] = 2 #BLUE

    #build new angle excess
    newboundary_1 = newboundary[:]
    newboundary_2 = newboundary[:]

    #new internal edge
    newboundary_1[newep[a]] = newboundary_2[newep[a]] = 0
    newboundary_1[m] = newboundary_2[m] = 0

    #determin colors of ``m+1`` and ``newep[a]``
    if b == newep[a]:
        newboundary_1[m+1] = newboundary[newep[a]]
    elif ((newcolouring_1[newep[a]], newcolouring_1[m+1], newcolouring_1[c]) == (BLUE, RED, BLUE)) or ((newcolouring_1[newep[a]], newcolouring_1[m+1], newcolouring_1[c]) == (RED, BLUE, RED)):
        newboundary_1[m+1] = newboundary[newep[a]] - 1
    elif ((newcolouring_1[newep[e]], newcolouring_1[m], newcolouring_1[b]) == (RED, BLUE, RED)) or ((newcolouring_1[newep[e]], newcolouring_1[m], newcolouring_1[b]) == (BLUE, RED, BLUE)):
        newboundary_1[m+1] = newboundary[newep[a]]
        newboundary_1[b] = newboundary[b] - 1
    else:
        newboundary_1[m+1] = newboundary[newep[a]]

    if b == newep[a]:
        newboundary_2[m+1] = newboundary[newep[a]]
    elif ((newcolouring_2[newep[a]], newcolouring_2[m+1], newcolouring_2[c]) == (BLUE, RED, BLUE)) or ((newcolouring_2[newep[a]], newcolouring_2[m+1], newcolouring_2[c]) == (RED, BLUE, RED)):
        newboundary_2[m+1] = newboundary[newep[a]] - 1
    elif ((newcolouring_2[newep[e]], newcolouring_2[m], newcolouring_2[b]) == (RED, BLUE, RED)) or ((newcolouring_2[newep[e]], newcolouring_2[m], newcolouring_2[b]) == (BLUE, RED, BLUE)):
        newboundary_2[m+1] = newboundary[newep[a]]
        newboundary_2[b] = newboundary[b] - 1
    else:
        newboundary_2[m+1] = newboundary[newep[a]]

    assert perm_check(newvp), (T, newvp)
    assert perm_check(newep), (T, newep)
    assert perm_check(newfp), (T, newfp)

    # TODO: here we use the same arrays newvp, newep, newfp for two distinct
    # data. This can lead to two veering triangulations sharing data which
    # must be avoided when mutable=True.
    T1 = [newvp, newep, newfp, newcolouring_1, newboundary_1, newbdryedge]
    T2 = [newvp, newep, newfp, newcolouring_2, newboundary_2, newbdryedge]

    return [T1, T2]


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

    def base_ring(self):
        from sage.rings.integer_ring import ZZ
        return ZZ

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
            UserWarning: the method StrebelGraph.from_face_edge_perms is deprecated; use the classmethod from_permutations instead
            StrebelGraph("(0:1,1,~0:1,~1)")
        """
        import warnings
        warnings.warn('the method StrebelGraph.from_face_edge_perms is deprecated; use the classmethod from_permutations instead')

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

    def as_linear_family(self):
        r"""
        Return this Strebel graph as a linear family.

        EXAMPLES::

            sage: from veerer import StrebelGraph
            sage: G = StrebelGraph("(0,1,2)(~0,~1:1,~2:2)")
            sage: G.as_linear_family()
            StrebelGraphLinearFamily("(0,1,2)(~2:2,~0,~1:1)", [(1, 0, 0), (0, 1, 0), (0, 0, 1)])
        """
        from sage.matrix.special import identity_matrix
        from .linear_family import StrebelGraphLinearFamily
        return StrebelGraphLinearFamily(self, identity_matrix(ZZ, self.num_edges()))

    def angle_excess(self, colouring, slope=VERTICAL):
        r"""
        Return the angle excess of the corners of the associated colouring.

        EXAMPLES::

            sage: from veerer import StrebelGraph
            sage: G = StrebelGraph("(0,1,2)(~0,~1:1,~2:2)")
            sage: for colouring in G.colourings():
            ....:     print(colouring, G.angle_excess(colouring))
            array('i', [1, 1, 1, 1, 1, 1]) array('i', [1, 1, 1, 3, 2, 1])
            array('i', [1, 1, 2, 2, 1, 1]) array('i', [0, 1, 1, 3, 2, 0])
            array('i', [1, 2, 1, 1, 2, 1]) array('i', [1, 1, 0, 2, 2, 1])
            array('i', [1, 2, 2, 2, 2, 1]) array('i', [0, 1, 1, 3, 2, 0])
            array('i', [2, 1, 1, 1, 1, 2]) array('i', [1, 0, 1, 3, 1, 1])
            array('i', [2, 1, 2, 2, 1, 2]) array('i', [1, 0, 1, 3, 1, 1])
            array('i', [2, 2, 1, 1, 2, 2]) array('i', [1, 1, 0, 2, 2, 1])
            array('i', [2, 2, 2, 2, 2, 2]) array('i', [1, 1, 1, 3, 2, 1])
        """
        # remark: red-red corners with angle excess 0 and 1 both correspond
        # to zero half-plane excess.
        # claim: it is impossible to have a red-red corner with 0 angle
        # excess in a strebel graph.

        n = self._n
        vp = self._vp
        alpha = array('i', self._bdry)

        for e in range(n):
            e1 = vp[e]
            if (slope == VERTICAL and (colouring[e] != RED or colouring[e1] != BLUE)) or \
               (slope == HORIZONTAL and (colouring[e] != BLUE or colouring[e1] != RED)):
                alpha[e] += 1

        return alpha

    def colourings(self):
        r"""
        Run through the red, blue colourings of this Strebel graph.

        Each colouring consists of an array of length the number of half-edges.
        """
        n = self._n
        m = n // 2
        ep = self._ep

        for colouring in itertools.product([RED, BLUE], repeat=m):
            colouring = array('i', colouring)
            colouring.extend([colouring[ep[e]] for e in range(m, n)])
            yield colouring

    def veering_triangulations(self, colouring, slope=VERTICAL, mutable=False):
        r"""
        Run through Strebel-veering triangulations obtained by completing this
        Strebel graph given the ``colouring`` of its edges.

        EXAMPLES::

            sage: from veerer import *
            sage: G = StrebelGraph("(~1:1,~0,1:1,0)")
            sage: G.stratum()  # optional - surface_dynamics
            H_1(2, -2)
            sage: for colouring in G.colourings():
            ....:     vts = G.veering_triangulations(colouring)
            ....:     assert all(vt.strebel_graph() == G for vt in vts)
            ....:     assert all(vt.stratum() == G.stratum() for vt in vts)  # optional - surface_dynamics
            ....:     print(colouring, list(G.veering_triangulations(colouring)))
            array('i', [1, 1, 1, 1]) [VeeringTriangulation("", boundary="(0:1,~1:2,~0:1,1:2)", colouring="RR")]
            array('i', [1, 2, 2, 1]) [VeeringTriangulation("(0,2,1)(3,~1,~0)", boundary="(~3:1,~2:2)", colouring="RBBR"), VeeringTriangulation("(0,2,1)(3,~1,~0)", boundary="(~3:2,~2:2)", colouring="RBBB"), VeeringTriangulation("(0,2,1)(3,~1,~0)", boundary="(~3:2,~2:2)", colouring="RBRR"), VeeringTriangulation("(0,2,1)(3,~1,~0)", boundary="(~3:2,~2:1)", colouring="RBRB")]
            array('i', [2, 1, 1, 2]) [VeeringTriangulation("", boundary="(0:1,~1:1,~0:1,1:1)", colouring="BR")]
            array('i', [2, 2, 2, 2]) [VeeringTriangulation("", boundary="(0:1,~1:2,~0:1,1:2)", colouring="BB")]
        """
        def is_complete(T):
            vp, ep, fp, colouring, boundary, bdryedge = T
            return not any(b1 and b2 == 0 for (b1, b2) in zip(bdryedge, boundary))

        n = self._n
        bdryedge =  [1] * n
        veering = []
        boundary = self.angle_excess(colouring, slope=slope)
        T = (self._vp[:], self._ep[:], self._fp[:], colouring, boundary, bdryedge)
        complete = []
        incomplete = []
        if is_complete(T):
            complete.append(T)
        else:
            incomplete.append(T)

        while incomplete:
            T = incomplete.pop()
            for TT in one_edge_completion(T):
                if is_complete(TT):
                    complete.append(TT)
                else:
                    incomplete.append(TT)

        from .veering_triangulation import VeeringTriangulation
        for vp, ep, fp, colouring, boundary, bdryedge in complete:
            nn = len(fp)
            yield VeeringTriangulation.from_permutations(vp, ep, fp, (boundary, colouring), mutable=mutable, check=True)
