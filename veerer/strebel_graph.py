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

from .permutation import perm_cycles, perm_cycles_to_string, str_to_cycles_and_data
from .triangulation import face_edge_perms_init, boundary_init
from .constellation import Constellation
from .constants import *


def one_edge_completion(T):
        r'''
        Return two colored triangulations (represented by a list of attributes) that adds one edge to T
        with different colors.

        T is a list: [face permutation, vertex permutation, colouring, boundary, bdryedge],
        where bdryedge is a list consisting of ``0`` or ``1`` and bdryedge[e] = 1 if and only
        if e is a boundary edge of T
        '''

        fp, vp, colouring, boundary, bdryedge = T
        n = len(bdryedge) #current number of the half-edges
        m = n // 2

        for e in range(n):
            if bdryedge[e]: #if e is a boundary half-edge
                ang = boundary[e] #angle excess of the half-edge
                if ang == 0:
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
                    newbdryedge.insert(m, 0) #the half-edge n/2 is always internal
                    newbdryedge.insert(m+1, 1) #the half-edge n/2 + 1 is always boundary

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
                    newfp[newep[c]] = m + 1
                    newfp[m + 1] = b

                    newvp[newep[a]] = m + 1
                    newvp[m + 1] = c
                    newvp[b] = m
                    newvp[m] = newep[e]

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

                    if ((newcolouring_1[newep[a]], newcolouring_1[m+1], newcolouring_1[c]) == (BLUE, RED, BLUE)) or ((newcolouring_1[newep[a]], newcolouring_1[m+1], newcolouring_1[c]) == (RED, BLUE, RED)):
                        newboundary_1[m+1] = newboundary_1[newep[a]] - 1
                        newboundary_1[newep[a]] = 0
                        newboundary_1[m] = 0
                    elif ((newcolouring_1[newep[e]], newcolouring_1[m], newcolouring_1[b]) == (RED, BLUE, RED)) or ((newcolouring_1[newep[e]], newcolouring_1[m], newcolouring_1[b]) == (BLUE, RED, BLUE)):
                        newboundary_1[m+1] = newboundary_1[newep[a]]
                        newboundary_1[newep[a]] = 0
                        newboundary_1[b] = newboundary_1[b] - 1
                        newboundary_1[m] = 0
                    else:
                        newboundary_1[m+1] = newboundary_1[newep[a]]
                        newboundary_1[newep[a]] = 0
                        newboundary_1[m] = 0

                    if ((newcolouring_2[newep[a]], newcolouring_2[m+1], newcolouring_2[c]) == (BLUE, RED, BLUE)) or ((newcolouring_2[newep[a]], newcolouring_2[m+1], newcolouring_2[c]) == (RED, BLUE, RED)):
                        newboundary_2[m+1] = newboundary_2[newep[a]] - 1
                        newboundary_2[newep[a]] = 0
                        newboundary_2[m] = 0
                    elif ((newcolouring_2[newep[e]], newcolouring_2[m], newcolouring_2[b]) == (RED, BLUE, RED)) or ((newcolouring_2[newep[e]], newcolouring_2[m], newcolouring_2[b]) == (BLUE, RED, BLUE)):
                        newboundary_2[m+1] = newboundary_2[newep[a]]
                        newboundary_2[newep[a]] = 0
                        newboundary_2[b] = newboundary_2[b] - 1
                        newboundary_2[m] = 0
                    else:
                        newboundary_2[m+1] = newboundary_2[newep[a]]
                        newboundary_2[newep[a]] = 0
                        newboundary_2[m] = 0


                    T1 = [newfp, newvp, newcolouring_1, newboundary_1, newbdryedge]
                    T2 = [newfp, newvp, newcolouring_2, newboundary_2, newbdryedge]

                    return [T1, T2]

        raise ValueError("T is already complete")


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
        Return the associated angle excess of the corners of a coloured strebel graph
        """

        #remark: red-red corners with angle excess 0 and 1 both correspond
        # to zero half-plane excess.
        #claim: it is impossible to have a red-red corner with 0 angle
        # excess in a strebel graph.

        n = self._n
        vp = self._vp
        beta = self._bdry
        alpha = array('i', [-1] * n)

        for e in range(n):
            e1 = vp[e]
            if slope == VERTICAL:
                if colouring[e] == RED and colouring[e1] == BLUE:
                    alpha[e] = beta[e]
                else:
                    alpha[e] = beta[e] + 1
            elif slope == HORIZONTAL:
                if colouring[e] == BLUE and colouring[e1] == RED:
                    alpha[e] = beta[e]
                else:
                    alpha[e] = beta[e] + 1

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

    def veering_triangulations(self, slope=VERTICAL, inclusion=False):
        r"""
        return a list of pairs: (veering, inclusion), where:
        ``veering`` is a slope-Strebel veering triangulation associated to the given Strebel graph,
        and ``inclusion`` is the map from the half-edges of the strebel graph into the half-edges of ``veering``

        EXAMPLES::

            sage: from veerer import *
            sage: G = StrebelGraph("(~1:1,~0,1:1,0)")
            sage: G.stratum()  # optional - surface_dynamics
            H_1(2, -2)
            sage: G.veering_triangulations(inclusion=True)
            [(VeeringTriangulation("", boundary="(0:1,~1:2,~0:1,1:2)", colouring="RR"),
              [0, 1, 2, 3]),
             (VeeringTriangulation("", boundary="(0:1,~1:1,~0:1,1:1)", colouring="RB"),
              [2, 0, 3, 1]),
             (VeeringTriangulation("", boundary="(0:1,~1:2,~0:1,1:2)", colouring="BB"),
              [0, 1, 2, 3]),
             (VeeringTriangulation("(0,~3,2)(1,3,~2)", boundary="(~1:2,~0:2)", colouring="BBRB"),
              [2, 4, 3, 5]),
             (VeeringTriangulation("(0,~3,2)(1,3,~2)", boundary="(~1:2,~0:2)", colouring="RRRB"),
              [2, 4, 3, 5])]
            sage: len(set(G.veering_triangulations(slope=HORIZONTAL)))
            5
        """
        from .veering_triangulation import VeeringTriangulation

        def is_complete(T):
            fp, vp, colouring, boundary, bdryedge = T
            return not any(b1 and b2 == 0 for (b1, b2) in zip(bdryedge, boundary))

        n = self._n
        ep_st = self._ep

        complete = []
        incomplete = []
        bdryedge =  [1] * n
        for colouring in self.colourings():
            boundary = self.angle_excess(colouring, slope=slope)
            T = (self._fp, self._vp, colouring, boundary, bdryedge)
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

        veering = []
        for fp, vp, colouring, boundary, bdryedge in complete:

            nn = len(fp)
            ep = array('i', list(reversed(range(nn))))

            #assert False, (ep,fp,boundary,colouring)
            T = VeeringTriangulation.from_permutations(vp, ep, fp, (boundary, colouring), mutable=True)
            if T.is_delaunay():
                veering.append(T)

        #keep different veering triangulations by considering canonical labels,
        # and record the map from the half-edges of the strebel graph into veering triangulations
        if not inclusion:
            for vt in veering:
                vt.set_canonical_labels()
                vt.set_immutable()
            return list(set(veering))
        else:
            inclusion_veering = {}
            for vt in veering:
                r, _ = vt.best_relabelling()
                vt.relabel(r, check=False)
                vt.set_immutable()
                if vt not in inclusion_veering:
                    #build the inclusion
                    ep = vt._ep
                    inclusion = []
                    for e in range(n):
                        if e < ep_st[e]:
                            inclusion.append(r[e])
                        else:
                            inclusion.append(ep[r[ep_st[e]]])
                    inclusion_veering[vt] = inclusion
            return list(inclusion_veering.items())
