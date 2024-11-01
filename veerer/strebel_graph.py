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
from .triangulation import face_edge_perms_init, boundary_init, Triangulation
from .constellation import Constellation
from .constants import *
from veerer.polyhedron import *


def one_edge_completion(t, angle_excess, colouring):
    r"""
    Return a pair of triples ``(new_t, new_angle_excess, new_colouring)`` corresponding to
    adding an edge forming a triangle with the two allowed colours.

    INPUT:

    - ``t`` -- a :class:`Triangulation` (with boundary)

    - ``angle_excess`` -- an array of angle excess

    - ``colouring`` -- an array of RED, BLUE colouring

    EXAMPLES::

        sage: from array import array
        sage: from veerer import Triangulation
        sage: from veerer.strebel_graph import one_edge_completion
        sage: t = Triangulation("", "(0:1)(1:1)(~0:1,~1:1)")
        sage: colouring = [2,1,1,2]
        sage: angle_excess = array('i', [3, 3, 0, 3])
        sage: one_edge_completion(t, angle_excess, colouring)
        ((Triangulation("(2,~0,~1)", boundary="(0:1)(1:1)(~2:1)"),
          array('i', [3, 3, 0, 3, 0, 0]),
          array('i', [2, 1, 1, 1, 1, 2])),
         (Triangulation("(2,~0,~1)", boundary="(0:1)(1:1)(~2:1)"),
          array('i', [3, 3, 0, 3, 0, 0]),
          array('i', [2, 1, 2, 2, 1, 2])))
    """
    # We insert the (m, M)-edge as follows
    #
    #   x
    #   |
    #   |
    #   |c
    #   |
    #   x
    #   | \ M
    #   |A   \
    #   |       \
    #  a|  e      m\   b
    #   x-----------x-----------x
    #             E           B

    vp = t._vp
    ep = t._ep
    fp = t._fp
    bdry = t._bdry
    angle_excess = angle_excess
    n = t._n
    m = n // 2
    M = m + 1

    found = False
    for e in range(n):
        if bdry[e] and angle_excess[e] == 0:
            found = True
            break
    if not found:
        raise ValueError('already complete')


    # add an edge (m, M) and possibly shift e
    newvp = array('i', vp)
    newep = array('i', ep)
    newfp = array('i', fp)

    newvp.insert(m, -1)
    newvp.insert(M, -1)
    newep.insert(m, -1)
    newep.insert(M, -1)
    newfp.insert(m, -1)
    newfp.insert(M, -1)

    newep[m] = M
    newep[M] = m

    for ee in range(n + 2):
        if (ee != m) and (ee != M):
            vpe = newvp[ee]
            fpe = newfp[ee]
            epe = newep[ee]
            if fpe > m - 1:
                newfp[ee] = fpe + 2
            if vpe > m - 1:
                newvp[ee] = vpe + 2
            if epe > m - 1:
                newep[ee] = epe + 2

    if e > m - 1:
        e = e + 2
    E = newep[e]

    # build the new triangle face (e, m, A)
    a = newvp[e]
    A = newep[a]
    b = newfp[e]
    c = newvp[A]
    C = newep[c]
    newfp[e] = m
    newfp[m] = A
    newvp[A] = M
    newvp[m] = E

    is_bigon = c == E
    assert is_bigon == (b == A)

    if is_bigon:
        assert b == A
        newfp[M] = M
        newvp[M] = m
        newvp[b] = M
    else:
        newfp[C] = M
        newfp[M] = b
        newvp[M] = c
        newvp[b] = m

    # new bdry: m is internal, M is boundary and e and A becomes internal
    newbdry = bdry[:]
    newbdry.insert(m, 0)
    newbdry.insert(M, 1)

    assert newbdry[e] == 1
    assert newbdry[A] == 1
    newbdry[e] = newbdry[A] = 0

    # build new colourings
    colouring1 = array('i', colouring)
    colouring1.insert(m, -1)
    colouring1.insert(M, -1)

    # claim: a and e must have different colours
    # (so that both RED and BLUE are allowed for the new edge (m, M)
    assert colouring1[a] != colouring1[e]

    colouring2 = colouring1[:]
    colouring1[m] = colouring1[M] = RED
    colouring2[m] = colouring2[M] = BLUE

    # build new angle excesses
    angle_excess1 = angle_excess[:]
    angle_excess1.insert(m, -1)
    angle_excess1.insert(M, -1)
    angle_excess1[M] = angle_excess1[A]
    angle_excess1[A] = angle_excess1[m] = 0

    angle_excess2 = angle_excess1[:]

    if not is_bigon:
        col_a = colouring1[a]
        col_b = colouring1[b]
        col_c = colouring1[c]
        col_e = colouring1[e]

        if col_a == BLUE:
            assert col_e == RED
            if col_b == RED:
                angle_excess2[b] -= 1
            if col_c == BLUE:
                angle_excess1[M] -= 1
        else:
            assert col_a == RED and col_e == BLUE
            if col_b == BLUE:
                angle_excess1[b] -= 1
            if col_c == RED:
                angle_excess2[M] -= 1

    t = Triangulation.from_permutations(newvp, newep, newfp, (newbdry,), mutable=False, check=False)
    return ((t, angle_excess1, colouring1), (t, angle_excess2, colouring2))


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
    __slots__ = ['_excess']

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
        self._excess = self._data[0]

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
        bdry_cycles = perm_cycles_to_string(perm_cycles(self._fp, n=self._n), involution=self._ep, data=self._excess)
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
                if self._excess[e] % 2 == 0: #if even, we flip the edge
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
        hol = [len(v) - 2 + sum(self._excess[i] for i in v) for v in self.vertices()]

        # degrees of meromorphic part
        mer = [-2 - sum(self._excess[i] for i in f) for f in self.faces()]

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
        alpha = array('i', self._excess)

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

    def cone(self, backend=None):
        r"""
        EXAMPLES::

            sage: from veerer import StrebelGraph
            sage: sg = StrebelGraph("(~1:1,~0,1:1,0)")
            sage: sg.cone(backend='ppl')
            Cone of dimension 4 in ambient dimension 4 made of 2 facets (backend=ppl)
            sage: sg.cone(backend='sage')
            Cone of dimension 4 in ambient dimension 4 made of 2 facets (backend=sage)
            sage: sg.cone(backend='normaliz-QQ')  # optional - pynormaliz
            Cone of dimension 4 in ambient dimension 4 made of 2 facets (backend=normaliz-QQ)
        """
        from .polyhedron import LinearExpressions, ConstraintSystem
        L = LinearExpressions(self.base_ring())
        ne = self.num_edges()
        cs = ConstraintSystem(2 * ne)
        for e in range(ne):
            cs.insert(L.variable(e) >= 0)
        return cs.cone(backend=backend)

    def veering_triangulations(self, colouring, slope=VERTICAL, mutable=False):
        r"""
        Run through Strebel-veering triangulations obtained by completing this
        Strebel graph given the ``colouring`` of its edges.

        EXAMPLES::

            sage: from veerer import *
            sage: examples = []
            sage: examples.append(StrebelGraph("(~1:1,~0,1:1,0)"))
            sage: examples.append(StrebelGraph("(0,~1)(1)(~0)"))
            sage: examples.append(StrebelGraph("(0:2)(1:2)(~1,~0:2)"))
            sage: for G in examples:  # optional - surface_dynamics
            ....:     print(G.stratum())
            H_1(2, -2)
            H_0(1, -1^3)
            H_0(4, -2^3)

            sage: for G in examples:
            ....:     print(G)
            ....:     for colouring in G.colourings():
            ....:         vts = G.veering_triangulations(colouring)
            ....:         assert all(vt.strebel_graph() == G for vt in vts)
            ....:         assert all(vt.stratum() == G.stratum() for vt in vts)  # optional - surface_dynamics
            ....:         print(colouring, list(G.veering_triangulations(colouring)))
            StrebelGraph("(0,~1:1,~0,1:1)")
            array('i', [1, 1, 1, 1]) [VeeringTriangulation("", boundary="(0:1,~1:2,~0:1,1:2)", colouring="RR")]
            array('i', [1, 2, 2, 1]) [VeeringTriangulation("(0,2,1)(3,~1,~0)", boundary="(~3:1,~2:2)", colouring="RBBR"), VeeringTriangulation("(0,2,1)(3,~1,~0)", boundary="(~3:2,~2:2)", colouring="RBBB"), VeeringTriangulation("(0,2,1)(3,~1,~0)", boundary="(~3:2,~2:2)", colouring="RBRR"), VeeringTriangulation("(0,2,1)(3,~1,~0)", boundary="(~3:2,~2:1)", colouring="RBRB")]
            array('i', [2, 1, 1, 2]) [VeeringTriangulation("", boundary="(0:1,~1:1,~0:1,1:1)", colouring="BR")]
            array('i', [2, 2, 2, 2]) [VeeringTriangulation("", boundary="(0:1,~1:2,~0:1,1:2)", colouring="BB")]
            StrebelGraph("(0,~1)(1)(~0)")
            array('i', [1, 1, 1, 1]) [VeeringTriangulation("", boundary="(0:1,~1:1)(1:1)(~0:1)", colouring="RR")]
            array('i', [1, 2, 2, 1]) [VeeringTriangulation("(0,2,~1)", boundary="(1:1)(~2:1)(~0:1)", colouring="RBR"), VeeringTriangulation("(0,2,~1)", boundary="(1:1)(~2:1)(~0:1)", colouring="RBB")]
            array('i', [2, 1, 1, 2]) [VeeringTriangulation("(0,~1,2)", boundary="(1:1)(~2:1)(~0:1)", colouring="BRR"), VeeringTriangulation("(0,~1,2)", boundary="(1:1)(~2:1)(~0:1)", colouring="BRB")]
            array('i', [2, 2, 2, 2]) [VeeringTriangulation("", boundary="(0:1,~1:1)(1:1)(~0:1)", colouring="BB")]
            StrebelGraph("(0:2)(1:2)(~1,~0:2)")
            array('i', [1, 1, 1, 1]) [VeeringTriangulation("", boundary="(0:3)(1:3)(~1:1,~0:3)", colouring="RR")]
            array('i', [1, 2, 2, 1]) [VeeringTriangulation("", boundary="(0:3)(1:3)(~1:1,~0:2)", colouring="RB")]
            array('i', [2, 1, 1, 2]) [VeeringTriangulation("(2,~0,~1)", boundary="(0:3)(1:3)(~2:3)", colouring="BRR"), VeeringTriangulation("(2,~0,~1)", boundary="(0:3)(1:3)(~2:3)", colouring="BRB")]
            array('i', [2, 2, 2, 2]) [VeeringTriangulation("", boundary="(0:3)(1:3)(~1:1,~0:3)", colouring="BB")]
        """
        def is_complete(t, angle_excess, colouring):
            return not any(b1 and b2 == 0 for (b1, b2) in zip(t._bdry, angle_excess))

        n = self._n
        angle_excess = self.angle_excess(colouring, slope=slope)
        t0 = Triangulation.from_permutations(self._vp[:], self._ep[:], self._fp[:], (array('i', [1] * n),), mutable=False, check=True)
        T = (t0, angle_excess, colouring)
        complete = []
        incomplete = []
        if is_complete(*T):
            complete.append(T)
        else:
            incomplete.append(T)

        while incomplete:
            T = incomplete.pop()
            for TT in one_edge_completion(*T):
                if is_complete(*TT):
                    complete.append(TT)
                else:
                    incomplete.append(TT)

        from .veering_triangulation import VeeringTriangulation
        for t, angle_excess, colouring in complete:
            vp = t._vp
            ep = t._ep
            fp = t._fp
            cols = array('i', colouring)
            yield VeeringTriangulation.from_permutations(vp, ep, fp, (angle_excess, cols), mutable=mutable, check=True)

    def delaunay_triangulations(self, colouring, slope=VERTICAL, mutable=False, backend=None):
        r"""
        Run through the Delaunay Strebel-veering triangulations

        EXAMPLES::

            sage: from veerer import *
            sage: examples = []
            sage: examples.append(StrebelGraph("(~1:1,~0,1:1,0)"))
            sage: examples.append(StrebelGraph("(0,2,~1)(1)(~2,~0)"))
            sage: examples.append(StrebelGraph("(0:2,2,~1)(1,~0)(~2)"))
            sage: for G in examples:
            ....:     print(G)
            ....:     for colouring in G.colourings():
            ....:         print(colouring, sum(1 for _ in G.veering_triangulations(colouring)), sum(1 for _ in G.delaunay_triangulations(colouring)))
            StrebelGraph("(0,~1:1,~0,1:1)")
            array('i', [1, 1, 1, 1]) 1 1
            array('i', [1, 2, 2, 1]) 4 2
            array('i', [2, 1, 1, 2]) 1 1
            array('i', [2, 2, 2, 2]) 1 1
            StrebelGraph("(0,2,~1)(1)(~2,~0)")
            array('i', [1, 1, 1, 1, 1, 1]) 1 1
            array('i', [1, 1, 2, 2, 1, 1]) 6 5
            array('i', [1, 2, 1, 1, 2, 1]) 3 3
            array('i', [1, 2, 2, 2, 2, 1]) 6 5
            array('i', [2, 1, 1, 1, 1, 2]) 6 3
            array('i', [2, 1, 2, 2, 1, 2]) 3 3
            array('i', [2, 2, 1, 1, 2, 2]) 6 3
            array('i', [2, 2, 2, 2, 2, 2]) 1 1
            StrebelGraph("(0:2,2,~1)(1,~0)(~2)")
            array('i', [1, 1, 1, 1, 1, 1]) 1 1
            array('i', [1, 1, 2, 2, 1, 1]) 2 2
            array('i', [1, 2, 1, 1, 2, 1]) 2 2
            array('i', [1, 2, 2, 2, 2, 1]) 2 2
            array('i', [2, 1, 1, 1, 1, 2]) 6 5
            array('i', [2, 1, 2, 2, 1, 2]) 6 5
            array('i', [2, 2, 1, 1, 2, 2]) 2 2
            array('i', [2, 2, 2, 2, 2, 2]) 1 1
        """
        for vt in self.veering_triangulations(colouring, slope, mutable):
            if vt.is_delaunay(backend):
                yield vt

    def _set_strebel_constraints(self, cs, L, slope):
    
        nhe = self._n
        ne = ZZ(nhe // 2)
        
        x = [L.variable(i) for i in range(nhe)]
        
        if slope == VERTICAL:
            for i in range(ne):
                cs.insert(x[i] >= 0)
        elif slope == HORIZONTAL:
            for j in range(ne, nhe):
                cs.insert(x[j] >= 0)
    
    def _set_residue_constraints(self, cs, coes, L, slope):
        r"""
        Return the constraints on the residues.

        The ``coes`` is a list of lists of coefficient. The numbers in each list of coefficients corresponds to the boundary faces of the Strebel graph.
        """

        ep = self._ep
        list_faces = self.faces()
        nf = len(list_faces)
        nhe = self._n
        ne = ZZ(nhe // 2)
        
        x = [L.variable(i) for i in range(nhe)]
        
        if self.is_abelian():
            o = [1 if e else -1 for e in self.is_abelian(certificate=True)[1]]
        else:
            raise ValueError('{} is not Abelian'.format(self))
            
        #build expression of residues
        p = []
        
        if slope == VERTICAL:
            for f in list_faces:
                pp = 0
                for e in f:
                    if e < ep[e]:
                        pp = pp + o[e]*x[e]
                    else:
                        pp = pp + o[e]*x[ep[e]]
                p.append(pp)
        elif slope == HORIZONTAL:
            for f in list_faces:
                pp = 0
                for e in f:
                    if e < ep[e]:
                        pp = pp + o[e]*x[ep[e]]
                    else:
                        pp = pp + o[e]*x[e]
                p.append(pp)
        else:
            raise ValueError('bad slope parameter')
            
        #build constraints    
        for coe in coes:
            assert nf == len(coe)
            res = sum(coe[i]*p[i] for i in range(nf))
            if not res.is_zero():
                cs.insert(res == 0)
        
        d = cs.cone().affine_dimension() #the contraints might be redundant
        if d == 0:
            raise ValueError('the constraints on the residues of {} have zero-dimensional solution space'.format(self))

    def residue_constraint_cone(self, coes, K, slope):
            r"""
            return a list generators for the subspace with constraints on residues.
            
            Question: Are the generators irralevant to the horizontal or vertical Strebel decomposition considered?

            EXAMPLES::
            
                sage: from veerer import *
                sage: from sage.rings.rational_field import QQ

            Here are examples that residue constraints force some coordinate to be zero::
                sage: G0 = StrebelGraph("(0)(~0)")
                sage: G0.residue_constraint_cone([[1, -1]], QQ, VERTICAL)
                Traceback (most recent call last):
                ...
                ValueError: the constraints on the residues force some coordinates to be zero 
                in the Strebel polyhedron of StrebelGraph("(0:0)(~0:0)")
                sage: G1 = StrebelGraph("(0,~1)(1)(~0)")
                sage: G1.residue_constraint_cone([[1,1,0]], QQ, HORIZONTAL)
                Traceback (most recent call last):
                ...
                ValueError: the constraints on the residues force some coordinates to be zero in the Strebel polyhedron of StrebelGraph("(0:0,~1:0)(1:0)(~0:0)")
            
            Here are valid examples::
                sage: G0 = StrebelGraph("(0)(~0)")
                sage: G0.residue_constraint_cone([[1,1]], QQ, HORIZONTAL)
                Cone of dimension 2 in ambient dimension 2 made of 1 facets (backend=ppl)
                sage: G0.residue_constraint_cone([[1,1]], QQ, HORIZONTAL).rays()
                [[0, 1]]
                sage: G1 = StrebelGraph("(0,~1)(1)(~0)")
                sage: G1.residue_constraint_cone([[0,1,-1]],QQ,VERTICAL)
                Cone of dimension 2 in ambient dimension 4 made of 1 facets (backend=ppl)
                sage: G1.residue_constraint_cone([[0,1,1]],QQ,VERTICAL).rays()
                [[1, 1, 0, 0]]
                sage: G2 = StrebelGraph("(0,2,~3,~1)(1)(3,~0)(~2)")
                sage: G2.residue_constraint_cone([[1,2,0,0]], QQ, VERTICAL)
                Cone of dimension 6 in ambient dimension 8 made of 3 facets (backend=ppl)
                sage: G2.residue_constraint_cone([[1,2,0,0]], QQ, VERTICAL).rays()
                [[1, 1, 0, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 0, 0, 0]]
                sage: G2.residue_constraint_cone([[0,1,-1,0],[0,1,0,-1]], QQ, VERTICAL)
                Cone of dimension 4 in ambient dimension 8 made of 2 facets (backend=ppl)
                sage: G2.residue_constraint_cone([[0,1,-1,0],[0,1,0,-1]], QQ, VERTICAL).rays()
                [[1, 1, 1, 0, 0, 0, 0, 0], [0, 1, 1, 1, 0, 0, 0, 0]]
            """
            
            L = LinearExpressions(K)
            nhe = self._n
            cs = ConstraintSystem(dim = nhe)
            x = [L.variable(i) for i in range(nhe)]

            self._set_strebel_constraints(cs, L, slope)
            
            self._set_residue_constraints(cs, coes, L, VERTICAL)
            self._set_residue_constraints(cs, coes, L, HORIZONTAL) 
            
            #check if the constraints force some x[i] to be zero
            d = cs.cone().affine_dimension()
            for i in range(nhe):
                cs1 = cs.copy()
                cs1.insert(x[i] == 0)
                if cs1.cone().affine_dimension() == d:
                    raise ValueError('the constraints on the residues force some coordinates to be zero in the Strebel polyhedron of {}'.format(self))
            
            return cs.cone()