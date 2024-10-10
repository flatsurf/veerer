import numbers
from array import array

from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE, rich_to_bool

from .permutation import (perm_init, perm_check, perm_cycles, perm_dense_cycles,
                          perm_invert, perm_conjugate, perm_cycle_string, perm_cycles_lengths,
                          perm_cycles_to_string, perm_on_list,
                          perm_num_cycles, str_to_cycles, str_to_cycles_and_data, perm_compose, perm_from_base64_str,
                          uint_base64_str, uint_from_base64_str, perm_base64_str,
                          perms_are_transitive, triangulation_relabelling_from)

from surface_dynamics import Stratum
from sage.rings.all import ZZ


def face_edge_perms_init(faces):
    r"""
    INPUT:: ``faces`` - a list or a string encoding a permutation

    EXAMPLES::

        sage: from veerer.triangulation import face_edge_perms_init  # random output due to deprecation warnings from realalg

        sage: face_edge_perms_init('(0,1,2)(~0,~1,~2)')
        (array('i', [1, 2, 0, 5, 3, 4]), array('i', [5, 4, 3, 2, 1, 0]))

        sage: face_edge_perms_init('(0,1,2)')
        (array('i', [1, 2, 0]), array('i', [0, 1, 2]))

        sage: face_edge_perms_init('(0,1,2)(~0)(~1)(~2)')
        (array('i', [1, 2, 0, 3, 4, 5]), array('i', [5, 4, 3, 2, 1, 0]))

    TESTS:

    Check that edge permutation do not depend on the details of faces::

        sage: f1 = "(0,~5,4)(3,5,6)(1,2,~6)"
        sage: f2 = "(0,6,5)(1,2,~6)(3,4,~5)"
        sage: f3 = "(6,4,3)(~6,~5,1)(5,0,2)"
        sage: assert face_edge_perms_init(f1)[1] == face_edge_perms_init(f2)[1] == face_edge_perms_init(f3)[1]
    """
    if isinstance(faces, str):
        l = str_to_cycles(faces)
    else:
        l = [[int(i) for i in c] for c in faces]

    pos = []
    neg = []
    for c in l:
        for e in c:
            if e < 0:
                neg.append(e)
            else:
                pos.append(e)

    for e in neg:
        if ~e not in pos:
            raise ValueError("inconsistent permutation data")

    pos.sort()
    neg.sort(reverse=True)
    for i in range(len(pos) - 1):
        if pos[i] == pos[i+1]:
            raise ValueError("repeated edge label {}".format(pos[i]))
        elif pos[i+1] != pos[i] + 1:
            raise ValueError("missing edge label {}".format(pos[i] + 1))
    for i in range(len(neg) - 1):
        if neg[i] == neg[i+1]:
            raise ValueError("repeated edge label ~{}".format(~neg[i]))

    # number of half edges
    n = len(pos) + len(neg)

    # build the edge permutation ep
    ep = [-1] * n  # edge permutation
    m = n - 1   # label available at the back
    for e in neg:
        E = ~e
        e = m
        if ep[E] != -1:
            raise ValueError("inconsistent permutation data")
        m -= 1
        ep[E] = e
        ep[e] = E
    for i in range(n):
        if ep[i] == -1:
            ep[i] = i

    fp = [-1] * n  # face permutation
    for c in l:
        k = len(c)
        for i in range(k):
            e0 = c[i]
            e1 = c[(i + 1) % k]
            if e0 < 0:
                e0 = ep[~e0]
            if e1 < 0:
                e1 = ep[~e1]
            fp[e0] = e1

    return perm_init(fp), perm_init(ep)


def boundary_init(fp, ep, boundary):
    n = len(ep)

    if boundary is None:
        return array('i', [0] * n)
    elif isinstance(boundary, (array, tuple, list)):
        if len(boundary) != n:
            raise ValueError('invalid input argument')
        return array('i', boundary)
    elif isinstance(boundary, dict):
        output = array('i', [0] * n)
        for e, v in boundary.items():
            if isinstance(e, str):
                if not e:
                    raise ValueError('keys must be valid edges, got {!r}'.format(e))
                elif e[0] == '~':
                    e = ~int(e[1:])
                else:
                    e = int(e)
            elif isinstance(e, numbers.Integral):
                e = int(e)
            if e > n:
                raise ValueError('keys must be valid edges, got {!r}'.format(e))
            output[e] = v
        return output
    else:
        raise TypeError('invalid boundary data')

        
class StrebelGraph(object):
    
    __slots__ = ['_mutable', '_n', '_fp', '_ep', '_vp', '_bdry']
        
    def __init__(self, triangles, boundary=None, mutable=False, check=True):
        if isinstance(triangles, StrebelGraph):
            self._fp = triangles.face_permutation(copy=True)
            self._ep = triangles.edge_permutation(copy=True)
            self._bdry = triangles.boundary_vector(copy=True)
        else:
            if boundary is not None and isinstance(boundary, str):
                boundary_cycles, boundary = str_to_cycles_and_data(boundary)
                if isinstance(triangles, str):
                    triangles = str_to_cycles(triangles)
                else:
                    triangles = list(triangles)
                triangles.extend(boundary_cycles)
            self._fp, self._ep = face_edge_perms_init(triangles)
            self._bdry = boundary_init(self._fp, self._ep, boundary)

        self._mutable = mutable

        fp = self._fp
        ep = self._ep
        n = self._n = len(fp)

        vp = self._vp = array('i', [-1] * n)
        for i in range(n):
            vp[fp[ep[i]]] = i

        if check:
            self._check(ValueError)
    
        
    def _check(self, error=RuntimeError):
        
        n = self._n

        if not (hasattr(self, '_vp') and hasattr(self, '_ep') and hasattr(self, '_fp') and hasattr(self, '_bdry')):
            raise error('missing attributes: these must be _vp, _ep, _fp, _bdry')
        if not perm_check(self._fp, n):
            raise error('fp is not a permutation: {}'.format(self._fp))
        if not perm_check(self._ep, n):
            raise error('ep is not permutation: {}'.format(self._ep))
        if not perm_check(self._vp, n):
            raise error('vp is not a permutation: {}'.format(self._vp))
        if not perms_are_transitive([self._fp, self._ep, self._vp]):
            raise error('(fp, ep, vp) do not generate a transitive group')
        if not isinstance(self._bdry, array) or self._bdry.typecode != 'i' or len(self._bdry) != n:
            raise error('bdry must be an integer array of same length as the underlying permutations')

        for i in range(n):
            if self._ep[self._ep[i]] != i:
                raise error('invalid edge permutation at half-edge i={}'.format(self._half_edge_string(i)))
            if self._fp[self._ep[self._vp[i]]] != i:
                raise error('fev relation not satisfied at half-edge i={}'.format(self._half_edge_string(i)))
    
    def _check_half_edge(self, e):
        if not isinstance(e, numbers.Integral):
            raise TypeError('invalid half-edge {}'.format(e))
        e = int(e)
        if e < 0 or e >= self._n:
            raise ValueError('half-edge number out of range e={}'.format(e))
        return e

    def set_immutable(self):
        self._mutable = False

    def __hash__(self):
        if self._mutable:
            raise ValueError('mutable veering triangulation not hashable')

        x = 140737488617563
        x = ((x ^ hash(self._vp.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._ep.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._fp.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._bdry.tobytes())) * 2147483693) + 82520 + self._n + self._n

        return x
    
    def copy(self, mutable=None):
        
        if mutable is None:
            mutable = self._mutable

        if not self._mutable and not mutable:
            # avoid copies of immutable objects
            if type(self) is StrebelGraph:
                return self
            else:
                T = StrebelGraph.__new__(StrebelGraph)
                T._n = self._n
                T._fp = self._fp
                T._ep = self._ep
                T._vp = self._vp
                T._bdry = self._bdry
                T._mutable = mutable

                return T
        else:
            T = StrebelGraph.__new__(StrebelGraph)
            T._n = self._n
            T._fp = self._fp[:]
            T._ep = self._ep[:]
            T._vp = self._vp[:]
            T._bdry = self._bdry[:]
            T._mutable = mutable

            return T
    
        
    def flip(self, e, check=True):
        r"""
        Flip a half-edge of the Strebel graph.
        
        All half-plane excesses are required to be zero within this method. 
        
        EXAMPLES::
        
            sage: from strebelgraph import *
            
            sage: G = Strebelgraph("(0,1,2,3)(~0,4,5,6)(~1,~6,~5,~4,~3,~2)",mutable = True)
            sage: G.flip(0)
            sage: G1 = Strebelgraph("(0,2,3,4)(~0,5,6,1)(~1,~6,~5,~4,~3,~2)")
            sage: G._richcmp_(G1,1)
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
        
            sage: from strebelgraph import *
        
        Flip a half-edge with zero half-plane excess::  
        
            sage: G = Strebelgraph("(0,1,2)", boundary="(~0:0)(~1:1,~2:1)", mutable = True)
            sage: G.halfdualflip(0)
            sage: G1 = Strebelgraph("(0,~2,~1)(1,2)(~0)", boundary = {"~1":1, "~2":1})
            sage: G._richcmp_(G1, 1)
            True
            
        Flip a half-edge with non-zero half-plane excess::
        
            sage: G = Strebelgraph("(0,1,2)", boundary="(~0:0)(~1:1,~2:1)", mutable = True)
            sage: G.halfdualflip(4)
            sage: G1 = Strebelgraph("(0,1,2)(~1,~2)(~0)", boundary = {"~1":0, "~2":2})
            sage: G._richcmp_(G1, 1)
            True
            
        The graph remains the same if ``e`` is in a loop edge::
        
            sage: G = Strebelgraph("(0,1,2)", boundary="(~0:0)(~1:1,~2:1)", mutable = True)
            sage: G1 = G.copy()
            sage: G.halfdualflip(5)
            sage: G._richcmp_(G1, 1)
            True
            
        The graph remains the same if ``e`` is an end of the graph with half-plane excess::
        
            sage: G = Strebelgraph("(0,~1,1,~2,2,~0)", mutable = True)
            sage: G1 = G.copy()
            sage: G.halfdualflip(0)
            sage: G._richcmp_(G1, 1)
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
        
            sage: from strebelgraph import *
            
            sage: G = Strebelgraph("(0,1,2)", boundary="(~0:0)(~1:1,~2:1)", mutable = True)
            sage: G.dualflip(1)
            sage: G1 = Strebelgraph("(0,2)(~0,1)(~1,~2)", boundary = {"~1": 0, "~2":2})
            sage: G._richcmp_(G1, 1)
            True
        We obtain the same Strebel graphs no matter we input ``e`` 
        or the edge permutation of ``e``::
            sage: G = Strebelgraph("(0,1,2)", boundary="(~0:0)(~1:1,~2:1)", mutable = True)
            sage: G.dualflip(4)
            sage: G._richcmp_(S1,1)
            True
            
        Flip the edge where ``e`` is an end of the Strebel graph::
        
            sage: G = Strebelgraph("(0,~1,1,~2,2,~0)", mutable = True)
            sage: G.dualflip(0)
            sage: G1 = Strebelgraph("(0,~2,2,~1,1,~0)")
            sage: G._richcmp_(S1,1)
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
            
            sage: from strebelgraph import *
        
            sage: G = Strebelgraph("(0,1,2)", boundary="(~0:0)(~1:1,~2:1)", mutable = True)
            sage: S = [2,1]
            sage: G.fp_ordered(S)
            [[4], [3], [1, 2]]
         
        The following is an example when ``S`` contains an end of the graph::
        
            sage: G = Strebelgraph("(0,~1,1,~2,2,~0)", mutable = True)
            sage: S = [5]
            sage: G.fp_ordered(S)
            [[5,0]]
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
            sage: from strebelgraph import *
            
            sage: G = Strebelgraph("(0,1,2)", boundary="(~0:0)(~1:1,~2:1)", mutable = True)
            sage: S = [1,2]
            sage: G.multidualflip(S)
            sage: G1 = Strebelgraph("(~0,1,2)(0)(~1,~2)", boundary = {"~1":1, "~2":1})
            sage: G._richcmp_(G1, 1)
            Ture
            
        The following is an example when ``S`` contains an end of the graph::
            
            sage: G = Strebelgraph("(0,~1,1,~2,2,~0)", mutable = True)
            sage: S = [5]
            sage: G.multidualflip(S)
            sage: G1 = Strebelgraph("(0,~2,2,~1,1,~0)")
            sage: G._richcmp_(S1,1)
            True
        """
        
        ll = self.fp_ordered(S)
        
        for l in ll:
            for e in l:
                self.relhalfdualflip(e, S)              
                
    def faces(self):
        r"""
        Return the list of faces as tuples of half-edges
        """
        
        return [c for c in perm_cycles(self._fp, True, self._n)]
    
    def vertices(self):
        r"""
        Return the list of vertices as tuples of half-edges
        """
        return perm_cycles(self._vp, True, self._n) 
    
    def _richcmp_(self, other, op):
        
        c = (self._n > other._n) - (self._n < other._n)
        if c:
            return rich_to_bool(op, c)

        c = (self._bdry > other._bdry) - (self._bdry < other._bdry)
        if c:
            return rich_to_bool(op, c)

        c = (self._fp > other._fp) - (self._fp < other._fp)
        if c:
            return rich_to_bool(op, c)
        
        c = (self._vp > other._vp) - (self._vp < other._vp)
        if c:
            return rich_to_bool(op, c)

        c = (self._ep > other._ep) - (self._ep < other._ep)
        return rich_to_bool(op, c)

        
        
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
            
   
    def num_vertices(self):
        r"""
        Return the number of vertices.
        """
        return perm_num_cycles(self._vp)
                
    
    def strebel_in_stratum(self):
        
        k = self.is_abelian(certificate=False)
        mu = []
        
        lv = self.vertices() #list of the vertives
        lf = self.faces() #list of faces
        
        if k%2: #abelian
            for l in lv:
                dv = len(l) #the degree of the vertex in G
                s = sum(self._bdry[i] for i in l)
                muv = ZZ((1/2) * s + (1/2) * dv - 1)
                mu.append(muv)
            for f in lf:
                s = sum(self._bdry[i] for i in f)
                muv = ZZ(-(1/2) * s - 1)
                mu.append(muv)
            k = 1
        else:
            for l in lv:
                dv = len(l) #the degree of the vertex in G
                muv = ZZ(sum(self._bdry[i] for i in l) + dv - 2)
                mu.append(muv)
            for f in lf:
                muv = ZZ(-sum(self._bdry[i] for i in f) - 2)
                mu.append(muv)
            k = 2
        
        return Stratum(mu,k)
        
        
    def _half_edge_string(self, e):
        f = self._ep[e]
        return '~%d' % f if f < e else '%d' % e
    
    def strebel_from_face_edge_perms(vp, ep, fp=None, boundary=None, mutable=False, check=True):
        r"""
        INPUT:

        - ``fp``, ``ep``, ``vp`` -- the face, edge and vertex permutation

        - ``check`` - boolean (default: ``True``) - if set to ``False`` no
          check are performed

        EXAMPLES::

            sage: from strebelgraph import *
            sage: from array import array
            sage: vp = array('i', [1, 3, 0, 2])
            sage: ep = array('i', [3, 2, 1, 0])
            sage: Strebelgraph.strebel_from_face_edge_perms(vp, ep, boundary = array('i', [1, 0, 0, 1]))
            
            ?????
            sage: fp = array('i', [1, 2, 0, 4, 8, 6, 7, 5, 3])
            sage: ep = array('i', [8, 7, 2, 3, 4, 5, 6, 1, 0])
            sage: vp = array('i', [2, 8, 7, 0, 3, 1, 5, 6, 4])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            Triangulation("(0,1,2)(3,4,~0)(5,6,~1)")
        """
        G = StrebelGraph.__new__(StrebelGraph)
        n = G._n = len(vp)
        
        G._vp = vp
        G._ep = ep
        
        if fp is None:
            vp = G._vp
            ep = G._ep
            fp = array('i', [-1] * n)
            for i in range(n):
                fp[ep[vp[i]]] = i
        G._fp = fp
        if boundary is None:
            G._bdry = array('i', [0] * n)
        else:
            G._bdry = array('i', boundary)

        G._mutable = mutable

        if check:
            G._check(ValueError)

        return G
    
    def relabel(self, p, check=True):

        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        n = self._n
        if check and not perm_check(p, n):
            # if the input is not a valid permutation, we assume that half-edges
            # are not separated
            p = perm_init(p, n, self._ep)
            if not perm_check(p, n, self._ep):
                raise ValueError('invalid relabeling permutation')

        self._vp = perm_conjugate(self._vp, p)
        self._ep = perm_conjugate(self._ep, p)
        self._fp = perm_conjugate(self._fp, p)
        perm_on_list(p, self._bdry, n)

    def _relabelling_from(self, start_edge):

        if start_edge < 0 or start_edge >= len(self._vp):
            raise ValueError
        return triangulation_relabelling_from(self._vp, self._ep, start_edge)
    
    def _automorphism_good_starts(self):

        return range(self._n)

    def best_relabelling(self, all=False):
        n = self._n
        fp = self._fp
        ep = self._ep

        best = None
        if all:
            relabellings = []

        for start_edge in self._automorphism_good_starts():
            relabelling = self._relabelling_from(start_edge)

            fp_new = perm_conjugate(fp, relabelling)
            ep_new = perm_conjugate(ep, relabelling)
            bdry_new = self._bdry[:]
            perm_on_list(relabelling, bdry_new, self._n)

            T = (fp_new, ep_new, bdry_new)
            if best is None or T < best:
                best_relabelling = relabelling
                best = T
                if all:
                    del relabelllings[:]
                    relabellings.append(relabelling)
            elif all and T == best:
                relabellings.append(relabelling)

        return (relabellings, best) if all else (best_relabelling, best)

    def set_canonical_labels(self):

        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        r, _ = self.best_relabelling()
        self.relabel(r, check=False)
       
        
        
        
        
        
        
        
        