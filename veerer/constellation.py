r"""
Constellations (possibly with boundary data)

Common base claas for class:~veerer.triangulation.Triangulations and :class:~veerer.strebel_graph.StrebelGraph
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2018 Mark Bell
#                     2018-2024 Vincent Delecroix
#                     2024 Kai Fu
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

import numbers
from array import array

from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE, rich_to_bool

from .permutation import (perm_init, perm_check, perm_cycles, perm_dense_cycles,
                          perm_invert, perm_conjugate, perm_cycle_string, perm_cycles_lengths,
                          perm_cycles_to_string, perm_on_list,
                          perm_num_cycles, str_to_cycles, str_to_cycles_and_data, perm_compose, perm_from_base64_str,
                          uint_base64_str, uint_from_base64_str, perm_base64_str,
                          perms_are_transitive, triangulation_relabelling_from)

class Constellation:
    __slots__ = ['_mutable', '_n', '_fp', '_ep', '_vp', '_bdry']

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

    def __getstate__(self):
        r"""
        TESTS::

            sage: from veerer import Triangulation
            sage: t = Triangulation("(0,1,2)")
            sage: dumps(t)  # indirect doctest
            b'x\x9ck`J.KM-J-\xd2+)\xcaL\xccK/\xcdI,\xc9\xcc\xcf\xe3\nA\xe1\x152h6\x162\xc6\x162ix3z3y3\x00!\x8cf\xe8LM\xd2\x03\x00_u\x15?'
        """
        a = list(self._fp)
        a.extend(self._ep)
        a.extend(self._bdry)
        a.append(self._mutable)
        return a

    def __setstate__(self, arg):
        r"""
            sage: from veerer import Triangulation
            sage: t0 = Triangulation("(0,1,2)", mutable=False)
            sage: t1 = Triangulation("(0,1,2)", mutable=True)
            sage: s0 = loads(dumps(t0))  # indirect doctest
            sage: assert s0 == t0
            sage: s0._mutable
            False
            sage: s0._check()

            sage: s1 = loads(dumps(t1))  # indirect doctest
            sage: assert s1 == t1
            sage: s1._mutable
            True
            sage: s1._check()
        """
        n = (len(arg) - 1) // 3
        self._n = n
        self._fp = array('i', arg[:n])
        self._ep = array('i', arg[n:2*n])
        self._bdry = array('i', arg[2*n:3*n])
        self._mutable = arg[-1]
        self._vp = array('i', [-1] * n)
        for i in range(n):
            self._vp[self._fp[self._ep[i]]] = i

    def set_immutable(self):
        self._mutable = False

    def __hash__(self):
        r"""
        TESTS::

            sage: from itertools import permutations, combinations
            sage: from veerer import Triangulation

            sage: triangulations = []

            sage: t = Triangulation("(0, 1, 2)")
            sage: triangulations.append(t)

            sage: for p in permutations(["1", "~1", "2", "~2"]):
            ....:     t = Triangulation("(0, {}, {})(~0, {}, {})".format(*p))
            ....:     triangulations.append(t)

            sage: for i, j in combinations([0, 1, 2, 3], 2):
            ....:     for k, l in permutations(set(range(4)).difference([i, j])):
            ....:         vars = {'i': i, 'j': j, 'k': k, 'l': l}
            ....:         t = Triangulation("({i}, {j}, {k})(~{i}, ~{j}, {l})".format(**vars))
            ....:         triangulations.append(t)
            ....:         t = Triangulation("({i}, ~{j}, {k})(~{i}, {j}, {l})".format(**vars))
            ....:         triangulations.append(t)
            ....:         t = Triangulation("({i}, {k}, {j})(~{i}, ~{j}, {l})".format(**vars))
            ....:         triangulations.append(t)
            ....:         t = Triangulation("({i}, {k}, ~{j})(~{i}, {j}, {l})".format(**vars))
            ....:         triangulations.append(t)
            ....:         t = Triangulation("({i}, {j}, {k})(~{i}, {l}, ~{j})".format(**vars))
            ....:         triangulations.append(t)
            ....:         t = Triangulation("({i}, ~{j}, {k})(~{i}, {l}, {j})".format(**vars))
            ....:         triangulations.append(t)
            ....:         t = Triangulation("({i}, {k}, {j})(~{i}, {l}, ~{j})".format(**vars))
            ....:         triangulations.append(t)
            ....:         t = Triangulation("({i}, {k}, ~{j})(~{i}, {l}, {j})".format(**vars))
            ....:         triangulations.append(t)

            sage: for i in range(len(triangulations)):
            ....:     for j in range(len(triangulations)):
            ....:         assert (triangulations[i] == triangulations[j]) == (i == j), (i, j)
            ....:         assert (triangulations[i] != triangulations[j]) == (i != j), (i, j)

            sage: hashes = {}
            sage: for t in triangulations:
            ....:     h = hash(t)
            ....:     if h in hashes:
            ....:         print('collision: {} {}'.format(hashes[h], t))
            ....:     else:
            ....:         hashes[h] = t
            sage: assert len(hashes) == len(triangulations), (len(hashes), len(triangulations))
        """
        if self._mutable:
            raise ValueError('mutable veering triangulation not hashable')

        x = 140737488617563
        x = ((x ^ hash(self._vp.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._ep.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._fp.tobytes())) * 2147483693) + 82520 + self._n + self._n

        return x

    def _check_half_edge(self, e):
        if not isinstance(e, numbers.Integral):
            raise TypeError('invalid half-edge {}'.format(e))
        e = int(e)
        if e < 0 or e >= self._n:
            raise ValueError('half-edge number out of range e={}'.format(e))
        return e

    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n == other._n and self._fp == other._fp and self._ep == other._ep and self._bdry == other._bdry

    def __ne__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n != other._n or self._fp != other._fp or self._ep != other._ep or self._bdry != other._bdry

    def _richcmp_(self, other, op):
        c = (self._n > other._n) - (self._n < other._n)
        if c:
            return rich_to_bool(op, c)

        c = (self._fp > other._fp) - (self._fp < other._fp)
        if c:
            return rich_to_bool(op, c)

        c = (self._ep > other._ep) - (self._ep < other._ep)
        if c:
            return rich_to_bool(op, c)

        c = (self._bdry > other._bdry) - (self._bdry < other._bdry)
        return rich_to_bool(op, c)

    def __lt__(self, other):
        return self._richcmp_(other, op_LT)

    def __le__(self, other):
        return self._richcmp_(other, op_LE)

    def __gt__(self, other):
        return self._richcmp_(other, op_GT)

    def __ge__(self, other):
        return self._richcmp_(other, op_GE)

    def copy(self, mutable=None, cls=None):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation([[0,1,2],[-1,-2,-3]], mutable=True)
            sage: U = T.copy()
            sage: T == U
            True
            sage: T.flip(0)
            sage: T == U
            False

            sage: T.set_immutable()
            sage: T.copy() is T
            True

            sage: U = T.copy(mutable=True)
            sage: U.flip(0)
            sage: T
            Triangulation("(0,2,~1)(1,~0,~2)")

        TESTS::

            sage: from veerer import Triangulation
            sage: T = Triangulation("(0,1,2)(~0,~1,~2)", mutable=True)
            sage: U = T.copy(mutable=False)
            sage: _ = hash(U)
        """
        if mutable is None:
            mutable = self._mutable
        if cls is None:
            cls = self.__class__

        if not self._mutable and not mutable:
            # avoid copies of immutable objects
            if type(self) is cls:
                return self
            else:
                T = cls.__new__(cls)
                T._n = self._n
                T._fp = self._fp
                T._ep = self._ep
                T._vp = self._vp
                T._bdry = self._bdry
                T._mutable = mutable

                return T
        else:
            T = cls.__new__(cls)
            T._n = self._n
            T._fp = self._fp[:]
            T._ep = self._ep[:]
            T._vp = self._vp[:]
            T._bdry = self._bdry[:]
            T._mutable = mutable

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

    def boundary_vector(self, copy=True):
        if copy:
            return self._bdry[:]
        else:
            return self._bdry

    def next_at_vertex(self, i, check=True):
        r"""
        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: T.next_at_vertex(0)
            9
            sage: T.next_at_vertex(9)
            8
            sage: T.next_at_vertex(5)
            4
        """
        if check:
            i = self._check_half_edge(i)
        return self._vp[i]

    def next_in_face(self, i, check=True):
        r"""
        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: T.next_in_face(0)
            10
            sage: T.next_in_face(9)
            1
            sage: T.next_in_face(5)
            3
        """
        if check:
            i = self._check_half_edge(i)
        return self._fp[i]

    def previous_at_vertex(self, i, check=True):
        r"""
        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: T.previous_at_vertex(9)
            0
            sage: T.previous_at_vertex(8)
            9
            sage: T.previous_at_vertex(4)
            5
        """
        if check:
            i = self._check_half_edge(i)
        return self._fp[self._ep[i]]

    def previous_in_face(self, i, check=True):
        r"""
        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: T.previous_in_face(10)
            0
            sage: T.previous_in_face(1)
            9
            sage: T.previous_in_face(3)
            5
        """
        if check:
            i = self._check_half_edge(i)
        return self._ep[self._vp[i]]

    def num_half_edges(self):
        r"""
        Return the number of half edges.
        """
        return self._n

    def folded_edges(self):
        r"""
        Return the list of darts that belong to folded edges.
        """
        n = self._n
        ep = self._ep
        return [i for i in range(n) if ep[i] == i]

    def num_folded_edges(self):
        r"""
        Return the number of folded edges.
        """
        n = self._n
        ep = self._ep
        return sum(ep[i] == i for i in range(n))

    def num_edges(self):
        r"""
        Return the number of edges.
        """
        return (self._n + self.num_folded_edges()) // 2

    def _edge_rep(self, e):
        f = self._ep[e]
        if f < e:
            return '~%d' % f
        else:
            return str(e)

    def _norm(self, e):
        f = self._ep[e]
        return f if f < e else e

    def _half_edge_string(self, e):
        f = self._ep[e]
        return '~%d' % f if f < e else '%d' % e

    def edges(self):
        r"""
        Return the list of edges as tuples of half-edges

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.edges()
            [[0, 8], [1], [2], [3, 7], [4], [5], [6]]
        """
        return perm_cycles(self._ep, True, self._n)

    def vertices(self):
        r"""
        Return the list of vertices as tuples of half-edges

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.vertices()
            [[0, 2, 1, 8, 6, 3, 5, 4, 7]]
        """
        return perm_cycles(self._vp, True, self._n)

    def num_vertices(self):
        r"""
        Return the number of vertices.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.num_vertices()
        """
        return perm_num_cycles(self._vp, self._n)

    def faces(self):
        r"""
        Return the list of edges as tuples of half-edges

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.faces()
            [[0, 1, 2], [3, 4, 5], [6, 8, 7]]
        """
        return perm_cycles(self._fp, True, self._n)

    def num_faces(self):
        r"""
        Return the number of faces.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.num_faces()
            3
        """
        return perm_num_cycles(self._fp, self._n)

    def swap(self, e, check=True):
        r"""
        Change the orientation of the edge ``e``.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)", mutable=True)
            sage: T.swap(0)
            sage: T
            Triangulation("(0,~1,~2)(1,2,~0)")
            sage: T.swap(1)
            sage: T
            Triangulation("(0,1,~2)(2,~0,~1)")
            sage: T.swap(2)
            sage: T
            Triangulation("(0,1,2)(~2,~0,~1)")

            sage: T = Triangulation("(0,~5,4)(3,5,6)(1,2,~6)", mutable=True)
            sage: T.swap(0)
            sage: T
            Triangulation("(0,~5,4)(1,2,~6)(3,5,6)")
            sage: T.swap(5)
            sage: T
            Triangulation("(0,5,4)(1,2,~6)(3,~5,6)")

        Also works for veering triangulations::

            sage: from veerer import VeeringTriangulation

            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: V = VeeringTriangulation(fp, cols, mutable=True)
            sage: V.swap(0)
            sage: V.swap(10)
            sage: V
            VeeringTriangulation("(0,1,~3)(2,~0,~1)(3,4,~5)(5,~9,~7)(6,~2,~4)(7,~6,8)(9,~10,~11)(10,11,~8)", "BRBBBRRBBBBR")

        One can alternatively use ``relabel``::

            sage: T = Triangulation("(0,~5,4)(3,5,6)(1,2,~6)", mutable=True)
            sage: T1 = T.copy()
            sage: T1.swap(3)
            sage: T1.swap(6)
            sage: T2 = T.copy()
            sage: T2.relabel("(3,~3)(6,~6)")
            sage: T1 == T2
            True
        """
        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        if check:
            e = self._check_half_edge(e)

        vp = self._vp
        ep = self._ep
        fp = self._fp
        E = ep[e]

        if e == E:
            return

        # images/preimages by vp
        e_vp = vp[e]
        E_vp = vp[E]
        e_vp_inv = fp[E]
        E_vp_inv = fp[e]
        assert vp[e_vp_inv] == e
        assert vp[E_vp_inv] == E

        # images/preimages by fp
        e_fp = fp[e]
        E_fp = fp[E]
        e_fp_inv = ep[e_vp]
        E_fp_inv = ep[E_vp]
        assert fp[e_fp_inv] == e
        assert fp[E_fp_inv] == E

        fp[e_fp_inv] = E
        fp[E_fp_inv] = e
        vp[e_vp_inv] = E
        vp[E_vp_inv] = e
        fp[e] = E_fp
        fp[E] = e_fp
        vp[e] = E_vp
        vp[E] = e_vp

        self._bdry[e], self._bdry[E] = self._bdry[E], self._bdry[e]

    def relabel(self, p, check=True):
        r"""
        Relabel this triangulation inplace according to the permutation ``p``.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)", mutable=True)
            sage: T.relabel("(0,~0)")
            sage: T
            Triangulation("(0,~1,~2)(1,2,~0)")
            sage: T.relabel("(0,1,~2)")
            sage: T
            Triangulation("(0,1,2)(~2,~0,~1)")

            sage: T.set_immutable()
            sage: T.relabel("(0,~1)")
            Traceback (most recent call last):
            ...
            ValueError: immutable triangulation; use a mutable copy instead

        An example of a flip sequence which forms a loop after non-trivial relabelling::

            sage: T0 = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)", mutable=True)
            sage: T = T0.copy()
            sage: T.flip_back(1)
            sage: T.flip_back(3)
            sage: T.flip_back(0)
            sage: T.flip_back(2)
            sage: T.relabel("(0,2)(1,3)")
            sage: T == T0
            True

        An example with boundary::

            sage: t = Triangulation("(0,1,2)(~0,3,4)(~4,~3,~2,~1)", {"~1": 1, "~2": 1, "~3": 1, "~4": 1}, mutable=True)
            sage: t.relabel("(0,3)(1,~2)")
            sage: t
            Triangulation("(0,4,~3)(3,~2,~1)", boundary="(1:1,2:1,~4:1,~0:1)")
        """
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
        r"""
        Return a canonical relabelling map obtained from walking
        along the triangulation starting at ``start_edge``.

        The returned relabelling array maps the current edge to the new
        labelling.

        EXAMPLES::

            sage: from veerer import *
            sage: from array import array

        The torus example (6 symmetries)::

            sage: fp = array('i', [1, 5, 4, 2, 3, 0])
            sage: ep = array('i', [4, 3, 5, 1, 0, 2])
            sage: vp = array('i', [2, 4, 1, 0, 5, 3])
            sage: T = Triangulation.from_face_edge_perms(fp, ep, vp, mutable=True)
            sage: T._relabelling_from(3)
            array('i', [4, 5, 3, 0, 1, 2])

            sage: p = T._relabelling_from(0)
            sage: T.relabel(p)
            sage: for i in range(6):
            ....:     p = T._relabelling_from(i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     assert S == T

        The sphere example (3 symmetries)::

            sage: fp = array('i', [1, 2, 0])
            sage: ep = array('i', [0, 1, 2])
            sage: vp = array('i', [2, 0, 1])
            sage: T = Triangulation.from_face_edge_perms(fp, ep, vp, mutable=True)
            sage: T._relabelling_from(1)
            array('i', [1, 0, 2])
            sage: p = T._relabelling_from(0)
            sage: T.relabel(p)
            sage: for i in range(3):
            ....:     p = T._relabelling_from(i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     assert S == T

        An example with no automorphism::

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)", mutable=True)
            sage: p = T._relabelling_from(0)
            sage: T.relabel(p)
            sage: for i in range(1, 9):
            ....:     p = T._relabelling_from(i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     S._check()
            ....:     assert S != T
        """
        if start_edge < 0 or start_edge >= len(self._vp):
            raise ValueError
        return triangulation_relabelling_from(self._vp, self._ep, start_edge)

    def _automorphism_good_starts(self):
        # TODO: we should also discriminate based on boundaries!
        # we discriminate based on the lengths of cycles in
        #   vp, ep, fp and vp[ep[fp]]
        n = self._n

        lv = perm_cycles_lengths(self._vp, n)
        le = perm_cycles_lengths(self._ep, n)
        lf = perm_cycles_lengths(self._fp, n)
        vef = perm_compose(perm_compose(self._fp, self._ep, n), self._vp, n)
        lvef = perm_cycles_lengths(vef, n)

        d = {}
        for i in range(n):
            s = (((self._bdry[i] * n + lv[i]) * n + le[i]) * n + lf[i]) * n + lvef[i]
            if s in d:
                d[s].append(i)
            else:
                d[s] = [i]

        winner = None
        values = None
        for s, v in d.items():
            if winner is None:
                winner = s
                values = v
            elif len(v) < len(values):
                winner = s
                values = v
            elif len(v) == len(values) and s < winner:
                winner = s

        return d[winner]

    def automorphisms(self):
        r"""
        Return the list of automorphisms of this triangulation.

        The output is a list of arrays that are permutations acting on the set
        of half edges.

        For triangulations with boundaries, we allow automorphism to permute
        boundaries. Though, boundary edge have to be mapped on boundary edge.

        EXAMPLES::

            sage: from veerer import *

        An example with 4 symmetries in genus 2::

            sage: T = Triangulation("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)")
            sage: A = T.automorphisms()
            sage: len(A)
            4

        And the "sphere octagon" has 8::

            sage: s  = "(0,8,~7)(1,9,~0)(2,10,~1)(3,11,~2)(4,12,~3)(5,13,~4)(6,14,~5)(7,15,~6)"
            sage: len(Triangulation(s).automorphisms())
            8

        A veering triangulation with 4 symmetries in genus 2::

            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: V = VeeringTriangulation(fp, cols)
            sage: A = V.automorphisms()
            sage: len(A)
            4
            sage: S = V.copy(mutable=True)
            sage: for a in A:
            ....:     S.relabel(a)
            ....:     assert S == V

        Examples with boundaries::

            sage: t = Triangulation("(0,1,2)", boundary="(~0:1)(~1:1)(~2:1)")
            sage: len(t.automorphisms())
            3
            sage: t = Triangulation("(0,1,2)", boundary="(~0:1,~1:1,~2:1)")
            sage: len(t.automorphisms())
            3
            sage: t = Triangulation("(0,1,2)", boundary="(~0:1,~1:1,~2:2)")
            sage: len(t.automorphisms())
            1
        """
        fp = self._fp
        ep = self._ep
        bdry = self._bdry

        best = None
        best_relabellings = []

        for start_edge in self._automorphism_good_starts():
            relabelling = self._relabelling_from(start_edge)

            fp_new = perm_conjugate(fp, relabelling)
            ep_new = perm_conjugate(ep, relabelling)
            bdry_new = perm_on_list(relabelling, self._bdry, self._n)

            T = (fp_new, ep_new, bdry_new)
            if best is None or T == best:
                best_relabellings.append(relabelling)
                best = T
            elif T < best:
                del best_relabellings[:]
                best_relabellings.append(relabelling)
                best = T

        p0 = perm_invert(best_relabellings[0])
        return [perm_compose(p, p0) for p in best_relabellings]

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
        r"""
        Set labels in a canonical way in its automorphism class.

        EXAMPLES::

            sage: from veerer import *
            sage: from veerer.permutation import perm_random, perm_random_centralizer

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: T = Triangulation(t, mutable=True)
            sage: T.set_canonical_labels()
            sage: S = T.copy(mutable=False)
            sage: for _ in range(10):
            ....:     p = perm_random(24)
            ....:     T.relabel(p)
            ....:     T.set_canonical_labels()
            ....:     assert T == S

        Veering triangulation::

            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = VeeringTriangulation(t, cols, mutable=True)
            sage: T.set_canonical_labels()
            sage: S = T.copy(mutable=False)
            sage: for _ in range(10):
            ....:     p = perm_random(24)
            ....:     T.relabel(p)
            ....:     T.set_canonical_labels()
            ....:     assert T == S

        Veering triangulation with boundary::

            sage: vt = VeeringTriangulation("", boundary="(0:2,1:1,~1:2,~0:1)", colouring="RR", mutable=True)
            sage: vt.set_canonical_labels()
            sage: vt
            VeeringTriangulation("", boundary="(0:2,1:1,~1:2,~0:1)", colouring="RR")

        Veering triangulation with a linear subspace constraint::

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(2, 3, 4, 5, 1, 2)
            sage: L = VeeringTriangulationLinearFamily(T, [s, t], mutable=True)
            sage: L.set_canonical_labels()
            sage: L0 = L.copy(mutable=False)
            sage: for _ in range(10):
            ....:     p = perm_random_centralizer(L.edge_permutation(copy=False))
            ....:     L.relabel(p)
            ....:     L.set_canonical_labels()
            ....:     assert L == L0, (L, L0)

        Examples with boundaries::

            sage: t = Triangulation("(0,1,2)(~0,3,4)", boundary="(~1:1)(~2:2)(~3:3)(~4:4)", mutable=True)
            sage: t.set_canonical_labels()
            sage: t
            Triangulation("(1,~3,~2)(~4,~1,~0)", boundary="(0:1)(2:4)(3:3)(4:2)")

            sage: t0 = t.copy(mutable=False)
            sage: for _ in range(10):
            ....:     p = perm_random_centralizer(t.edge_permutation(copy=False))
            ....:     t.relabel(p)
            ....:     t.set_canonical_labels()
            ....:     assert t == t0, (t, t0)
        """
        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        r, _ = self.best_relabelling()
        self.relabel(r, check=False)

    def is_isomorphic(self, other, certificate=False):
        r"""
        Check whether ``self`` is isomorphic to ``other``.

        TESTS::

            sage: from veerer import Triangulation, VeeringTriangulation
            sage: from veerer.permutation import perm_random_centralizer

            sage: T = Triangulation("(0,5,1)(~0,4,2)(~1,~2,~4)(3,6,~5)", mutable=True)
            sage: TT = T.copy()
            sage: for _ in range(10):
            ....:     rel = perm_random_centralizer(TT.edge_permutation(False))
            ....:     TT.relabel(rel)
            ....:     assert T.is_isomorphic(TT)

            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: V = VeeringTriangulation(fp, cols, mutable=True)
            sage: W = V.copy()
            sage: p = perm_random_centralizer(V.edge_permutation(copy=False))
            sage: W.relabel(p)
            sage: assert V.is_isomorphic(W) is True
            sage: ans, cert = V.is_isomorphic(W, True)
            sage: V.relabel(cert)
            sage: assert V == W
        """
        if type(self) is not type(other):
            raise TypeError("can only check isomorphisms between identical types")

        if self._non_isom_easy(other):
            return (False, None) if certificate else False

        r1, data1 = self.best_relabelling()
        r2, data2 = other.best_relabelling()

        if data1 != data2:
            return (False, None) if certificate else False
        elif certificate:
            return (True, perm_compose(r1, perm_invert(r2)))
        else:
            return True

