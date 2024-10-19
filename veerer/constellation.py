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

import collections
import numbers
from array import array

from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE, rich_to_bool

from .permutation import (perm_init, perm_check, perm_cycles, perm_dense_cycles,
                          perm_invert, perm_conjugate, perm_cycle_string, perm_cycles_lengths,
                          perm_cycles_to_string, perm_on_list, perm_cycle_type,
                          perm_num_cycles, str_to_cycles, str_to_cycles_and_data, perm_compose, perm_from_base64_str,
                          uint_base64_str, uint_from_base64_str, perm_base64_str,
                          perms_are_transitive, triangulation_relabelling_from)

# TODO: maybe do a class Constellation with boundary?
class Constellation:
    __slots__ = ['_mutable', '_n', '_fp', '_ep', '_vp', '_data']

    def __init__(self, n, vp, ep, fp, data, mutable=False, check=True):
        self._n = n

        if vp is None:
            vp = self._vp = array('i', [-1] * n)
            for i in range(n):
                vp[fp[ep[i]]] = i
        else:
            self._vp = vp

        self._ep = ep
        self._fp = fp

        self._data = data
        self._mutable = mutable
        self._set_data_pointers()

        if check:
            self._check(ValueError)

    def _set_data_pointers(self):
        pass

    def _check(self, error=RuntimeError):
        n = self._n

        if not (hasattr(self, '_vp') and hasattr(self, '_ep') and hasattr(self, '_fp') and hasattr(self, '_data')):
            raise error('missing attributes: these must be _vp, _ep, _fp, _data')
        if not perm_check(self._vp, n):
            raise error('vp is not a permutation: {}'.format(self._vp))
        if not perm_check(self._ep, n):
            raise error('ep is not permutation: {}'.format(self._ep))
        if not perm_check(self._fp, n):
            raise error('fp is not a permutation: {}'.format(self._fp))
        if not perms_are_transitive([self._vp, self._ep, self._fp]):
            raise error('(fp, ep, vp) do not generate a transitive group')
        for l in self._data:
            if not isinstance(l, collections.abc.Sequence) or len(l) != n:
                raise error('each data must be a sequence of same length as the underlying permutations got a {} of length {}'.format(type(l).__name__, len(l)))
            if self._mutable and not isinstance(l, collections.abc.MutableSequence):
                raise error('immutable data in mutable object')

        for i in range(n):
            if self._ep[self._ep[i]] != i:
                raise error('invalid edge permutation at half-edge i={} (vp={} ep={} fp={})'.format(self._half_edge_string(i), self._vp, self._ep, self._fp))
            if self._fp[self._ep[self._vp[i]]] != i:
                raise error('fev relation not satisfied at half-edge i={}'.format(self._half_edge_string(i)))

    def _check_alloc(self, n):
        if len(self._vp) < n or len(self._ep) < n or len(self._fp) < n:
            raise TypeError("reallocation needed")

    def _realloc(self, n_max):
        if n_max < self._n:
            return
        self._vp.extend([-1] * (n_max - self._n))
        self._ep.extend([-1] * (n_max - self._n))
        self._fp.extend([-1] * (n_max - self._n))

    def __getstate__(self):
        r"""
        TESTS::

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph

            sage: t = Triangulation("(0,1,2)")
            sage: dumps(t)  # indirect doctest
            b'x\x9ck`J.KM-J-\xd2+)\xcaL\xccK/\xcdI,\xc9\xcc\xcf\xe3\nA\xe1\x152h6\x162\xc6\x162ix3{3vz3z3y3\x00!\x8cfHM\xd2\x03\x00\xb9\xd6\x15\xd9'

            sage: t = VeeringTriangulation("(0,1,2)", "BBR")
            sage: dumps(t)  # indirect doctest
            b'x\x9ck`J.KM-J-\xd2\x03Q\x99y\xe9\xf1%E\x99\x89y\xe9\xa59\x89%\x99\xf9y\\a\x10\xd1\x10\x14\xc1B\x06\xcd\xc6B\xc6\xd8B&\rofo\xa6NoFo&o\x06 \x84\xd1\x0c@\x9a\xc9\x9b15I\x0f\x00Q\x1f\x1c\xdf'
            sage: t = VeeringTriangulation("(0,1,2)(~0,3,4)", "(~1:1)(~2:1)(~3:1)(~4:1)", "RBRBR")
            sage: dumps(t)  # indirect doctest
            b'x\x9ck`J.KM-J-\xd2\x03Q\x99y\xe9\xf1%E\x99\x89y\xe9\xa59\x89%\x99\xf9y\\a\x10\xd1\x10\x14\xc1B\x06\xcd\xc6B\xc6\xd8B&\ro.o\xa6NoFo&o\x06o\x16oNoVo6ovo\x0eof \x9b\x03\xc8b\x03\x8a\xb0\x00yL@5\x0cH\x90\x11\n\x19\xc0z!\x18\xceJM\xd2\x03\x00\x8a)%{'
        """
        a = [self._n]
        a.append(len(self._data))
        a.append(self._mutable)
        a.extend(self._fp)
        a.extend(self._ep)
        for l in self._data:
            a.extend(l)
        return a

    def __setstate__(self, arg):
        r"""
        TESTS::

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

            sage: from veerer import VeeringTriangulation
            sage: t0 = VeeringTriangulation("(0,1,2)", "BBR", mutable=False)
            sage: t1 = VeeringTriangulation("(0,1,2)", "BBR", mutable=True)
            sage: t2 = VeeringTriangulation("(0,1,2)(~0,3,4)", "(~1:1)(~2:1)(~3:1)(~4:1)", "RBRBR")

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

            sage: s2 = loads(dumps(t2))
            sage: assert s2 == t2
        """
        # We do not know how many slots we have in data
        n = self._n = arg[0]
        k = arg[1]  # length of data
        self._mutable = arg[2]
        self._fp = array('i', arg[3 : n + 3])
        self._ep = array('i', arg[n + 3 : 2 * n + 3])
        data = []
        for i in range(2, k + 2):
            data.append(array('i', arg[i * n + 3 : (i + 1) * n + 3]))

        self._vp = array('i', [-1] * n)
        for i in range(n):
            self._vp[self._fp[self._ep[i]]] = i

        self._data = tuple(data)
        self._set_data_pointers()

    def set_immutable(self):
        self._mutable = False

    def __hash__(self):
        r"""
        TESTS::

            sage: from itertools import permutations, combinations
            sage: from veerer import Triangulation, VeeringTriangulation

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

            sage: triangulations = []
            sage: for cols in ["RRB", "RBR", "BRR", "BBR", "BRB", "RBB"]:
            ....:     t = VeeringTriangulation("(0,1,2)", cols)
            ....:     triangulations.append(t)

            sage: for i in range(len(triangulations)):
            ....:     for j in range(len(triangulations)):
            ....:         assert (triangulations[i] == triangulations[j]) == (i == j), (i, j)
            ....:         assert (triangulations[i] != triangulations[j]) == (i != j), (i, j)

            sage: hashes1 = {}
            sage: hashes2 = {}
            sage: for t in triangulations:
            ....:     h1 = hash(t) % (2 ** 16)
            ....:     h2 = (hash(t) >> 16) % (2 ** 16)
            ....:     if h1 in hashes1:
            ....:         print('collision 1: {} {}'.format(hashes1[h1], t))
            ....:     else:
            ....:         hashes1[h1] = t
            ....:     if h2 in hashes2:
            ....:         print('collision 2: {} {}'.format(hashes2[h2], t))
            ....:     else:
            ....:         hashes2[h2] = t
            sage: assert len(hashes1) == len(hashes2) == len(triangulations), (len(hashes1), len(hashes2), len(triangulations))

            sage: t = VeeringTriangulation("(0,1,2)", cols, mutable=True)
            sage: hash(t)
            Traceback (most recent call last):
            ...
            ValueError: mutable veering triangulation not hashable
        """
        if self._mutable:
            raise ValueError('mutable veering triangulation not hashable')

        x = 140737488617563
        x = ((x ^ hash(self._vp.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._ep.tobytes())) * 2147483693) + 82520 + self._n + self._n

        for l in self._data:
            x = ((x ^ hash(l.tobytes())) * 2147483693) + 82520 + self._n + self._n

        return x

    def _check_half_edge(self, e):
        if not isinstance(e, numbers.Integral):
            raise TypeError('invalid half-edge {}'.format(e))
        e = int(e)
        if e < 0 or e >= self._n:
            raise ValueError('half-edge number out of range e={}'.format(e))
        return e

    def to_string(self):
        r"""
        Serialize this triangulation as a string.

        EXAMPLES::

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph

            sage: Triangulation("(0,1,2)(~0,~1,~2)").to_string()
            '6_354102_543210_120534_000000'
            sage: Triangulation("(0,1,2)", boundary="(~0:1)(~1:1,~2:1)").to_string()
            '6_354120_543210_120435_000111'

            sage: VeeringTriangulation("(0,1,2)", "RRB").to_string()
            '3_201_012_120_000_112'

            sage: StrebelGraph("(0,1,2)(~0,~1:1,~2:2)").to_string()
            '6_354102_543210_120534_000210'
        """
        return uint_base64_str(self._n) + '_' + perm_base64_str(self._vp) + '_' + perm_base64_str(self._ep) + '_' + perm_base64_str(self._fp) + '_' + '_'.join(perm_base64_str(l) for l in self._data)

    def from_face_edge_perms(self, fp, ep, data=(), mutable=False, check=True):
        raise ValueError

    @classmethod
    def from_permutations(cls, vp, ep, fp, data=(), mutable=False, check=True):
        r"""
        INPUT:

        - ``vp``, ``ep``, ``fp`` -- the vertex, edge and face permutations

        - ``data``

        - ``check`` - boolean (default: ``True``) - if set to ``False`` no
          check are performed

        EXAMPLES::

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph
            sage: from array import array

            sage: vp = array('i', [2, 8, 7, 0, 3, 1, 5, 6, 4])
            sage: ep = array('i', [8, 7, 2, 3, 4, 5, 6, 1, 0])
            sage: fp = array('i', [1, 2, 0, 4, 8, 6, 7, 5, 3])
            sage: Triangulation.from_permutations(vp, ep, fp, (array('i', [0] * 9),))
            Triangulation("(0,1,2)(3,4,~0)(5,6,~1)")
            sage: Triangulation.from_permutations(None, ep, fp, (array('i', [0] * 9),))
            Triangulation("(0,1,2)(3,4,~0)(5,6,~1)")
            sage: Triangulation.from_permutations(vp, None, fp, (array('i', [0] * 9),))
            Triangulation("(0,1,2)(3,4,~0)(5,6,~1)")
            sage: Triangulation.from_permutations(vp, ep, None, (array('i', [0] * 9),))
            Triangulation("(0,1,2)(3,4,~0)(5,6,~1)")

            sage: vp = array('i', [1, 3, 0, 2])
            sage: ep = array('i', [3, 2, 1, 0])
            sage: StrebelGraph.from_permutations(vp, ep, None, data=(array('i', [1, 0, 0, 1]),))
            StrebelGraph("(0:1,1,~0:1,~1)")
        """
        if (vp is None) + (ep is None) + (fp is None) > 1:
            raise ValueError('at most one of vp, ep, fp could be None')

        C = cls.__new__(cls)
        if vp is not None:
            n = len(vp)
        elif ep is not None:
            n = len(ep)

        if vp is None:
            vp = array('i', [-1] * n)
            for i in range(n):
                vp[fp[ep[i]]] = i
        elif ep is None:
            ep = array('i', [-1] * n)
            for i in range(n):
                ep[vp[fp[i]]] = i
        elif fp is None:
            fp = array('i', [-1] * n)
            for i in range(n):
                fp[ep[vp[i]]] = i

        C._n = n
        C._vp = vp
        C._ep = ep
        C._fp = fp
        C._mutable = mutable
        C._data = data
        C._set_data_pointers()

        if check:
            C._check(ValueError)

        return C

    @classmethod
    def from_string(cls, s, mutable=False, check=True):
        r"""
        Deserialization from string.

        EXAMPLES::

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph

            sage: T = Triangulation("(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: Triangulation.from_string(T.to_string()) == T
            True
        """
        parts = s.split('_')
        n = parts[0]
        vp = parts[1]
        ep = parts[2]
        fp = parts[3]
        data = parts[4:]
        n = uint_from_base64_str(n)
        vp = perm_from_base64_str(vp, n)
        ep = perm_from_base64_str(ep, n)
        fp = perm_from_base64_str(fp, n)
        data = tuple(perm_from_base64_str(ss, n) for ss in data)
        return cls.from_permutations(vp, ep, fp, data, mutable, check)

    def __eq__(self, other):
        r"""
        Return whether ``self`` and ``other`` are equal.

        EXAMPLES::

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph

            sage: Triangulation("(0,1,2)(~0,~1,~2)") == Triangulation("(0,1,2)(~0,~1,~2)")
            True
            sage: Triangulation("(0,1,2)(~0,~1,~2)") == Triangulation("(0,~0,1)(~1,2,~2)")
            False

            sage: StrebelGraph("(0,1,2,~0,~1,~2)") == StrebelGraph("(0,1,2,~0,~1,~2)")
            True
            sage: StrebelGraph("(0,1,2,~0,~1,~2)") == StrebelGraph("(0,1:1,2,~0,~1,~2)")
            False
        """
        if type(self) != type(other):
            raise TypeError
        return self._n == other._n and self._fp == other._fp and self._ep == other._ep and self._data == other._data

    def __ne__(self, other):
        r"""
        Return whether ``self`` and ``other`` are different.

        EXAMPLES::

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph

            sage: Triangulation("(0,1,2)(~0,~1,~2)") != Triangulation("(0,1,2)(~0,~1,~2)")
            False
            sage: Triangulation("(0,1,2)(~0,~1,~2)") != Triangulation("(0,~0,1)(~1,2,~2)")
            True

            sage: StrebelGraph("(0,1,2,~0,~1,~2)") != StrebelGraph("(0,1,2,~0,~1,~2)")
            False
            sage: StrebelGraph("(0,1,2,~0,~1,~2)") != StrebelGraph("(0,1:1,2,~0,~1,~2)")
            True
        """
        if type(self) != type(other):
            raise TypeError
        return self._n != other._n or self._fp != other._fp or self._ep != other._ep or self._data != other._data

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other`` according to the operator ``op``.
        """
        if type(self) != type(other):
            raise TypeError

        c = (self._n > other._n) - (self._n < other._n)
        if c:
            return rich_to_bool(op, c)

        c = (self._fp > other._fp) - (self._fp < other._fp)
        if c:
            return rich_to_bool(op, c)

        c = (self._ep > other._ep) - (self._ep < other._ep)
        if c:
            return rich_to_bool(op, c)

        c = (self._data > other._data) - (self._data < other._data)
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

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], "RRB", mutable=True)
            sage: S1 = T.copy()
            sage: S2 = T.copy()
            sage: T == S1 == S2
            True
            sage: S1.flip(1,BLUE)
            sage: T == S1
            False
            sage: T == S2
            True

        TESTS::

            sage: from veerer import Triangulation
            sage: T = Triangulation("(0,1,2)(~0,~1,~2)", mutable=True)
            sage: U = T.copy(mutable=False)
            sage: _ = hash(U)

            sage: from veerer import VeeringTriangulation
            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB", mutable=True)
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
                T._data = self._data
                T._mutable = mutable
        else:
            T = cls.__new__(cls)
            T._n = self._n
            T._fp = self._fp[:]
            T._ep = self._ep[:]
            T._vp = self._vp[:]
            T._data = tuple(l[:] for l in self._data)
            T._mutable = mutable

        T._set_data_pointers()
        return T

    def vertex_permutation(self, copy=True):
        if copy:
            return self._vp[:]
        else:
            return self._vp

    def next_at_vertex(self, e, check=True):
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
            e = self._check_half_edge(e)
        return self._vp[e]

    def previous_at_vertex(self, e, check=True):
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
            e = self._check_half_edge(e)
        return self._fp[self._ep[e]]

    def edge_permutation(self, copy=True):
        if copy:
            return self._ep[:]
        else:
            return self._ep

    def next_in_edge(self, e, check=True):
        if check:
            self._check_half_edge(e)
        return self._ep[e]

    def previous_in_edge(self, e, check=True):
        if check:
            self._check_half_edge(e)
        return self._vp[self._fp[e]]

    def face_permutation(self, copy=True):
        if copy:
            return self._fp[:]
        else:
            return self._fp

    def next_in_face(self, e, check=True):
        if check:
            e = self._check_half_edge(e)
        return self._fp[e]

    def previous_in_face(self, e, check=True):
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
            e = self._check_half_edge(e)
        return self._ep[self._vp[e]]

    def boundary_vector(self, copy=True):
        if copy:
            return self._data[0][:]
        else:
            return self._data[0]

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
            1
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

        for l in self._data:
            l[e], l[E] = l[E], l[e]

    def relabel(self, p, check=True):
        r"""
        Relabel this triangulation inplace according to the permutation ``p``.

        EXAMPLES::

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph, BLUE, RED

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

            sage: T0 = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)")
            sage: T = T0.copy(mutable=True)
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

        Veering triangulations::

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RBB", mutable=True)
            sage: T.relabel([0,1,3,2,5,4])
            sage: T
            VeeringTriangulation("(0,1,~2)(2,~0,~1)", "RBB")
            sage: T._check()

        Composing relabellings and permutation composition::

            sage: from veerer.permutation import perm_compose, perm_random_centralizer
            sage: fp = "(0,16,~15)(1,19,~18)(2,22,~21)(3,21,~20)(4,20,~19)(5,23,~22)(6,18,~17)(7,17,~16)(8,~1,~23)(9,~2,~8)(10,~3,~9)(11,~4,~10)(12,~5,~11)(13,~6,~12)(14,~7,~13)(15,~0,~14)"
            sage: cols = "RRRRRRRRBBBBBBBBBBBBBBBB"
            sage: T0 = VeeringTriangulation(fp, cols)
            sage: for _ in range(10):
            ....:     p1 = perm_random_centralizer(T0.edge_permutation(copy=False))
            ....:     p2 = perm_random_centralizer(T0.edge_permutation(copy=False))
            ....:     T1 = T0.copy(mutable=True)
            ....:     T1.relabel(p1)
            ....:     T1.relabel(p2)
            ....:     T2 = T0.copy(mutable=True)
            ....:     T2.relabel(perm_compose(p1, p2))
            ....:     assert T1  == T2

        TESTS:

        This example used to be wrong::

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE], mutable=True)
            sage: T.relabel([1,5,0,2,4,3])
            sage: T.edge_colour(0) == BLUE
            True
            sage: T.edge_colour(1) == RED
            True
            sage: T._check()

            sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE], mutable=True)
            sage: from veerer.permutation import perm_random
            sage: for _ in range(10):
            ....:     r = perm_random(6)
            ....:     T.relabel(r)
            ....:     T._check()
        """
        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        n = self._n
        if check and not perm_check(p, n):
            # if the input is not a valid permutation, we assume that half-edges
            # are not separated
            p = perm_init(p, n, self._ep)
            if not perm_check(p, n):
                raise ValueError('invalid relabeling permutation')

        # TODO: would better be inplace!!
        self._vp = perm_conjugate(self._vp, p)
        self._ep = perm_conjugate(self._ep, p)
        self._fp = perm_conjugate(self._fp, p)
        for l in self._data:
            perm_on_list(p, l, n)

        # TODO: remove check
        self._check()

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
            sage: T = Triangulation.from_permutations(vp, ep, fp, (array('i', [0]*6),), mutable=True)
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
            sage: T = Triangulation.from_permutations(vp, ep, fp, (array('i', [0]*3),), mutable=True)
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
        return triangulation_relabelling_from(self._vp, self._ep, start_edge)

    def _automorphism_good_starts(self):
        # TODO: should we try to discriminate using vp, ep, fp, data?
        return range(self._n)

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
        best_relabellings = self.best_relabelling(all=True)[0]
        p0 = perm_invert(best_relabellings[0])
        return [perm_compose(p, p0) for p in best_relabellings]

    def best_relabelling(self, all=False):
        r"""
        Return a pair ``(r, data)`` where ``r`` is a relabelling that
        brings this constellation to the canonical one.

        EXAMPLES::

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph
            sage: from veerer.permutation import perm_random_centralizer

            sage: examples = []
            sage: triangles = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: examples.append(Triangulation(triangles, mutable=True))
            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: examples.append(VeeringTriangulation(fp, cols, mutable=True))
            sage: examples.append(StrebelGraph("(0,6,~5,~3,~1,4,~4:3,2,~2:3)(1:2)(3:2,~0)(5:2)(~6)", mutable=True))

            sage: for G in examples:
            ....:     print(G)
            ....:     r, (fp, ep, data) = G.best_relabelling()
            ....:     for _ in range(10):
            ....:         p = perm_random_centralizer(G.edge_permutation(copy=False))
            ....:         G.relabel(p)
            ....:         r2, (fp2, ep2, data2) = G.best_relabelling()
            ....:         assert fp2 == fp, (G, fp, fp2)
            ....:         assert ep2 == ep, (G, ep, ep2)
            ....:         assert data2 == data, (G, data, data2)
            Triangulation("(0,~1,2)(1,~3,~0)(3,4,~5)(5,~9,~7)(6,~2,~4)(7,~6,8)(9,10,~11)(11,~8,~10)")
            VeeringTriangulation("(0,~1,2)(1,~3,~0)(3,4,~5)(5,~9,~7)(6,~2,~4)(7,~6,8)(9,10,~11)(11,~8,~10)", "BRBBBRRBBBBR")
            StrebelGraph("(0,6,~5,~3,~1,4,~4:3,2,~2:3)(1:2)(3:2,~0)(5:2)(~6)")
        """
        n = self._n
        fp = self._fp
        ep = self._ep

        best = None
        if all:
            relabellings = []

        for start_edge in range(self._n):
            relabelling = self._relabelling_from(start_edge)

            fp_new = perm_conjugate(fp, relabelling)
            ep_new = perm_conjugate(ep, relabelling)
            data_new = [l[:] for l in self._data]
            for l in data_new:
                perm_on_list(relabelling, l, self._n)

            T = (fp_new, ep_new, data_new)
            if best is None or T < best:
                best_relabelling = relabelling
                best = T
                if all:
                    del relabellings[:]
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
            sage: T
            Triangulation("(0,10,~9)(1,~8,9)(2,~6,7)(3,~4,5)(4,~3,~11)(6,~2,~5)(8,~1,~7)(11,~10,~0)")
            sage: T.set_canonical_labels()
            sage: T
            Triangulation("(0,1,8)(2,11,~3)(3,~11,~4)(4,10,~5)(5,~10,~6)(6,9,~7)(7,~9,~8)(~2,~1,~0)")
        """
        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        r, _ = self.best_relabelling()
        self.relabel(r, check=False)

    def iso_sig(self):
        r"""
        Return a canonical signature.

        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(0,3,1)(~0,4,2)(~1,~2,~4)")
            sage: T.iso_sig()
            '9_583764021_876543210_134052786_000000000'
            sage: TT = Triangulation.from_string(T.iso_sig())
            sage: TT
            Triangulation("(0,1,3)(2,4,~3)(~2,~1,~0)")
            sage: TT.iso_sig() == T.iso_sig()
            True

            sage: T = Triangulation("(0,10,~6)(1,12,~2)(2,14,~3)(3,16,~4)(4,~13,~5)(5,~1,~0)(6,~17,~7)(7,~14,~8)(8,13,~9)(9,~11,~10)(11,~15,~12)(15,17,~16)")
            sage: T.iso_sig()
            'A_f23456y89azcxesvromwphuiqbjl0dkgnt71_zyxwvutsrqponmlkjihgfedcba9876543210_a6cjfmxegokhwiruqnlv0dtbp987y54321zs_000000000000000000000000000000000000'
            sage: Triangulation.from_string(T.iso_sig())
            Triangulation("(0,10,~15)(1,6,~2)(2,12,~3)(3,~16,~4)(4,15,~5)(5,~13,~6)(7,14,~8)(8,16,~9)(9,~11,~10)(11,17,~12)(13,~17,~14)(~7,~1,~0)")

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = VeeringTriangulation(t, cols, mutable=True)
            sage: T.iso_sig()
            'o_fn345678mhjlkig9eadbc021_nmlkjihgfedcba9876543210_18bcad9e0gikjhf765432mnl_000000000000000000000000_122121212222222212121221'

        If we relabel the triangulation, the isomorphic signature does not change::

            sage: from veerer.permutation import perm_random
            sage: p = perm_random(24)
            sage: T.relabel(p)
            sage: T.iso_sig()
            'o_fn345678mhjlkig9eadbc021_nmlkjihgfedcba9876543210_18bcad9e0gikjhf765432mnl_000000000000000000000000_122121212222222212121221'

        An isomorphic triangulation can be reconstructed from the isomorphic
        signature via::

            sage: s = T.iso_sig()
            sage: T2 = VeeringTriangulation.from_string(s)
            sage: T == T2
            False
            sage: T.is_isomorphic_to(T2)
            True

        TESTS::

            sage: from veerer.veering_triangulation import VeeringTriangulation
            sage: from veerer.permutation import perm_random

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = VeeringTriangulation(t, cols, mutable=True)
            sage: iso_sig = T.iso_sig()
            sage: for _ in range(10):
            ....:     p = perm_random(24)
            ....:     T.relabel(p)
            ....:     assert T.iso_sig() == iso_sig

            sage: VeeringTriangulation("(0,1,2)(3,4,~1)(5,6,~4)", "RBGGRBG").iso_sig()
            '9_523706841_872345610_631270584_000000000_218281812'
        """
        T = self.copy(mutable=True)
        T.set_canonical_labels()
        return T.to_string()

    def _non_isom_easy(self, other):
        r"""
        A quick certificate of non-isomorphism that does not require relabellings.
        """
        return (perm_cycle_type(self._vp) != perm_cycle_type(other._vp) or
            perm_cycle_type(self._ep) != perm_cycle_type(other._ep) or
            perm_cycle_type(self._fp) != perm_cycle_type(other._fp) or
            any(sorted(l_self) != sorted(l_other) for l_self, l_other in zip(self._data, other._data)))

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
