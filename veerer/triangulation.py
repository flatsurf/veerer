r"""
Triangulation of surfaces.
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

import numbers
from array import array

from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE, rich_to_bool

from .permutation import (perm_init, perm_check, perm_cycles, perm_dense_cycles,
                          perm_invert, perm_conjugate, perm_cycle_string, perm_cycles_lengths,
                          perm_num_cycles, str_to_cycles, perm_compose, perm_from_base64_str,
                          uint_base64_str, uint_from_base64_str, perm_base64_str,
                          perms_are_transitive, triangulation_relabelling_from)


def face_edge_perms_init(data):
    r"""
    EXAMPLES::

        sage: from veerer.triangulation import face_edge_perms_init  # random output due to deprecation warnings from realalg

        sage: face_edge_perms_init('(0,1,2)(~0,~1,~2)')
        (array('i', [1, 2, 0, 5, 3, 4]), array('i', [5, 4, 3, 2, 1, 0]))

        sage: face_edge_perms_init('(0,1,2)')
        (array('i', [1, 2, 0]), array('i', [0, 1, 2]))

    TESTS:

    Check that edge permutation do not depend on the details of faces::

        sage: f1 = "(0,~5,4)(3,5,6)(1,2,~6)"
        sage: f2 = "(0,6,5)(1,2,~6)(3,4,~5)"
        sage: f3 = "(6,4,3)(~6,~5,1)(5,0,2)"
        sage: assert face_edge_perms_init(f1)[1] == face_edge_perms_init(f2)[1] == face_edge_perms_init(f3)[1]
    """
    if isinstance(data, str):
        l = str_to_cycles(data)
    else:
        l = [[int(i) for i in c] for c in data]

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

# NOTE: we don't really care that we have a triangulation here. When
# there is no restriction on the cycle decomposition of _fp we got
# a coding for any cell decomposition (with possibly folded edges).
# flipping here still makes sense
#
#  +                 +                 +               _
#   \b             a/                   \            _/ /
#    \      e      /                     \      e __/  /
#     +-----------+       --------->      +    __/
#    /             \                     /  __/        \
#   /c             d\                   /__/            \
#  +                 +                 +/
#
# and we have operations that consist in removing/adding edges.
# For that purpose, it would be convenient to allow partial
# permutations of {0, 1, ..., n-1}.
#
# TODO: implement: is_removable(), remove(), add()
#
# NOTE: Much of the above is already done in the flatsurf package
# with a quite different encoding, see
#
#    https://github.com/flatsurf/sage-flatsurf


class Triangulation(object):
    r"""
    A triangulation of an oriented surface

    with attributes:

    * _n  number of half-edges (an int)
    * _fp  face permutation (an array)
    * _ep  edge permutation (an array)
    * _vp  vertex permutation (an array)

    Each of fp, ep, vp is a permutation of oriented half-edges (aka
    "darts"), given as a function.  Our conventions for the
    permutations are set out in the following figure::

              ~b
            ----->
        w-----------*-----------v
         \             <-----  /
          \ \             b   / /
           \ \ c             / /~a
            \ \     F       / /
             \ v           / v
              *         ^ *
             ^ \       / /
              \ \   a / /
             ~c\ \   / /
                \ \   /
                   \ /
                    u

    Here the face permutation sends a to b, b to c, and c to a.  The
    vertex permutation send a to ~c, b to ~a, and c to ~b.  The edge
    permutation interchanges e and ~e for every edge; the edge e is
    folded if and only if e = ~e.  Thus folded edges are fixed by the
    edge permutation.

    EXAMPLES::

        sage: from veerer import *

        sage: T = Triangulation("(~2, 1, ~0)(~1, 0, 2)", mutable=True)
        sage: T.genus()
        1
        sage: T.num_faces()
        2
        sage: T.num_vertices()
        1
        sage: T.flip(0)
        sage: T
        Triangulation("(0,~1,~2)(1,2,~0)")
        sage: T.flip(0)
        sage: T
        Triangulation("(0,~2,1)(2,~1,~0)")
        sage: T.flip(0)
        sage: T
        Triangulation("(0,1,2)(~2,~0,~1)")
        sage: T.flip(0)
        sage: T
        Triangulation("(0,2,~1)(1,~0,~2)")

        sage: T.set_immutable()
        sage: T.flip(0)
        Traceback (most recent call last):
        ...
        ValueError: immutable triangulation; use a mutable copy instead

    The surface must be connected::

        sage: Triangulation("(0,1,2)(3,4,5)")
        Traceback (most recent call last):
        ...
        ValueError: (fp, ep, vp) do not generate a transitive group
    """
    __slots__ = ['_mutable', '_n', '_fp', '_ep', '_vp']

    def __init__(self, triangles, mutable=False, check=True):
        if isinstance(triangles, Triangulation):
            self._fp = triangles.face_permutation(copy=True)
            self._ep = triangles.edge_permutation(copy=True)
        else:
            self._fp, self._ep = face_edge_perms_init(triangles)

        self._mutable = mutable

        fp = self._fp
        ep = self._ep
        n = self._n = len(fp)

        # TODO: face labels are disabled for now. We should actually
        # make it a *partition* of the faces (two faces are allowed to
        # have the same label) The trivial partition {0,1,2,...,F-1}
        # would correspond to no labels and the atomic
        # {{0},{1},...,{F-1}} would correspond to everybody being
        # labelled
        # fl = self._fl = [None] * F # face labels

        # TODO: edge labels are disabled for now.
        # el = self._fl = [None] * E # edge labels

        vp = self._vp = array('i', [-1] * n)
        for i in range(n):
            vp[fp[ep[i]]] = i

        # TODO: vertex labels are disabled for now
        # vl = self._vl = perm_cycles(vp)[1]....

        if check:
            self._check(ValueError)

    def __getstate__(self):
        r"""
        TESTS::

            sage: from veerer import Triangulation
            sage: t = Triangulation("(0,1,2)")
            sage: dumps(t)  # indirect doctest
            b'x\x9ck`J.KM-J-\xd2+)\xcaL\xccK/\xcdI,\xc9\xcc\xcf\xe3\nA\xe1\x152h6\x162\xc6\x162ix3z3y3\x00!\x90\xeeLM\xd2\x03\x00\xe7\x1e\x14^'
        """
        a = list(self._fp)
        a.extend(self._ep)
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
        n = (len(arg) - 1) // 2
        self._n = n
        self._fp = array('i', arg[:n])
        self._ep = array('i', arg[n:2*n])
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

    @staticmethod
    def from_face_edge_perms(fp, ep, vp=None, mutable=False, check=True):
        r"""
        INPUT:

        - ``fp``, ``ep``, ``vp`` -- the face, edge and vertex permutation

        - ``check`` - boolean (default: ``True``) - if set to ``False`` no
          check are performed

        EXAMPLES::

            sage: from veerer import Triangulation
            sage: from array import array

            sage: fp = array('i', [1, 2, 0, 4, 8, 6, 7, 5, 3])
            sage: ep = array('i', [8, 7, 2, 3, 4, 5, 6, 1, 0])
            sage: vp = array('i', [2, 8, 7, 0, 3, 1, 5, 6, 4])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            Triangulation("(0,1,2)(3,4,~0)(5,6,~1)")
        """
        T = Triangulation.__new__(Triangulation)
        n = T._n = len(fp)
        T._fp = fp
        T._ep = ep
        if vp is None:
            fp = T._fp
            ep = T._ep
            vp = array('i', [-1] * n)
            for i in range(n):
                vp[fp[ep[i]]] = i
        T._vp = vp

        T._mutable = mutable

        if check:
            T._check(ValueError)

        return T

    @staticmethod
    def from_string(s, mutable=False, check=True):
        r"""
        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: Triangulation.from_string(T.to_string()) == T
            True
        """
        n, fp, ep = s.split('_')
        n = uint_from_base64_str(n)
        fp = perm_from_base64_str(fp, n)
        ep = perm_from_base64_str(ep, n)
        return Triangulation.from_face_edge_perms(fp, ep, mutable=mutable, check=check)

    def _check(self, error=RuntimeError):
        r"""
        TESTS::

            sage: from veerer import Triangulation

            sage: Triangulation("(0,1,3)")
            Traceback (most recent call last):
            ...
            ValueError: missing edge label 2

            sage: Triangulation("(0,1,~2)")
            Traceback (most recent call last):
            ...
            ValueError: inconsistent permutation data

            sage: Triangulation("(0)")
            Traceback (most recent call last):
            ...
            ValueError: broken face permutation

            sage: from array import array
            sage: fp = array('i', [1,2,0])
            sage: ep = array('i', [0,1,2])
            sage: vp = array('i', [1,2,0])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            Traceback (most recent call last):
            ...
            ValueError: fev relation not satisfied

            sage: fp = array('i', [1,2,0])
            sage: ep = array('i', [1,2,0])
            sage: vp = array('i', [1,2,0])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            Traceback (most recent call last):
            ...
            ValueError: broken edge permutation
        """
        n = self._n

        if not perm_check(self._fp, n):
            raise error('fp is not a permutation: {}'.format(self._fp))
        if not perm_check(self._ep, n):
            raise error('ep is not permutation: {}'.format(self._ep))
        if not perm_check(self._vp, n):
            raise error('vp is not a permutation: {}'.format(self._vp))
        if not perms_are_transitive([self._fp, self._ep, self._vp]):
            raise error('(fp, ep, vp) do not generate a transitive group')

        # The face, edge, and vertex permutations fp, ep, and vp must
        # satisfy the relations fp^3 = ep^2 = fp.ep.vp = Id.  Note
        # that this is a representation of PSL(2, \ZZ) \isom S^2(2, 3,
        # \infty) into the symmetric group on the half edges.

        for i in range(n):
            # NOTE: this first test is relevant only if we enforce triangulations
            # when we generalize to CW complex we can simply remove it
            if self._fp[i] == i or self._fp[self._fp[self._fp[i]]] != i:
                raise error('broken face permutation')
            if self._ep[self._ep[i]] != i:
                raise error('broken edge permutation')
            if self._fp[self._ep[self._vp[i]]] != i:
                raise error('fev relation not satisfied')

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
        return self._n == other._n and self._fp == other._fp and self._ep == other._ep

    def __ne__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n != other._n or self._fp != other._fp or self._ep != other._ep

    def _richcmp_(self, other, op):
        c = (self._n > other._n) - (self._n < other._n)
        if c:
            return rich_to_bool(op, c)

        c = (self._fp > other._fp) - (self._fp < other._fp)
        if c:
            return rich_to_bool(op, c)

        c = (self._ep > other._ep) - (self._ep < other._ep)
        return rich_to_bool(op, c)

    def __lt__(self, other):
        return self._richcmp_(other, op_LT)

    def __le__(self, other):
        return self._richcmp_(other, op_LE)

    def __gt__(self, other):
        return self._richcmp_(other, op_GT)

    def __ge__(self, other):
        return self._richcmp_(other, op_GE)

    def copy(self, mutable=None):
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

        if not self._mutable and not mutable:
            # avoid copies of immutable objects
            if type(self) is Triangulation:
                return self
            else:
                T = Triangulation.__new__(Triangulation)
                T._n = self._n
                T._fp = self._fp
                T._ep = self._ep
                T._vp = self._vp
                T._mutable = mutable

                return T
        else:
            T = Triangulation.__new__(Triangulation)
            T._n = self._n
            T._fp = self._fp[:]
            T._ep = self._ep[:]
            T._vp = self._vp[:]
            T._mutable = mutable

            return T

    def to_flipper(self):
        r"""
        Return the corresponding flipper triangulation

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.to_flipper()  # optional - flipper
            [(~2, ~0, ~1), (0, 1, 2)]

        Conversion of veering triangulation forgot the colors::

            sage: V = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: V.to_flipper()  # optional - flipper
            [(~2, ~0, ~1), (0, 1, 2)]
        """
        from .features import flipper_feature
        flipper_feature.require()

        import flipper
        ep = self._ep
        F = []
        for f in self.faces():
            face = []
            for e in f:
                if ep[e] == e:
                    raise ValueError("flipper do not accept folded edges")
                if ep[e] < e:
                    face.append(~int(ep[e]))
                else:
                    face.append(int(e))
            F.append(tuple(face))

        return flipper.create_triangulation(F)

    def to_curver(self):
        r"""
        Return the corresponding curver triangulation

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.to_curver()  # optional - curver
            [(~2, ~0, ~1), (0, 1, 2)]

        Conversion of veering triangulation forgot the colors::

            sage: V = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: V.to_curver()  # optional - curver
            [(~2, ~0, ~1), (0, 1, 2)]
        """
        from .features import curver_feature
        curver_feature.require()

        import curver
        ep = self._ep
        F = []
        for f in self.faces():
            face = []
            for e in f:
                if ep[e] == e:
                    raise ValueError("curver do not accept folded edges")
                if ep[e] < e:
                    face.append(~int(ep[e]))
                else:
                    face.append(int(e))
            F.append(tuple(face))

        return curver.create_triangulation(F)

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

    def previous_in_face(self, i):
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
        return self._ep[self._vp[i]]

    def num_half_edges(self):
        return self._n

    def folded_edges(self):
        n = self._n
        ep = self._ep
        return [i for i in range(n) if ep[i] == i]

    def num_folded_edges(self):
        n = self._n
        ep = self._ep
        return sum(ep[i] == i for i in range(n))

    def num_edges(self):
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

    def num_faces(self):
        return perm_num_cycles(self._fp)

    def faces(self):
        r"""
        Return the list of faces as triples of half-edges

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.faces()
            [[0, 1, 2], [3, 4, 5], [6, 8, 7]]
        """
        return perm_cycles(self._fp, True, self._n)

    def edges(self):
        r"""
        Return the list of faces as tuples of half-edges

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.edges()
            [[0, 8], [1], [2], [3, 7], [4], [5], [6]]
        """
        return perm_cycles(self._ep, True, self._n)

    def vertices(self):
        r"""
        Return the list of faces as tuples of half-edges

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.vertices()
            [[0, 2, 1, 8, 6, 3, 5, 4, 7]]
        """
        return perm_cycles(self._vp, True, self._n)

    def num_vertices(self):
        return perm_num_cycles(self._vp)

    def euler_characteristic(self):
        r"""
        Return the Euler characteristic of this triangulation.

        EXAMPLES::

            sage: from veerer import Triangulation

        A sphere::

            sage: T = Triangulation("(0,1,2)")
            sage: T.euler_characteristic()
            2

        A torus::

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.euler_characteristic()
            0

        A genus 2 surface::

            sage: T = Triangulation("(0,1,2)(~2,3,4)(~4,5,6)(~6,~0,7)(~7,~1,8)(~8,~3,~5)")
            sage: T.euler_characteristic()
            -2
        """
        return self.num_faces() - self.num_edges() + (self.num_vertices() + self.num_folded_edges())

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
        # chi = 2 - 2g
        return (2 - self.euler_characteristic()) // 2

    def __str__(self):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: str(T)
            'Triangulation("(0,1,2)(~2,~0,~1)")'
        """
        return 'Triangulation("%s")' % perm_cycle_string(self._fp, n=self._n, involution=self._ep)

    def __repr__(self):
        return str(self)

    def _check_homology_matrix(self, m):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: m = matrix([[1,1],[-1,0],[0,-1]])
            sage: T._check_homology_matrix(m)
        """
        ne = self.num_edges()
        ep = self._ep

        assert m.nrows() == ne
        try:
            V = m.row_ambient_module()
        except AttributeError:
            V = m._row_ambient_module()

        # boundary condition
        for F in self.faces():
            v = V.zero()
            for e in F:
                E = ep[e]
                if E > e:
                    v += m[e]
                elif E < e:
                    v -= m[E]
                else:
                    # folded edge condition
                    # NOTE: it is stupid to keep all these zeros in
                    # the matrix! We should label the folded edges
                    # after the non-folded ones...
                    # though we can allow non-zero stuff assuming
                    # that the half edge is oriented toward the pole
                    assert m[e].is_zero()
            assert v.is_zero()

    # TODO: also compute the intersection pairing!
    def homology_matrix(self):
        r"""
        Return a basis of homology as a matrix.

        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)")
            sage: T.homology_matrix().transpose().echelon_form()
            [ 1  0  1  0  1  1]
            [ 0  1  0  0 -1  0]
            [ 0  0  0  1  0 -1]

            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,6,0)")
            sage: T.homology_matrix().transpose().echelon_form()
            [ 1  0  1  0  1  1  0]
            [ 0  1  0  0 -1  0  0]
        """
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ

        ep = self._ep
        nf = self.num_faces()
        ne = self.num_edges()
        nfe = self.num_folded_edges()
        m = matrix(ZZ, nf + nfe, ne)

        # face equations
        for i, f in enumerate(self.faces()):
            for e in f:
                if ep[e] == e:
                    continue
                elif ep[e] < e:
                    m[i, ep[e]] -= 1
                else:
                    m[i, e] += 1

        # force the folded edge to have coefficient zero
        for e in range(ne):
            if ep[e] == e:
                m[i, e] = 1
                i += 1

        # compute and check
        h = m.right_kernel_matrix().transpose()
        self._check_homology_matrix(h)
        return h

    def flip_homological_action(self, e, m, twist=False):
        r"""
        Multiply the matrix ``m`` on the left by the homology action of
        the ``e``-flip.

        The matrix ``m`` must have ``ne`` rows and each column represents a
        vector in cohomology (possibly twisted for quadratic differentials).
        That is to say, each column has to satisfy the following conditions:

        - for each face `F` the evaluation on its boundary is zero (up to twist).

        - for each folded edge, the corresponding entry is zero (up to twist).

        INPUT:

        - ``e`` - a half edge

        - ``m`` - matrix

        - ``twist`` - whether we consider the twisted homology (mainly useful for
          veering triangulations)

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)", mutable=True)
            sage: A = matrix([[1,1],[-1,0],[0,-1]])
            sage: B = copy(A)
            sage: for e in [0,1,0,1]:
            ....:     T.flip_homological_action(e, B)
            ....:     T.flip(e)
            sage: B
            [ 1 -3]
            [-1  4]
            [ 0 -1]

        In the above example we are back to the initial triangulation and
        one can recover the homology matrix using pseudo inverses::

            sage: A.pseudoinverse() * B
            [ 1 -4]
            [ 0  1]

        Another torus example (a single Dehn twist along the curve
        w)::

            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)", mutable=True)
            sage: u = [0,1,0,0,-1,0]
            sage: v = [0,0,0,1,0,-1]
            sage: w = [1,1,1,1,0,0]
            sage: A = matrix(ZZ, 3, 6, [u,v,w]).transpose()
            sage: B = copy(A)
            sage: for i in (0,2,1,3):
            ....:     T.flip_homological_action(i, B)
            ....:     T.flip(i)
            ....:     T._check_homology_matrix(B)
            sage: p = "(0,2)(~0,~2)(1,3)(~1,~3)"
            sage: T.relabel_homological_action(p, B)
            sage: T.relabel(p)
            sage: T._check_homology_matrix(B)
            sage: A.pseudoinverse() * B
            [1 0 0]
            [0 1 0]
            [1 1 1]
        """
        ne = self.num_edges()
        assert m.nrows() == ne
        ep = self._ep

        if not twist and ep[e] == e:
            return
        elif ep[e] < e:
            e = ep[e]

        a, b, c, d = self.square_about_edge(e)
        # v_e use to be v_c + v_d and becomes v_d + v_a
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

        # NOTE: We might want to use optimized row operations such as
        #   swap_rows(i, j)
        #   add_multiple_of_row(i, j, s)

        A = ep[a]
        D = ep[d]
        if twist:
            m[e] = (m[a] if a < A else m[A]) + (m[d] if d < D else m[D])
        else:
            if a == A and d == D:
                try:
                    m[e] = m.row_ambient_module().zero()
                except AttributeError:
                    m[e] = m._row_ambient_module().zero()
            elif a == A:
                m[e] = m[d] if d < D else -m[D]
            elif d == D:
                m[e] = m[a] if a < A else -m[A]
            else:
                m[e] = (m[a] if a < A else -m[A]) + (m[d] if d < D else -m[D])

    def relabel_homological_action(self, p, m, twist=False, check=True):
        r"""
        Acts on the homological vectors ``m`` by the relabeling permutation ``p``.

        The matrix ``m`` must have ``ne`` rows and each column represents a
        vector in cohomology. That is to say, each column has to satisfy the
        following conditions:

        - for each face `F` the evaluation on its boundary is zero.

        - for each folded edge, the corresponding entry is zero.

        INPUT:

        - ``p`` - automorphism

        - ``m`` - matrix

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)", mutable=True)
            sage: u = [0,1,0,0,-1,0]
            sage: v = [0,0,0,1,0,-1]
            sage: w = [1,1,1,1,0,0]
            sage: A = matrix(ZZ, 3, 6, [u,v,w]).transpose()
            sage: T._check_homology_matrix(A)
            sage: perms = ["(0,1)(~0,~1)", "(3,~3)", "(3,~5)(~3,5)",
            ....:          "(0,1,2,~0,~1,~2)",
            ....:          "(0,~1,4,~0,1,~4)(2,3)(~2,~3)"]
            sage: for p in perms * 5:
            ....:     T.relabel_homological_action(p, A)
            ....:     T.relabel(p)
            ....:     T._check_homology_matrix(A)
        """
        n = self._n
        ne = self.num_edges()
        if check and not perm_check(p, n):
            p = perm_init(p, self._n, self._ep)
            if not perm_check(p, n):
                raise ValueError('invalid relabeling permutation')

        q = perm_invert(p)

        ep = self._ep
        seen = [False] * ne

        for e0 in range(ne):
            if seen[e0]:
                continue

            seen[e0] = True

            e = q[e0]
            E = ep[e]
            if E < e:
                is_e0_neg = True
                e, E = E, e
            else:
                is_e0_neg = False
            while not seen[e]:
                assert e < ne
                seen[e] = True

                ee = q[e]
                EE = ep[ee]

                if EE < ee:
                    is_neg = True
                    ee, EE = EE, ee
                else:
                    is_neg = False
                m.swap_rows(e, ee)
                if is_neg and not twist:
                    m[e] *= -1

                e = ee

            # one more sign change?
            assert e == e0
            if is_e0_neg and not twist:
                m[e] *= -1

    def is_flippable(self, e):
        r"""
        Check whether the half-edge e is flippable.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~2,4)(~1,3,~3)")
            sage: T.is_flippable(0)
            True
            sage: T.is_flippable(1)
            True
            sage: T.is_flippable(3)
            False
            sage: T.is_flippable(4)
            True
        """
        e = int(e)
        E = self._ep[e]
        a = self._fp[e]
        b = self._fp[a]
        return a != E and b != E

    def flippable_edges(self):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.flippable_edges()
            [0, 1, 2]
            sage: V = VeeringTriangulation(T, [RED, RED, BLUE])
            sage: V.flippable_edges()
            [0, 1]
        """
        n = self._n
        ep = self._ep
        return [e for e in range(n) if e <= ep[e] and self.is_flippable(e)]

    def square_about_edge(self, e):
        r"""
        Return the four edges that makes ``e`` the diagonal of a quadrilateral.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.square_about_edge(0)
            (1, 2, 4, 3)

            sage: T = Triangulation("(0,1,2)")
            sage: T.square_about_edge(0)
            (1, 2, 1, 2)
        """
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
        E = self._ep[e]

        a = self._fp[e]
        b = self._fp[a]
        c = self._fp[E]
        d = self._fp[c]

        return a, b, c, d

    def swap(self, e):
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

    def flip(self, e, check=True):
        r"""
        Flip the edge ``e``.

        EXAMPLES::

            sage: from veerer import Triangulation

        Flipping a folded edge is an involution::

            sage: T = Triangulation("(0,1,2)", mutable=True)
            sage: T
            Triangulation("(0,1,2)")
            sage: T.flip(0)
            sage: T
            Triangulation("(0,2,1)")
            sage: T.flip(0)
            sage: T
            Triangulation("(0,1,2)")

            sage: T == Triangulation("(0,1,2)")
            True
        """
        # v<----------u     v<----------u
        # |     a    ^^     |^    a     ^
        # |         / |     | \         |
        # |  F     /  |     |  \     G  |
        # |       /   |     |   \       |
        # |b    e/   d| --> |b   \     d|
        # |     /     |     |     \     |
        # |    /      |     |     e\    |
        # |   /       |     |       \   |
        # |  /     G  |     | F      \  |
        # | /         |     |         \ |
        # v/    c     |     v     c    \|
        # w---------->x     w---------->x

        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

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

        # Disabled for now
        # F = self._fl[e]
        # G = self._fl[E]

        # v = self._vl[b]
        # x = self._vl[d]

        # fix face perm and cycles
        self._fp[e] = b
        self._fp[b] = c
        self._fp[c] = e
        self._fp[a] = E
        self._fp[E] = d
        self._fp[d] = a

        # Face labels
        # self._fl[a] = G
        # self._fl[c] = F

        # fix vertex perm
        self._vp[a] = D
        self._vp[b] = E
        self._vp[E] = A
        self._vp[c] = B
        self._vp[d] = e
        self._vp[e] = C

        # Vertex labels
        # self._vl[e] = x
        # self._vl[E] = v

    def flip_back(self, e, check=True):
        r"""
        Flip back the edge ``e``.

        EXAMPLES::

            sage: from veerer import *
            sage: T0 = Triangulation([(0,1,2),(-1,-2,-3)], mutable=True)
            sage: T = T0.copy()
            sage: T.flip(0)
            sage: T.flip_back(0)
            sage: T == T0
            True

            sage: T.flip(1)
            sage: T.flip(2)
            sage: T.flip_back(2)
            sage: T.flip_back(1)
            sage: T == T0
            True

            sage: T.set_immutable()
            sage: T.flip(0)
            Traceback (most recent call last):
            ...
            ValueError: immutable triangulation; use a mutable copy instead
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

        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

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

        # Disabled for now
        # F = self._fl[e]
        # G = self._fl[E]

        # v = self._vl[b]
        # x = self._vl[d]

        # fix face perm and cycles
        self._fp[e] = d
        self._fp[d] = a
        self._fp[a] = e
        self._fp[b] = c
        self._fp[c] = E
        self._fp[E] = b

        # Face labels
        # Disabled for now
        # self._fl[a] = G
        # self._fl[c] = F

        # fix vertex perm
        self._vp[a] = D
        self._vp[b] = e
        self._vp[e] = A
        self._vp[c] = B
        self._vp[d] = E
        self._vp[E] = C

        # Vertex labels
        # Disabled for now
        # self._vl[e] = x
        # self._vl[E] = v

    def to_string(self):
        r"""
        Serialize this triangulation as a string.

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.to_string()
            '6_120534_543210'
        """
        return uint_base64_str(self._n) + '_' + perm_base64_str(self._fp) + '_' + perm_base64_str(self._ep)

    def conjugate(self):
        r"""
        Conjugate this triangulation.

        The face permutation is replaced by its inverse, conjugated by
        the edge permutation. The vertex permutation is replaced by
        its inverse.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)", mutable=True)
            sage: T.conjugate()
            sage: T
            Triangulation("(0,2,4)(1,3,5)(~5,~4,~3)(~2,~1,~0)")

            sage: T = Triangulation("(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)", mutable=False)
            sage: T.conjugate()
            Traceback (most recent call last):
            ...
            ValueError: immutable triangulation; use a mutable copy instead
        """
        # for reference
        #
        # x<----------x     x
        # |     a    ^^     ^\
        # |         /       | \
        # |        /        |  \
        # |       /         |   \
        # |b    c/      --> |    \
        # |     /           |     \
        # |    /            |B    C\
        # |   /             |       \
        # |  /              |        \
        # | /               |         \
        # v/                |     A    v
        # x                 x<----------x
        #
        # (a, b, c)     -->  (C, B, A)

        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        self._fp = perm_conjugate(perm_invert(self._fp), self._ep)
        self._vp = perm_invert(self._vp)

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
            s = ((lv[i] * n + le[i]) * n + lf[i]) * n + lvef[i]
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
        """
        fp = self._fp
        ep = self._ep

        best = None
        best_relabellings = []

        for start_edge in self._automorphism_good_starts():
            relabelling = self._relabelling_from(start_edge)

            fp_new = perm_conjugate(fp, relabelling)
            ep_new = perm_conjugate(ep, relabelling)

            T = (fp_new, ep_new)
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

            T = (fp_new, ep_new)
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
            '9_345876210_087654321'
            sage: TT = Triangulation.from_string(T.iso_sig())
            sage: TT
            Triangulation("(0,3,~1)(1,4,~2)(2,~4,~3)")
            sage: TT.iso_sig() == T.iso_sig()
            True

            sage: T = Triangulation("(0,10,~6)(1,12,~2)(2,14,~3)(3,16,~4)(4,~13,~5)(5,~1,~0)(6,~17,~7)(7,~14,~8)(8,13,~9)(9,~11,~10)(11,~15,~12)(15,17,~16)")
            sage: T.iso_sig()
            'A_lkp9ijfmcvegq0osnyutxdrbaw87654321zh_zyxwvutsrqponmlkjihgfedcba9876543210'
            sage: Triangulation.from_string(T.iso_sig())
            Triangulation("(0,~14,13)(1,~15,~2)(2,~10,~3)(3,9,~4)(4,~17,~5)(5,~16,~6)(6,15,~7)(7,~13,~8)(8,12,~9)(10,14,~11)(11,16,~12)(17,~1,~0)")

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = VeeringTriangulation(t, cols, mutable=True)
            sage: T.iso_sig()
            'RBRBRBBBRBBBBBBRBBBRBRBR_7ad9e8fbhjl0mkig654321nc_nmlkjihgfedcba9876543210'

        If we relabel the triangulation, the isomorphic signature does not change::

            sage: from veerer.permutation import perm_random
            sage: p = perm_random(24)
            sage: T.relabel(p)
            sage: T.iso_sig()
            'RBRBRBBBRBBBBBBRBBBRBRBR_7ad9e8fbhjl0mkig654321nc_nmlkjihgfedcba9876543210'

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
            'GBRGRGBRB_731264580_082375641'
        """
        T = self.copy(mutable=True)
        T.set_canonical_labels()
        return T.to_string()

    def _non_isom_easy(self, other):
        return self._n != other._n or self.num_folded_edges() != other.num_folded_edges()

    def is_isomorphic_to(self, other, certificate=False):
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
            ....:     assert T.is_isomorphic_to(TT)

            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: V = VeeringTriangulation(fp, cols, mutable=True)
            sage: W = V.copy()
            sage: p = perm_random_centralizer(V.edge_permutation(copy=False))
            sage: W.relabel(p)
            sage: assert V.is_isomorphic_to(W) is True
            sage: ans, cert = V.is_isomorphic_to(W, True)
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

    def cover(self, c, mutable=False, check=True):
        from .cover import TriangulationCover
        return TriangulationCover(self, c, mutable=mutable, check=check)

    def _check_xy(self, x, y):
        if len(x) == self.num_edges():
            x = [x[self._norm(i)] for i in range(self._n)]
        elif len(x) != self._n or any(a < 0 for a in x):
            raise ValueError('invalid argument x')
        if len(y) == self.num_edges():
            y = [y[self._norm(i)] for i in range(self._n)]
        elif len(y) != self._n or any(a < 0 for a in y):
            raise ValueError('invalid argument y')
        return (x, y)

    def colouring_from_xy(self, x, y, check=True):
        r"""
        Return the veering colouring associated with the holonomy data ``x`` and ``y``.

        EXAMPLES::

            sage: from veerer import Triangulation
            sage: t = "(0,1,2)(~0,~1,~2)"
            sage: t = Triangulation(t)
            sage: x = [1, 2, 1]
            sage: y = [1, 1, 2]
            sage: t.colouring_from_xy(x, y)
            array('i', [1, 2, 2, 2, 2, 1])
            sage: t = "(0,1,2)(~0,~1,4)(~2,5,3)(~3,~4,~5)"
            sage: t = Triangulation(t)
            sage: x = [1, 2, 1, 2, 1, 1]
            sage: y = [1, 1, 2, 1, 2, 3]
            sage: t.colouring_from_xy(x, y)
            array('i', [1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1])
        """
        from .constants import BLUE, RED, PURPLE, GREEN

        if check:
            x, y = self._check_xy(x, y)

        def check_set_colouring(colouring, e, ep, col):
            if colouring[e] is None:
                colouring[e] = colouring[ep[e]] = col
            elif colouring[e] != col:
                raise ValueError('inconsistent colouring between x and y for edge e={}'.format(e))

        colouring = [None] * self._n
        faces = perm_cycles(self._fp, True, self._n)
        ep = self._ep
        for face in faces:
            a, b, c = face
            x_degenerate = (x[a] == 0) + (x[b] == 0) + (x[c] == 0)
            if x_degenerate == 0:
                if x[a] == x[b] + x[c]:
                    xlarge = 0
                elif x[b] == x[c] + x[a]:
                    xlarge = 1
                elif x[c] == x[a] + x[b]:
                    xlarge = 2
                else:
                    raise ValueError('inconsistent x data for triangle {}'.format(face))
                check_set_colouring(colouring, face[(xlarge + 1) % 3], ep, BLUE)
                check_set_colouring(colouring, face[(xlarge + 2) % 3], ep, RED)
            elif x_degenerate == 1:
                if x[a] == 0:
                    xvert = 0
                elif x[b] == 0:
                    xvert = 1
                elif x[c] == 0:
                    xvert = 2
                check_set_colouring(colouring, face[xvert], ep, GREEN)
                check_set_colouring(colouring, face[(xvert + 1) % 3], ep, RED)
                check_set_colouring(colouring, face[(xvert + 2) % 3], ep, BLUE)
            else:
                raise ValueError('inconsistent x data for triangle {}'.format(face))

            y_degenerate = (y[a] == 0) + (y[b] == 0) + (y[c] == 0)
            if y_degenerate == 0:
                if y[a] == y[b] + y[c]:
                    ylarge = 0
                elif y[b] == y[c] + y[a]:
                    ylarge = 1
                elif y[c] == y[a] + y[b]:
                    ylarge = 2
                else:
                    raise ValueError('inconsistent y data for triangle {}'.format(face))
                check_set_colouring(colouring, face[(ylarge + 1) % 3], ep, RED)
                check_set_colouring(colouring, face[(ylarge + 2) % 3], ep, BLUE)
            elif y_degenerate == 1:
                if y[a] == 0:
                    yhor = 0
                elif y[b] == 0:
                    yhor = 1
                elif y[c] == 0:
                    yhor = 2
                check_set_colouring(colouring, face[yhor], ep, PURPLE)
                check_set_colouring(colouring, face[(yhor + 1) % 3], ep, BLUE)
                check_set_colouring(colouring, face[(yhor + 2) % 3], ep, RED)
            else:
                raise ValueError('inconsistent y data for triangle {}'.format(face))

        return array('i', colouring)
