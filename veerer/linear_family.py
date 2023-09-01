r"""
Linear family of coordinates on a veering triangulation
"""
######################################################################
# This file is part of veering.
#
#       Copyright (C) 2023 Vincent Delecroix
#
# veerer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# veerer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with veerer. If not, see <https://www.gnu.org/licenses/>.
######################################################################

from array import array
from copy import copy
import collections
import itertools
import numbers
from random import choice, shuffle

from .env import sage, ppl
from .constants import VERTICAL, HORIZONTAL, BLUE, RED
from .permutation import perm_cycle_string, perm_cycles, perm_check, perm_conjugate, perm_on_list
from .polyhedron import LinearExpressions, ConstraintSystem
from .veering_triangulation import VeeringTriangulation


if sage is not None:
    from sage.structure.element import get_coercion_model
    from sage.rings.integer_ring import ZZ
    from sage.rings.rational_field import QQ
    from sage.matrix.constructor import matrix
    from sage.geometry.polyhedron.constructor import Polyhedron
    from sage.arith.misc import gcd
    from sage.categories.number_fields import NumberFields

    cm = get_coercion_model()
    _NumberFields = NumberFields()
else:
    matrix = None
    Polyhedron = None
    gcd = None
    ZZ = None
    QQ = None

    cm = None
    _NumberFields = None


def subspace_are_equal(subspace1, subspace2, check=True):
    r"""
    Test whether the subspaces generated by the rows of ``subspace1`` and
    ``subspace2`` are equal.

    INPUT:

    - ``subspace1``, ``subspace2`` -- full rank matrices with the same number
      of columns

    - ``check`` -- boolean (default ``True``)

    EXAMPLES::

        sage: from veerer.linear_family import subspace_are_equal

        sage: m1 = random_matrix(ZZ, 3, 5)
        sage: m2 = copy(m1)
        sage: m2.add_multiple_of_row(0, 1, -2)
        sage: m2.add_multiple_of_row(0, 2, 1)
        sage: m2.add_multiple_of_row(1, 2, 1)
        sage: subspace_are_equal(m1, m2)
        True

        sage: m1 = matrix(ZZ, [[1, 1, 0]])
        sage: m2 = matrix(ZZ, [[1, 1, -1]])
        sage: subspace_are_equal(m1, m2)
        False
    """
    if check:
        if subspace1.ncols() != subspace2.ncols():
            raise ValueError('subspace1 and subspace2 of different ambient dimensions')
        if subspace1.rank() != subspace1.nrows():
            raise ValueError('subspace1 not full rank')
        if subspace2.rank() != subspace2.nrows():
            raise ValueErrror('subspace2 not full rank')

    n = subspace1.nrows()
    if n != subspace2.nrows():
        return False

    base_ring = cm.common_parent(subspace1.base_ring(), subspace2.base_ring())
    mat = matrix(base_ring, n + 1, subspace1.ncols())
    mat[:n] = subspace1
    for v in subspace2.rows():
        mat[n] = v
        r = mat.rank()
        if r < n:
            raise RuntimeError('matrices where expected to be full rank')
        if r > n:
            return False
    return True


def subspace_cmp(subspace1, subspace2, check=True):
    if check:
        if subspace1.ncols() != subspace2.ncols():
            raise ValueError('subspace1 and subspace2 of different ambient dimensions')
        if subspace1.rank() != subspace1.nrows():
            raise ValueError('subspace1 not full rank')
        if subspace2.rank() != subspace2.nrows():
            raise ValueErrror('subspace2 not full rank')

    n = subspace1.nrows()
    if n != subspace2.nrows():
        return False

    base_ring = cm.common_parent(subspace1.base_ring(), subspace2.base_ring())
    subspace1 = subspace1.echelon_form()
    subspace2 = subspace2.echelon_form()
    for r1, r2 in zip(subspace1, subspace2):
        c = (r1 > r2) - (r1 < r2)
        if c:
            return c
    return 0


def relabel_on_edges(ep, r, n, m):
    r"""
    INPUT:

    - ep - edge permutation

    - r - relabelling permutation on half edges (list of length n)

    - n - num half edges

    - m - num edges

    OUTPUT: list of length m

    EXAMPLES::

        sage: from array import array
        sage: from veerer.linear_family import relabel_on_edges

        sage: ep = array('i', [8, 1, 2, 7, 4, 5, 6, 3, 0])
        sage: r = array('i', [3, 0, 5, 4, 6, 2, 1, 8, 7])
        sage: relabel_on_edges(ep, r, 9, 7)
        array('i', [3, 0, 5, 4, 6, 2, 1])
    """
    rr = array('i', [-1] * m)
    for i in range(m):
        if ep[i] < i:
            raise ValueError("not in canonical form")
        j = r[i]
        k = r[ep[i]]
        if (j >= m and k >= m):
            raise ValueError("relabelling not preserving canonical form")
        if j < k:
            rr[i] = j
        else:
            rr[i] = k
    return rr


def matrix_permutation(mat, perm):
    r"""
    Permute in place the columnns of ``mat`` according to the permutation ``perm``

    EXAMPLES::

        sage: from array import array
        sage: from veerer.linear_family import matrix_permutation

        sage: p = array('i', [2, 0, 1])
        sage: m = matrix([[0, 1, 2]])
        sage: matrix_permutation(m, p)
        sage: print(m)
        [1 2 0]
    """
    m = mat.ncols()
    for c in perm_cycles(perm, False, m):
        for i in range(1, len(c)):
            mat.swap_columns(c[0], c[i])


class VeeringTriangulationLinearFamily(VeeringTriangulation):
    r"""
    Veering triangulation together with a subspace of H^1(S, Sigma; \bR) that
    describes a (piece of a) linear GL(2,R)-invariant immersed sub-orbifold.
    """
    __slots__ = ['_subspace']

    def __init__(self, *args, mutable=False, check=True):
        if len(args) == 2:
            vt, subspace = args
            t = vt
            colouring = vt._colouring
        elif len(args) == 3:
            t, colouring, subspace = args
        VeeringTriangulation.__init__(self, t, colouring, mutable=mutable, check=False)

        if not isinstance(subspace, sage.structure.element.Matrix):
            subspace = matrix(subspace)

        self._subspace = subspace
        self._subspace.echelonize()
        if not mutable:
            self._subspace.set_immutable()

        if check:
            self._check(ValueError)

    def __getstate__(self):
        R = self.base_ring()
        a = list(self._fp[:])
        a.extend(self._ep)
        a.extend(self._colouring)
        a.append(self._mutable)
        # NOTE: here we know that all entries have the same parent
        if R is ZZ:
            entries = list(map(int, self._subspace.list()))
        elif R is QQ:
            entries = []
            for x in self._subspace.list():
                entries.append(int(x.numerator()))
                entries.append(int(x.denominator()))
        elif _NumberFields is not None and R in _NumberFields:
            entries = []
            for x in self._subspace.list():
                v = x.vector()
                d = v.denominator()
                entries.extend(int(x) for x in (v * d))
                entries.append(int(d))
        else:
            entries = self._subspace.list()
        return a, R, entries

    def __setstate__(self, arg):
        r"""
        TESTS:

        An example over ZZ::

            sage: from veerer import VeeringTriangulation
            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: lf = vt.as_linear_family()
            sage: lf2 = loads(dumps(lf))  # indirect doctest
            sage: assert lf == lf2
            sage: lf2._check()

        An example over number fields::

            sage: from veerer.linear_family import VeeringTriangulationLinearFamilies
            sage: lf = VeeringTriangulationLinearFamilies.triangle_3_4_13_unfolding_orbit_closure()
            sage: lf2 = loads(dumps(lf))  # indirect doctest
            sage: assert lf == lf2
            sage: lf2._check()
        """
        a, R, raw_entries = arg
        n = (len(a) - 1) // 3
        self._n = n
        self._fp = array('i', a[:n])
        self._ep = array('i', a[n:2*n])
        self._colouring = array('i', a[2*n:3*n])
        self._mutable = a[-1]
        self._vp = array('i', [-1] * n)
        for i in range(n):
            self._vp[self._fp[self._ep[i]]] = i

        if R is ZZ:
            entries = [ZZ(x) for x in raw_entries]
        elif R is QQ:
            entries = []
            for i in range(0, len(raw_entries), 2):
                entries.append(QQ((raw_entries[i], raw_entries[i + 1])))
        elif _NumberFields is not None and R in _NumberFields:
            entries = []
            d = R.degree()
            for i in range(0, len(raw_entries), d + 1):
                v = tuple(x / raw_entries[i + d] for x in raw_entries[i:i + d])
                entries.append(R(v))
        self._subspace = matrix(R, len(entries) // self.num_edges(), self.num_edges(), entries)
        if not self._mutable:
            self._subspace.set_immutable()

    def veering_triangulation(self, mutable=False):
        r"""
        Return the underlying veering triangulation.

        EXAMPLES::

            sage: from veerer import VeeringTriangulation
            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: vt.as_linear_family().veering_triangulation()
            VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB")
        """
        return VeeringTriangulation.copy(self, mutable)

    def _horizontal_subspace(self):
        mat = copy(self._subspace)
        ne = self.num_edges()
        ep = self._ep
        for j in range(ne):
            if ep[j] < j:
                raise ValueError('not in standard form')
            if self._colouring[j] == BLUE:
                for i in range(mat.nrows()):
                    mat[i, j] *= -1
        return mat

    def as_linear_family(self):
        return self

    def conjugate(self):
        raise NotImplementedError

    def rotate(self):
        r"""
        Conjugate this family.

        EXAMPLES::

            sage: from veerer import *

            sage: fp = "(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)"
            sage: cols = "BRRBRR"
            sage: f = VeeringTriangulation(fp, cols).as_linear_family(mutable=True)
            sage: f.rotate()
            sage: f
            VeeringTriangulationLinearFamily("(0,1,2)(3,4,5)(~5,~3,~1)(~4,~2,~0)", "RBBRBB", [(1, 0, -1, 0, 0, 0), (0, 1, 1, 0, 1, 1), (0, 0, 0, 1, 0, -1)])

            sage: fp = "(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)"
            sage: cols = "RRRRRRBBBBBBBBBBBB"
            sage: f = VeeringTriangulation(fp, cols).as_linear_family(mutable=True)
            sage: f.rotate()
            sage: f
            VeeringTriangulationLinearFamily("(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)", "BBBBBBRRRRRRRRRRRR", [(1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 2, 2, 1, 1, 1, 0, 0), (0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0), (0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1), (0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0), (0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)])
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation family; use a mutable copy instead')

        subspace = self._horizontal_subspace()
        subspace.echelonize()
        VeeringTriangulation.rotate(self)
        self._subspace = subspace

        # TODO: remove check
        self._check()

    def _set_subspace_constraints(self, insert, x, slope):
        ambient_dim = self._subspace.ncols()
        if slope == VERTICAL:
            subspace = self._subspace
        elif slope == HORIZONTAL:
            subspace = self._horizontal_subspace()
        for row in subspace.right_kernel_matrix():
            insert(sum(row[i] * x[i] for i in range(ambient_dim)) == 0)

    def copy(self, mutable=None):
        r"""
        Return a copy of this linear family.
        """
        if mutable is None:
            mutable = self._mutable

        if not self._mutable and not mutable:
            # avoid copies of immutable objects
            return self

        L = VeeringTriangulationLinearFamily.__new__(VeeringTriangulationLinearFamily)
        L._n = self._n
        L._vp = self._vp[:]
        L._ep = self._ep[:]
        L._fp = self._fp[:]
        L._colouring = self._colouring[:]
        L._subspace = copy(self._subspace)
        L._mutable = mutable
        if not mutable:
            L._subspace.set_immutable()
        return L

    def base_ring(self):
        return self._subspace.base_ring()

    def set_immutable(self):
        VeeringTriangulation.set_immutable(self)
        self._subspace.set_immutable()

    def __hash__(self):
        r"""
        TESTS::

            sage: from veerer import *
            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: h = hash(vt.as_linear_family())
        """
        if self._mutable:
            raise ValueError('mutable veering triangulation linear family not hashable')

        x = 140737488617563
        x = ((x ^ hash(self._vp.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._ep.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._fp.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._colouring.tobytes())) * 2147483693) + 82520 + self._n + self._n
        x = ((x ^ hash(self._subspace) * 2147483693)) + 82520 + self._n + self._n

        return x

    def __str__(self):
        r"""
        Return a string representation.

        TESTS::

            sage: from veerer import *
            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: str(vt.as_linear_family())
            'VeeringTriangulationLinearFamily("(0,1,2)(~2,~0,~1)", "RRB", [(1, 0, -1), (0, 1, 1)])'
        """
        return "VeeringTriangulationLinearFamily(\"{}\", \"{}\", {})".format(
               perm_cycle_string(self._fp, False, self._n, self._ep),
               self._colouring_string(short=True),
               self._subspace.rows())

    def __repr__(self):
        return str(self)

    def _check(self, error=ValueError):
        subspace = self._subspace
        VeeringTriangulation._check(self, error)
        if subspace.ncols() != self.num_edges():
            raise error('subspace matrix has wrong dimension')
        if subspace.rank() != subspace.nrows():
            raise error('subspace matrix is not of full rank')
        # test that elements satisfy the switch condition
        for v in subspace.rows():
            self._set_switch_conditions(self._tt_check, v, VERTICAL)
        if subspace != subspace.echelon_form():
            raise error('subspace not in echelon form')
        if self._mutable != self._subspace.is_mutable():
            raise error('incoherent mutability states')

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from veerer import *
            sage: vt, s, t = VeeringTriangulations.L_shaped_surface(1,1,1,1)
            sage: s = vector(QQ, s)
            sage: t = vector(QQ, t)
            sage: f1 = VeeringTriangulationLinearFamily(vt, [s, t])
            sage: f2 = VeeringTriangulationLinearFamily(vt, [s + 2*t, -s - t])
            sage: f1 == f2
            True
            sage: from veerer import *
            sage: vt2, s2, t2 = VeeringTriangulations.L_shaped_surface(1,2,1,3)
            sage: f3 = VeeringTriangulationLinearFamily(vt2, [s2, t2])
            sage: f1 == f3
            False
            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: vt.as_linear_family() == f1
            False
        """
        if type(self) is not type(other):
            raise TypeError
        return VeeringTriangulation.__eq__(self, other) and self._subspace == other._subspace

    def __ne__(self, other):
        r"""
        TESTS::

            sage: from veerer import *
            sage: vt, s, t = VeeringTriangulations.L_shaped_surface(1,1,1,1)
            sage: s = vector(QQ, s)
            sage: t = vector(QQ, t)
            sage: f1 = VeeringTriangulationLinearFamily(vt, [s, t])
            sage: f2 = VeeringTriangulationLinearFamily(vt, [s + 2*t, -s - t])
            sage: f1 != f2
            False
            sage: from veerer import *
            sage: vt2, s2, t2 = VeeringTriangulations.L_shaped_surface(1,2,1,3)
            sage: f3 = VeeringTriangulationLinearFamily(vt2, [s2, t2])
            sage: f1 != f3
            True
            sage: vt = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: vt.as_linear_family() != f1
            True
        """
        if type(self) is not type(other):
            raise TypeError
        return VeeringTriangulation.__ne__(self, other) or self._subspace != other._subspace

    def _richcmp_(self, other, op):
        c = (self._n > other._n) - (self._n < other._n)
        if c:
            return rich_to_bool(op, c)

        c = (self._colouring > other._colouring) - (self._colouring < other._colouring)
        if c:
            return rich_to_bool(op, c)

        c = (self._fp > other._fp) - (self._fp < other._fp)
        if c:
            return rich_to_bool(op, c)

        c = (self._ep > other._ep) - (self._ep < other._ep)
        if c:
            return rich_to_bool(op, c)

        c = subspace_cmp(self._subspace, other._subspace)
        return rich_to_bool(op, c)

    def train_track_polytope(self, slope=VERTICAL, low_bound=0, backend=None):
        r"""
        Return the polytope of non-negative elements in the subspace.

        EXAMPLES::

            sage: from veerer import *
            sage: vt, s, t = VeeringTriangulations.L_shaped_surface(1, 3, 1, 1)
            sage: f = VeeringTriangulationLinearFamily(vt, [s, t])
            sage: f.train_track_polytope(VERTICAL)
            Cone of dimension 2 in ambient dimension 7 made of 2 facets (backend=ppl)
            sage: f.train_track_polytope(HORIZONTAL)
            Cone of dimension 2 in ambient dimension 7 made of 2 facets (backend=ppl)

            sage: sorted(f.train_track_polytope(VERTICAL).rays())
            [[0, 1, 3, 3, 1, 1, 0], [1, 0, 0, 1, 1, 1, 1]]
            sage: sorted(f.train_track_polytope(HORIZONTAL).rays())
            [[1, 0, 0, 1, 1, 1, 1], [3, 1, 3, 0, 2, 2, 3]]
        """
        ne = self.num_edges()
        L = LinearExpressions(self.base_ring())
        cs = ConstraintSystem()
        for i in range(ne):
            cs.insert(L.variable(i) >= low_bound)
        self._set_subspace_constraints(cs.insert, [L.variable(i) for i in range(ne)], slope)
        return cs.cone(backend)

    def dimension(self):
        r"""
        Return the dimension of the linear family.
        """
        return self._subspace.nrows()

    def relabel(self, p, check=True):
        r"""
        Relabel inplace the veering triangulation linear family according to the permutation ``p``.

        Relabelling the subspace as well::

            sage: from veerer import VeeringTriangulations, VeeringTriangulationLinearFamily
            sage: from veerer.permutation import perm_random_centralizer

            sage: vt, s, t = VeeringTriangulations.L_shaped_surface(2, 3, 4, 5, 1, 2)
            sage: f = VeeringTriangulationLinearFamily(vt, [s, t], mutable=True)
            sage: for _ in range(10):
            ....:     p = f._relabelling_from(choice(range(9)))
            ....:     f.relabel(p)
            ....:     f._check()

            sage: for _ in range(10):
            ....:     p = perm_random_centralizer(f.edge_permutation(copy=False))
            ....:     f.relabel(p)
            ....:     f._check()
        """
        n = self.num_half_edges()
        m = self.num_edges()
        ep = self._ep
        if check and not perm_check(p, n):
            p = perm_init(p, n, ep)
            if not perm_check(p, n):
                raise ValueError('invalid relabelling permutation')

        rr = relabel_on_edges(self._ep, p, n, m)
        matrix_permutation(self._subspace, rr)
        self._subspace.echelonize()
        VeeringTriangulation.relabel(self, p, False)

        # TODO: remove check
        self._check()

    def best_relabelling(self, all=False):
        n = self.num_half_edges()
        m = self.num_edges()

        best = None
        if all:
            relabellings = []

        for start_edge in self._automorphism_good_starts():
            relabelling = self._relabelling_from(start_edge)
            rr = relabel_on_edges(self._ep, relabelling, n, m)

            fp = perm_conjugate(self._fp, relabelling)
            ep = perm_conjugate(self._ep, relabelling)
            cols = self._colouring[:]
            perm_on_list(relabelling, cols)
            subspace = copy(self._subspace)
            matrix_permutation(subspace, rr)
            subspace.echelonize()
            subspace.set_immutable()

            T = (cols, fp, ep, subspace)
            if best is None or T < best:
                best = T
                best_relabelling = relabelling
                if all:
                    del relabellings[:]
                    relabellings.append(relabelling)
            elif all and T == best:
                relabellings.append(relabelling)

        return (relabellings, best) if all else (best_relabelling, best)

    # TODO: change to canonicalize ? Since we also need to canonicalize the subspace
    # it is not only about labels
    def _non_isom_easy(self, other):
        return (VeeringTriangulation._non_isom_easy(self, other) or
                self._subspace.nrows() != other._subspace.nrows())

    def flip(self, e, col, check=True):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)

            sage: L = VeeringTriangulationLinearFamily(T, [s, t], mutable=True)
            sage: T = T.copy(mutable=True)

            sage: T.flip(3, 2)
            sage: L.flip(3, 2)
            sage: T
            VeeringTriangulation("(0,3,2)(1,4,~0)(5,6,~1)", "BRRBBBB")
            sage: L
            VeeringTriangulationLinearFamily("(0,3,2)(1,4,~0)(5,6,~1)", "BRRBBBB", [(1, 0, 0, 1, 1, 1, 1), (0, 1, 1, -1, 1, 1, 0)])

            sage: L.flip(4, 2)
            sage: T.flip(4, 2)
            sage: T
            VeeringTriangulation("(0,3,2)(1,~0,4)(5,6,~1)", "BRRBBBB")
            sage: L
            VeeringTriangulationLinearFamily("(0,3,2)(1,~0,4)(5,6,~1)", "BRRBBBB", [(1, 0, 0, 1, 1, 1, 1), (0, 1, 1, -1, -1, 1, 0)])

            sage: T.flip(5, 2)
            sage: L.flip(5, 2)
            sage: T
            VeeringTriangulation("(0,3,2)(1,~0,4)(5,~1,6)", "BRRBBBB")
            sage: L
            VeeringTriangulationLinearFamily("(0,3,2)(1,~0,4)(5,~1,6)", "BRRBBBB", [(1, 0, 0, 1, 1, 1, 1), (0, 1, 1, -1, -1, -1, 0)])
        """
        VeeringTriangulation.flip(self, e, col, Gx=self._subspace, check=check)
        self._subspace.echelonize()

    def flip_back(self, e, col, check=True):
        VeeringTriangulation.flip_back(self, e, col, Gx=self._subspace, check=check)
        self._subspace.echelonize()

    def geometric_polytope(self, x_low_bound=0, y_low_bound=0, hw_bound=0, backend=None):
        r"""
        Return the geometric polytope.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
            sage: T.geometric_polytope()
            Cone of dimension 4 in ambient dimension 6 made of 6 facets (backend=ppl)
            sage: T.as_linear_family().geometric_polytope(backend='ppl')
            Cone of dimension 4 in ambient dimension 6 made of 6 facets (backend=ppl)
            sage: T.as_linear_family().geometric_polytope(backend='sage')
            Cone of dimension 4 in ambient dimension 6 made of 6 facets (backend=sage)

        An example in genus 2 involving a linear constraint::

            sage: vt, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: f = VeeringTriangulationLinearFamily(vt, [s, t])
            sage: PG = f.geometric_polytope(backend='ppl')
            sage: PG
            Cone of dimension 4 in ambient dimension 14 made of 6 facets (backend=ppl)
            sage: sorted(PG.rays())
            [[0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1],
             [0, 1, 1, 1, 1, 1, 0, 2, 0, 0, 2, 2, 2, 2],
             [0, 1, 1, 1, 1, 1, 0, 2, 2, 2, 0, 0, 0, 2],
             [0, 2, 2, 2, 2, 2, 0, 1, 1, 1, 0, 0, 0, 1],
             [1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
             [1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1],
             [2, 0, 0, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, 1]]
        """
        ne = self._subspace.ncols()
        L = LinearExpressions(self.base_ring())
        x = [L.variable(i) for i in range(ne)]
        y = [L.variable(ne + i) for i in range(ne)]
        cs = ConstraintSystem()
        for i in range(ne):
            cs.insert(x[i] >= x_low_bound)
        for i in range(ne):
            cs.insert(y[i] >= y_low_bound)
        self._set_subspace_constraints(cs.insert, x, VERTICAL)
        self._set_subspace_constraints(cs.insert, y, HORIZONTAL)
        self._set_geometric_constraints(cs.insert, x, y, hw_bound=hw_bound)
        return cs.cone(backend)

    def geometric_automaton(self, run=True, backend=None):
        from .automaton import GeometricAutomatonSubspace
        A = GeometricAutomatonSubspace(backend=backend)
        A.set_seed(self)
        if run:
            A.run()
        return A

    def random_forward_flip(self, repeat=1):
        r"""
        Apply a random forward flip randomly among the ones that keeps the family core.

        INPUT:

        - ``repeat`` - integer (default 1) - if provided, make ``repeat`` flips instead
          of 1

        EXAMPLES::

            sage: from veerer.linear_family import VeeringTriangulationLinearFamilies
            sage: X = VeeringTriangulationLinearFamilies.prototype_H1_1(0, 4, 1, 0, mutable=True)
            sage: X.is_core()
            True
            sage: X.is_geometric()
            False
            sage: while not X.is_geometric():
            ....:     X.random_forward_flip()
        """
        if not self._mutable:
            raise ValueError('immutable veering triangulation family; use a mutable copy instead')

        cols = [RED, BLUE]
        for _ in range(repeat):
            e = choice(self.forward_flippable_edges())
            old_col = self._colouring[e]
            shuffle(cols)
            for c in cols:
                self.flip(e, c)
                # NOTE: for veering triangulation, the following check is made much faster
                # with edge_has_curve. However for families, we have additional linear
                # constraints on the curve that makes it harder to build.
                if self.is_core():
                    break
                else:
                    self.flip_back(e, old_col)

class VeeringTriangulationLinearFamilies:
    r"""
    A collection of linear families.
    """
    @staticmethod
    def H2_prototype_parameters(D, spin=None):
        r"""
        Return the list of prototypes in a given discriminant ``D``.

        If ``spin`` is provided, then return only the prototypes with the given
        spin parameter.

        EXAMPLES::

            sage: from veerer.linear_family import VeeringTriangulationLinearFamilies
            sage: list(VeeringTriangulationLinearFamilies.H2_prototype_parameters(5))
            [(0, 1, 1, -1)]
            sage: list(VeeringTriangulationLinearFamilies.H2_prototype_parameters(8))
            [(0, 2, 1, 0), (0, 2, 1, 0), (0, 1, 1, -2)]
            sage: list(VeeringTriangulationLinearFamilies.H2_prototype_parameters(9))
            [(0, 2, 1, -1)]
            sage: list(VeeringTriangulationLinearFamilies.H2_prototype_parameters(12))
            [(0, 3, 1, 0), (0, 3, 1, 0), (0, 1, 2, -2), (0, 2, 1, -2)]
            sage: list(VeeringTriangulationLinearFamilies.H2_prototype_parameters(13))
            [(0, 3, 1, 1), (0, 3, 1, -1), (0, 1, 1, -3)]
            sage: list(VeeringTriangulationLinearFamilies.H2_prototype_parameters(16))
            [(0, 4, 1, 0), (0, 4, 1, 0), (0, 3, 1, -2)]

            sage: list(VeeringTriangulationLinearFamilies.H2_prototype_parameters(17, spin=0))
            [(1, 2, 2, -1), (0, 4, 1, 1), (0, 2, 1, -3)]
            sage: list(VeeringTriangulationLinearFamilies.H2_prototype_parameters(17, spin=1))
            [(0, 2, 2, -1), (0, 4, 1, -1), (0, 1, 2, -3)]
        """
        # NOTE: for the spin computation, decompose D = E f^2 and the spin is
        #  (e - f) / 2 + (c + 1)(a + b + a*b) mod 2

        from sage.rings.integer_ring import ZZ

        if not isinstance(D, numbers.Integral):
            raise ValueError('D must be an integer')
        D = ZZ(D)
        if D <= 4 or not D % 4 in [0, 1]:
            raise ValueError('D must be an integer > 4 congruent to either 0 or 1 mod 4')

        if spin is not None:
            if D % 8 != 1:
                raise ValueError('spin can only be specified if D is congruent to 1 mod 8')
            if not isinstance(spin, numbers.Integral):
                raise ValueError('spin must be an integer')
            spin = ZZ(spin)
            if spin not in [0, 1]:
                raise ValueError('spin must either be 0 or 1')

        E = D.squarefree_part()
        f = (D // E).sqrtrem()[0]

        e = D % 4
        while e**2 < D:
            assert (D - e**2) % 4 == 0
            bc = (D - e**2) // 4
            for b in bc.divisors():
                c = bc // b
                # need c + e < b
                for s in (1,-1):
                    if c + s*e < b:
                        for a in range(gcd(b, c)):
                            if gcd([a, b, c, e]) == 1:
                                if spin is None or ((s*e - f) // 2 + (c + 1) * (a + b + a*b)) % 2 == spin:
                                    yield (a, b, c, s*e)
            e += 2

    @staticmethod
    def prototype_H2(a, b, c, e, mutable=False):
        r"""
        Return the Mcmullen prototype with parameters ``(a, b, c, e)`` or rather
        its quotient in `Q(1, -1^5)`.

        EXAMPLES::

            sage: from veerer.linear_family import VeeringTriangulationLinearFamilies
            sage: from veerer.automaton import GeometricAutomatonSubspace

            sage: X9 = VeeringTriangulationLinearFamilies.prototype_H2(0, 2, 1, -1)
            sage: X9.base_ring()
            Rational Field
            sage: X9.is_geometric()
            True
            sage: GeometricAutomatonSubspace(X9)
            Geometric veering linear constraint automaton with 6 vertices

            sage: X17 = VeeringTriangulationLinearFamilies.prototype_H2(0, 2, 2, -1)
            sage: X17.base_ring()
            Number Field in sqrt17 with defining polynomial x^2 - 17 with sqrt17 = 4.123105625617660?
            sage: X17.is_geometric()
            True
            sage: GeometricAutomatonSubspace(X17)  # long time
            Geometric veering linear constraint automaton with 210 vertices

        We check below part of McMullen theorem about connectedness::

            sage: a0, b0, c0, e0 = next(VeeringTriangulationLinearFamilies.H2_prototype_parameters(17, spin=0))
            sage: X17_0 = VeeringTriangulationLinearFamilies.prototype_H2(a0, b0, c0, e0)
            sage: A0 = GeometricAutomatonSubspace(X17_0, backend='sage')  # long time
            sage: for a, b, c, e in VeeringTriangulationLinearFamilies.H2_prototype_parameters(17, spin=0):  # long time
            ....:     X = VeeringTriangulationLinearFamilies.prototype_H2(a, b, c, e, mutable=True)
            ....:     X.set_canonical_labels()
            ....:     assert X in A0
            sage: a1, b1, c1, e1 = next(VeeringTriangulationLinearFamilies.H2_prototype_parameters(17, spin=1))
            sage: X17_1 = VeeringTriangulationLinearFamilies.prototype_H2(a1, b1, c1, e1)
            sage: A1 = GeometricAutomatonSubspace(X17_1, backend='sage')  # long time
            sage: for a, b, c, e in VeeringTriangulationLinearFamilies.H2_prototype_parameters(17, spin=1):  # long time
            ....:     X = VeeringTriangulationLinearFamilies.prototype_H2(a, b, c, e, mutable=True)
            ....:     X.set_canonical_labels()
            ....:     assert X in A1
        """
        #         (a,c)
        #           x-------x---------x  (a+b,c)
        #          /                 /
        #         /                 /
        # (0,0)  x-------x---------x  (b,0)
        #        |       |
        #        |       |
        #        |       |
        #        x-------x
        #     (0,-lbda)    (lbda,-lbda)
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.number_field.number_field import NumberField
        from sage.rings.qqbar import AA
        if not all(isinstance(x, numbers.Integral) for x in (a, b, c, e)):
            raise ValueError('a, b, c, e must be integers')
        a = ZZ(a)
        b = ZZ(b)
        c = ZZ(c)
        e = ZZ(e)
        if not (0 < b and 0 < c and 0 <= a < gcd(b, c) and c + e < b and gcd([a, b, c, e]) == 1):
            raise ValueError('a, b, c, e must satisfy: 0 < b, 0 < c, 0 <= a < gcd(b, c), c + e < b and gcd(a, b, c, e) = 1')
        D = e**2 + 4*b*c
        if D.is_square():
            K = QQ
            sqrtD = D.sqrtrem()[0]
        else:
            E = D.squarefree_part()
            f = (D // E).sqrtrem()[0]
            x = QQ['x'].gen()
            K = NumberField(x**2 - E, 'sqrt%d' % E, embedding=AA(E).sqrt())
            sqrtD = f * K.gen()

        assert sqrtD**2 == D

        lbda = (e + sqrtD) / 2

        fp = "(~0,5,6)(0,2,4)(~2,3,1)"
        cols = "BBRRRRR"
        vt = VeeringTriangulation(fp, cols, mutable=mutable)

        sx = (lbda, b-lbda, a, b+a-lbda, lbda+a, 0, lbda)
        sy = (0, 0, c, c, c, lbda, lbda)

        return VeeringTriangulationLinearFamily(vt, [sx, sy], mutable=mutable)

    @staticmethod
    def prototype_H1_1(a, b, c, e, mutable=False):
        r"""
        EXAMPLES::

            sage: from veerer.linear_family import VeeringTriangulationLinearFamilies
            sage: from veerer.automaton import GeometricAutomatonSubspace

            sage: X9 = VeeringTriangulationLinearFamilies.prototype_H1_1(0, 2, 1, -1)
            sage: X9.base_ring()
            Rational Field
            sage: X9.geometric_automaton()  # long time
            Geometric veering linear constraint automaton with 1244 vertices
        """
        #         (a+r,c)         (a+b,c)
        #           x-------x------o--x  (a+b+r,c)
        #          /                 /
        #         /                 /
        # (0,0)  o-------o--x------o  (b+r,0)
        #        \        \   (lbda+r,0)
        #         \        \
        #          \        \
        #           x--------x
        #     (r,-lbda)    (lbda+r,-lbda)
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.number_field.number_field import NumberField
        from sage.rings.qqbar import AA
        if not all(isinstance(x, numbers.Integral) for x in (a, b, c, e)):
            raise ValueError('a, b, c, e must be integers')
        a = ZZ(a)
        b = ZZ(b)
        c = ZZ(c)
        e = ZZ(e)
        if not (0 < b and 0 < c and 0 <= a < gcd(b, c) and c + e < b and gcd([a, b, c, e]) == 1):
            raise ValueError('a, b, c, e must satisfy: 0 < b, 0 < c, 0 <= a < gcd(b, c), c + e < b and gcd(a, b, c, e) = 1')
        D = e**2 + 4*b*c
        if D.is_square():
            K = QQ
            sqrtD = D.sqrtrem()[0]
        else:
            E = D.squarefree_part()
            f = (D // E).sqrtrem()[0]
            x = QQ['x'].gen()
            K = NumberField(x**2 - E, 'sqrt%d' % E, embedding=AA(E).sqrt())
            sqrtD = f * K.gen()

        assert sqrtD**2 == D

        lbda = (e + sqrtD) / 2

        fp = "(0,3,8)(~0,5,6)(~3,4,2)(~4,1,7)"
        cols = "BBBRRRRRR"
        vt = VeeringTriangulation(fp, cols, mutable=mutable)

        sx = (lbda, 0, b-lbda, a, a+b-lbda, 0, lbda, a+b-lbda, lbda+a)
        sy = (0, 0, 0, c, c, lbda, lbda, c, c)
        sr = (0, 1, -1, 1, 0, -1, -1, -1, 1)

        return VeeringTriangulationLinearFamily(vt, [sx, sy, sr], mutable=mutable)

    @staticmethod
    def L_shaped_surface(a1, a2, b1, b2, t1=0, t2=0):
        vt, s, t = VeeringTriangulations.L_shaped_surface(a1, a2, b1, b2, t1, t2)
        return VeeringTriangulationLinearFamily(vt, matrix([s, t]))

    @staticmethod
    def triangle_3_4_13_unfolding_orbit_closure():
        r"""
        EXAMPLES::

            sage: from veerer.linear_family import VeeringTriangulationLinearFamilies
            sage: f = VeeringTriangulationLinearFamilies.triangle_3_4_13_unfolding_orbit_closure()
            sage: f
            VeeringTriangulationLinearFamily("(0,9,~8)(1,8,2)(3,11,~10)(4,~14,15)(5,~15,12)(6,~16,13)(7,~0,22)(10,14,~9)(16,~12,~4)(17,20,~18)(18,~5,~23)(19,~22,~21)(21,~13,23)(~20,~6,~17)(~19,~7,~1)(~11,~2,~3)", "BBRBRBRRRRRRRRBRBBRRRBRR", [(1, phi, 0, 0, 0, 1, 0, 0, -phi, -phi - 1, 0, 0, -phi, -phi, phi + 1, -phi - 1, phi, 0, 0, -phi, 0, phi - 1, -1, -1), (0, 0, 1, 0, 0, 0, phi - 1, 0, 1, 1, 1, 1, 0, phi - 1, 0, 0, 0, 0, phi - 1, 0, phi - 1, 0, 0, phi - 1), (0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, phi - 1, 0, 0, -phi + 1, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0, phi - 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, phi - 1, 0, 0, phi - 1, 0)])

            sage: vt.stratum()  # optional - surface_dynamics
            Q_4(11, 1)
        """
        from sage.rings.rational_field import QQ
        from sage.rings.number_field.number_field import NumberField
        from sage.rings.qqbar import AA

        x = QQ['x'].gen()

        fp = "(0,9,~8)(1,8,2)(10,14,~9)(3,11,~10)(~2,~3,~11)(~14,15,4)(12,5,~15)(~4,16,~12)(13,6,~16)(23,21,~13)(18,~5,~23)(17,20,~18)(~6,~17,~20)(~0,22,7)(~21,19,~22)(~7,~1,~19)"
        cols = "BBRBRBRRRRRRRRBRBBRRRBRR"
        vt = VeeringTriangulation(fp, cols)

        # equations (beyond switches)
        # B = phi A
        # S = phi T
        # C = phi D
        # U = phi V
        # where A = 0, B = 1, C = 2, D = 6, S = 3, T = 17, U = 4, V = 7
        K = NumberField(x**2 - x - 1, 'phi', embedding=(AA(5).sqrt() + 1)/2)
        phi = K.gen()
        L = LinearExpressions(K)
        cs = ConstraintSystem()
        vt._set_switch_conditions(cs.insert, [L.variable(e) for e in range(24)])
        A = L.variable(0)
        B = L.variable(1)
        C = L.variable(2)
        C = L.variable(2)
        D = L.variable(6)
        S = L.variable(3)
        T = L.variable(17)
        U = L.variable(4)
        V = L.variable(7)
        cs.insert(B == phi * A)
        cs.insert(S == phi * T)
        cs.insert(C == phi * D)
        cs.insert(U == phi * V)

        return VeeringTriangulationLinearFamily(vt, cs.linear_generators_matrix())

