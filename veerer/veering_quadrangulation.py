r"""
Veering quadrangulations

A Veering quadrangulation is a purple (or square-tiled) veering triangulation.
This module implements the Abelian case which can conveniently be encoded by
a pair of permutations ``pr`` (for right) and ``pl`` (for left) of
`\{1, ..., n\}` (where `n` is the number of quadrilaterals). In order to keep
the structure of a quadrangulation, only simultaneous flips inside a right or
left staircase (topological cylinders) are allowed. As shown by Ferenczi and
Zamboni [FeZa10]_, this is not a restriction in the stratum `H^{hyp}(2g-2)` or
`H^{hyp}(g-1,g-1)`.

Veering quadrangulations generalize Penner construction.

.. TODO:

    The general case of quadratic differential can be encoded with a triple
    of involutions without fixed points on symmetric group on n symbols
    (tau_R, tau_B, tau_P).
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2020-2021 Vincent Delecroix
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

import functools
import numbers

from sage.rings.all import ZZ, QQ, AA
from sage.rings.qqbar import number_field_elements_from_algebraics
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.matrix.special import identity_matrix, zero_matrix, block_matrix
from sage.libs.gap.libgap import libgap

from veerer.permutation import *

RIGHT = 1
LEFT = 2


@functools.total_ordering
class VeeringQuadrangulation:
    r"""
    A Veering quadrangulation.

    The class stores three attributes

    - ``_n`` -- the number of quadrilaterals

    - ``_pr`` -- the right permutation

    - ``_pl`` -- the left permutation
    """
    __slots__ = ['_n', '_pr', '_pl']

    def __init__(self, pr, pl, n=-1, check=True):
        self._pr = perm_init(pr, n)  # right permutation

        if n == -1:
            n = len(self._pr)
        elif isinstance(n, numbers.Integral):
            n = int(n)
        else:
            raise ValueError('invalid n')

        self._pl = perm_init(pl, n)  # left permutation
        self._n = n                  # number of quadrilaterals
        if check:
            self._check()

    def _check(self):
        if not isinstance(self._n, int) or self._n < 0:
            raise ValueError('invalid n')
        if not perm_check(self._pr, self._n):
            raise ValueError('invalid right permutation')
        if not perm_check(self._pl, self._n):
            raise ValueError('invalid left permutation')

    def copy(self):
        r"""
        Return a copy of this veering quadrangulation.
        """
        V = VeeringQuadrangulation.__new__(VeeringQuadrangulation)
        V._n = self._n
        V._pr = self._pr[:]
        V._pl = self._pl[:]
        return V

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation  # random output due to deprecation warnings from realalg
            sage: VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)") == VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            True
            sage: VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)") == VeeringQuadrangulation("(0,1)(2)", "(0,2,1)")
            False
            sage: VeeringQuadrangulation("(0,1,2)", "(0,2)(1)") == VeeringQuadrangulation("(0,1)(2)", "(0,2,1)")
            False
        """
        if type(self) is not type(other):
            raise TypeError
        return self._n == other._n and self._pr == other._pr and self._pl == other._pl

    def __lt__(self, other):
        r"""
        Comparison operator <

        We first compare the number of cycles in the right permutation `pr`,
        then the number of cycles in the left permutation `pl`. If they are
        equal we do a plain comparison of the permutations.
        """
        if type(self) is not type(other):
            raise TypeError

        if self._n < other._n:
            return True
        elif self._n > other._n:
            return False

        s_rct = perm_cycle_type(self._pr)
        o_rct = perm_cycle_type(other._pr)
        if len(s_rct) < len(o_rct):
            return True
        elif len(s_rct) > len(o_rct):
            return False

        s_lct = perm_cycle_type(self._pl)
        o_lct = perm_cycle_type(other._pl)
        if len(s_lct) < len(o_lct):
            return True
        elif len(s_lct) > len(o_rct):
            return False

        return (self._pr, self._pl) < (other._pr, other._pl)

    def __repr__(self):
        return "VeeringQuadrangulation(\"{}\", \"{}\")".format(
                 perm_cycle_string(self._pr, n=self._n),
                 perm_cycle_string(self._pl, n=self._n))

    def to_origami(self):
        from .features import surface_dynamics_feature
        surface_dynamics_feature.require()

        from surface_dynamics import Origami
        return Origami(list(self._pr), list(self._pl), as_tuple=True, check=False)

    def stratum(self):
        r"""
        Return the ambient stratum of Abelian differential.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: q = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: q.stratum() # optional - surface_dynamics
            H_2(2)
            sage: q = VeeringQuadrangulation("(0,1,2,3,4)", "(0)(1,2)(3,4)")
            sage: q.stratum() # optional - surface_dynamics
            H_3(4)
        """
        return self.to_origami().stratum()

    def stratum_component(self):
        r"""
        Return the ambient stratum component.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: q = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: q.stratum_component() # optional - surface_dynamics
            H_2(2)^hyp
            sage: q = VeeringQuadrangulation("(0,1,2,3,4)", "(0)(1,2)(3,4)")
            sage: component = q.stratum_component() # optional - surface_dynamics  # random output due to deprecation warnings from surface dynamics
            sage: component  # optional - surface_dynamics
            H_3(4)^odd
        """
        return self.to_origami().stratum_component()

    def to_string(self):
        r"""
        String encoding of a veering quadrangulation.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)").to_string()
            '3_102_210'
        """
        return uint_base64_str(self._n) + '_' + perm_base64_str(self._pr) + '_' + perm_base64_str(self._pl)

    @staticmethod
    def from_string(s, check=True):
        r"""
        Build a veering quadrangulation from a string encoding.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: VeeringQuadrangulation.from_string('3_102_210')
            VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
        """
        n, pr, pl = s.split('_')
        n = uint_from_base64_str(n)
        pr = perm_from_base64_str(pr, n)
        pl = perm_from_base64_str(pl, n)
        return VeeringQuadrangulation(pr, pl, n)

    def __hash__(self):
        return hash((tuple(self._pr), tuple(self._pl)))

    def r_cycles(self):
        r"""
        Return the cycles of the right permutation.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: q = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: q.r_cycles()
            [[0, 1], [2]]
        """
        return perm_cycles(self._pr, self._n)

    def l_cycles(self):
        r"""
        Return the cycles of the left permutation.
        """
        return perm_cycles(self._pl, self._n)

    def r_orbit(self, i):
        r"""
        Return the orbit of ``i`` under the right permutation.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: q = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: q.r_orbit(1)
            [1, 0]
        """
        return perm_orbit(self._pr, i)

    def r_orbit_size(self, i):
        r"""
        Return the size of the orbit of ``i`` under the right permutation.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: q = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: q.r_orbit_size(1)
            2
        """
        return perm_orbit_size(self._pr, i)

    def l_orbit(self, i):
        r"""
        Return the orbit of ``i`` under the left permutation.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: q = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: q.l_orbit(2)
            [2, 0]
        """
        return perm_orbit(self._pl, i)

    def r_cycles_representatives(self, representatives=None):
        r"""
        Return the cycle representatives and cycle lengths of the right
        permutation.

        If the argument ``representatives`` is given, then only return the
        indices whose ``pr``-orbit intersect the ``pl``-orbit of
        ``representatives``.

        OUTPUT: a pair of lists ``(elements, sizes)`` where ``elements``
        are the minimum in the cycle of the right permutation and ``sizes``
        are the size of the corresponding orbit.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation

            sage: Q = VeeringQuadrangulation('(0,1)(2)', '(0,1,2)')
            sage: Q.r_cycles_representatives([0])
            ([0, 2], [2, 1])

            sage: Q = VeeringQuadrangulation('(0,1)(2,3)(4,5)(6,7)', '(0,3,5,7)')
            sage: Q.r_cycles_representatives()
            ([0, 3, 5, 7], [2, 2, 2, 2])

            sage: Q = VeeringQuadrangulation('(0,1)(2,3)', '(0,1,2)(3)')
            sage: Q.r_cycles_representatives([0])
            ([0, 2], [2, 2])
            sage: Q.r_cycles_representatives([3])
            ([3], [2])
            sage: Q.symmetry()
            sage: Q.r_cycles_representatives([0])
            ([0], [3])
            sage: Q.r_cycles_representatives([2])
            ([2, 3], [3, 1])
        """
        n = self._n
        pr = self._pr
        pl = self._pl
        ans = []
        sizes = []
        seen = [0] * n
        if representatives is None:
            representatives = range(n)
        for i in representatives:
            while seen[i] != -1:
                seen[i] = -1
                j = pr[i]
                m = 1
                while not seen[j]:
                    seen[j] = 1
                    j = pr[j]
                    m += 1
                if i == j:
                    ans.append(i)
                    sizes.append(m)
                i = pl[i]
        return ans, sizes

    def symmetry(self):
        r"""
        Exchange inplace the right and left permutations.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: Q = VeeringQuadrangulation('(0,1)(2,3,4)(5,6,7)', '(0,5,3)(1)(2,6,7)')
            sage: Q
            VeeringQuadrangulation("(0,1)(2,3,4)(5,6,7)", "(0,5,3)(1)(2,6,7)(4)")
            sage: Q.symmetry()
            sage: Q
            VeeringQuadrangulation("(0,5,3)(1)(2,6,7)(4)", "(0,1)(2,3,4)(5,6,7)")
        """
        self._pl, self._pr = self._pr, self._pl

    def relabel(self, p):
        r"""
        Relabel according to the permutation ``p``.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation

            sage: Q1 = VeeringQuadrangulation('(0,1)(2,3,4)(5,6,7)', '(0,5,3)(1)(2,6,7)')
            sage: Q2 = VeeringQuadrangulation('(0,1,7)(2,6,4)(3,5)', '(0,4,2)(1,3,6)(5)(7)')
            sage: Q1.is_isomorphic(Q2)  # optional - surface_dynamics
            True
            sage: ans,p = Q1.is_isomorphic(Q2, True)  # optional - surface_dynamics
            sage: Q1.relabel(p)  # optional - surface_dynamics
            sage: Q1 == Q2  # optional - surface_dynamics
            True
        """
        n = self._n
        if not perm_check(p, n):
            p = perm_init(p, n)
            if not perm_check(p, n):
                raise ValueError('invalid relabeling permutation')

        self._pr = perm_conjugate(self._pr, p)
        self._pl = perm_conjugate(self._pl, p)

    def r_staircase_move(self, i, m=1):
        r"""
        Perform a right staircase move in the cylinder containing the quadrilateral number ``i``.

        INPUT:

        - ``i`` -- index of a quadrilateral

        - ``m`` (integer, default 1) -- multiplicity of the staircase move to perform

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation

            sage: Q = VeeringQuadrangulation('(0,1)', '(0,2)', n=3)
            sage: Q.r_staircase_move(0)
            sage: Q
            VeeringQuadrangulation("(0,1)(2)", "(0,1,2)")
            sage: Q.r_staircase_move(0)
            sage: Q
            VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
        """
        pr = self._pr
        pl = self._pl

        # compute orbit size
        s = perm_orbit_size(pr, i)

        # not very efficient... we would better go by step of size m
        # and count how many of them we need (gcd computation)
        for _ in range(m % s):
            pli = pl[i]
            j = i
            jj = pr[j]
            while jj != i:
                pl[j] = pl[jj]
                j = jj
                jj = pr[j]
            pl[j] = pli

    def r_matrix(self, representatives, multiplicities):
        r"""
        Return the matrix associated with the right staircase moves given
        by ``representatives`` and ``multiplicities``.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: Q = VeeringQuadrangulation('(0,1)', '(0,2)', n=3)
            sage: Q.r_matrix([0, 2], [1, 0])
            [0 0 1]
            [0 1 0]
            [0 0 0]
            sage: Q.r_matrix([0, 2], [0, 1])
            [0 0 0]
            [0 0 0]
            [1 0 0]
            sage: Q.r_matrix([0, 2], [1, 1])
            [0 0 1]
            [0 1 0]
            [1 0 0]
        """
        pr = self._pr
        pl = self._pl
        A = matrix(ZZ, self._n)
        for i,m in zip(representatives, multiplicities):
            if not m:
                continue

            # compute the pr-orbit size
            s = perm_orbit_size(pr, i)
            a = m // s
            mm = m - a * s
            assert mm >= 0

            # multiplicative step
            if a:
                j = i
                while True:
                    jj = j
                    while True:
                        A[j, pl[jj]] = a
                        jj = pr[jj]
                        if jj == j:
                            break
                    j = pr[j]
                    if j == i:
                        break

            # additive step
            if mm:
                j = i
                while True:
                    jj = j
                    for k in range(mm):
                        A[j, pl[jj]] += 1
                        jj = pr[jj]
                    j = pr[j]
                    if j == i:
                        break
        return A

    def l_staircase_move(self, i):
        pr = self._pr
        pl = self._pl
        pri = pr[i]
        j = i
        jj = pl[j]
        while jj != i:
            pr[j] = pr[jj]
            j = jj
            jj = pl[j]
        pr[j] = pri

    def is_isomorphic(self, other, certificate=False):
        r"""
        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation

            sage: Q = VeeringQuadrangulation('(0,1)(2,3,4)(5,6,7)', '(0,5,3)(1)(2,6,7)')
            sage: QQ = Q.copy()
            sage: QQ.relabel('(0,5,2,4)(3,6)')
            sage: Q.is_isomorphic(QQ)  # optional - surface_dynamics
            True
        """
        if type(self) is not type(other):
            raise TypeError
        if self._n != other._n:
            return (False, None) if certificate else False

        sctr = perm_cycle_type(self._pr, self._n)
        octr = perm_cycle_type(other._pr, other._n)
        sctr.sort()
        octr.sort()
        if sctr != octr:
            return (False, None) if certificate else False

        sctl = perm_cycle_type(self._pl, self._n)
        octl = perm_cycle_type(other._pl, other._n)
        sctl.sort()
        octl.sort()
        if sctl != octl:
            return (False, None) if certificate else False

        o1 = self.to_origami()
        o2 = other.to_origami()

        if certificate:
            ans, cert = o1.is_isomorphic(o2, True)
            if ans:
                return ans, perm_init(cert, self._n)
            else:
                return ans, None
        else:
            return o1.is_isomorphic(o2)


class FlatVeeringQuadrangulation:
    r"""
    A flat structure on a veering quadrangulation.

    EXAMPLES::

        sage: from veerer.veering_quadrangulation import FlatVeeringQuadrangulation

        sage: x = polygen(QQ)
        sage: K.<sqrt5> = NumberField(x^2 - 5, embedding=AA(5).sqrt())
        sage: phi = (1+sqrt5) / 2
        sage: zr = [vector(K, 2, [1,phi])] * 4
        sage: zl = [vector(K, 2, [-phi,1])] * 4
        sage: pir = '(0,1)(2,3)'
        sage: pil = '(0,1,2,3)'
        sage: FlatVeeringQuadrangulation(pir, pil, zr, zl, 4)
        FlatVeeringQuadrangulation("(0,1)(2,3)", "(0,1,2,3)", [(1, 1/2*sqrt5 + 1/2), (1, 1/2*sqrt5 + 1/2), (1, 1/2*sqrt5 + 1/2), (1, 1/2*sqrt5 + 1/2)], [(-1/2*sqrt5 - 1/2, 1), (-1/2*sqrt5 - 1/2, 1), (-1/2*sqrt5 - 1/2, 1), (-1/2*sqrt5 - 1/2, 1)])

        sage: zr = [(5,7), (3,2), (5,7)]
        sage: zl = [(-3,2), (-3,2), (-2,5)]
        sage: pir = '(0,1)(2)'
        sage: pil = '(0,2)(1)'
        sage: F = FlatVeeringQuadrangulation(pir, pil, zr, zl, 3)
    """
    __slots__ = ['_n', '_pr', '_pl', '_zr', '_zl', '_V', '_K']

    def __init__(self, pr, pl, zr, zl, n=-1, base_ring=None, check=True):
        self._pr = perm_init(pr, n)
        if n == -1:
            n = len(self._pr)
        elif isinstance(n, numbers.Integral):
            n = int(n)
        else:
            raise TypeError('n must be integral')
        self._pl = perm_init(pl, n)
        if base_ring is None:
            from sage.structure.sequence import Sequence
            S = Sequence([vector(v) for v in zr] + [vector(v) for v in zl])
            self._V = S.universe()
            self._K = self._V.base_ring()
            z = list(S)
            zr = [v.__copy__() for v in z[:n]]
            zl = [v.__copy__() for v in z[n:]]
        else:
            self._V = zl[0].parent()
            self._K = self._V.base_ring()
            zr = [v.__copy__() for v in zr]
            zl = [v.__copy__() for v in zl]
        self._zr = zr
        self._zl = zl
        self._n = n
        if check:
            self._check()

    def _check(self):
        if not isinstance(self._n, int) or self._n < 0:
            raise ValueError('invalid n={}'.format(self._n))
        if not perm_check(self._pr, self._n):
            raise ValueError('invalid right permutation')
        if not perm_check(self._pl, self._n):
            raise ValueError('invalid left permutation')
        if len(self._zl) != self._n or len(self._zr) != self._n:
            raise ValueError('invalid vectors')
        for i in range(self._n):
            if self._zr[i][0] <= 0:
                raise ValueError('non-positive real part of right vector')
            if self._zl[i][0] >= 0:
                raise ValueError('non-negative real part of left vector')
            if self._zr[i][1] <= 0 or self._zl[i][1] <= 0:
                raise ValueError('non-positive imaginary part')
            d1 = self._zr[i] + self._zl[self._pr[i]]
            d2 = self._zl[i] + self._zr[self._pl[i]]
            if d1 != d2:
                raise ValueError('incompatible vectors for quadrilateral {}'.format(i))

    def copy(self):
        V = FlatVeeringQuadrangulation.__new__(FlatVeeringQuadrangulation)
        V._n = self._n
        V._pr = self._pr[:]
        V._pl = self._pl[:]
        V._V = self._V
        V._K = self._K
        V._zl = [v.__copy__() for v in self._zl]
        V._zr = [v.__copy__() for v in self._zr]
        return V

    def __eq__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._n == other._n and self._pr == other._pr and self._pl == other._pl and \
               self._zl == other._zl and self._zr == other._zr

    def __ne__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._n != other._n or self._pr != other._pr or self._pl != other._pl and \
               self._zl != other._zl or self._zr != other._zr

    def __repr__(self):
        return "FlatVeeringQuadrangulation(\"{}\", \"{}\", {}, {})".format(
                 perm_cycle_string(self._pr, n=self._n),
                 perm_cycle_string(self._pl, n=self._n),
                 self._zr,
                 self._zl)

    def is_r_slanted(self, i):
        j = i
        while True:
            if self._zr[i][0] + self._zl[self._pr[i]][0] >= 0:
                return False
            j = self._pr[j]
            if j == i:
                break
        return True

    def relabel(self, p):
        r"""
        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation
            sage: Q1 = VeeringQuadrangulation('(0,1)(2,3,4)(5,6,7)', '(0,5,3)(1)(2,6,7)')
            sage: Q2 = VeeringQuadrangulation('(0,1,7)(2,6,4)(3,5)', '(0,4,2)(1,3,6)(5)(7)')
            sage: Q1.is_isomorphic(Q2)  # optional - surface_dynamics
            True
            sage: ans,p = Q1.is_isomorphic(Q2, True)  # optional - surface_dynamics
            sage: Q1.relabel(p)  # optional - surface_dynamics
            sage: Q1 == Q2  # optional - surface_dynamics
            True

            sage: from veerer.veering_quadrangulation import FlatVeeringQuadrangulation
            sage: zr = [(5,7), (3,2), (5,7)]
            sage: zl = [(-3,2), (-3,2), (-2,5)]
            sage: pir = '(0,1)(2)'
            sage: pil = '(0,2)(1)'
            sage: F = FlatVeeringQuadrangulation(pir, pil, zr, zl, 3)
            sage: F.relabel('(0,1,2)')
            sage: F._check()
        """
        n = self._n
        if not perm_check(p, n):
            p = perm_init(p, n)
            if not perm_check(p, n):
                raise ValueError('invalid relabeling permutation')

        self._pr = perm_conjugate(self._pr, p)
        self._pl = perm_conjugate(self._pl, p)
        perm_on_list(p, self._zr, self._n)
        perm_on_list(p, self._zl, self._n)

    def well_slanted_r_staircases(self):
        return [i[0] for i in perm_cycles(self._pr, self._n) if self.is_r_slanted(i[0])]

    def flip(self):
        r"""
        Perform one step of induction and return the combinatorics.
        """
        representatives = []
        multiplicities = []
        for c in perm_cycles(self._pr, self._n):
            i = c[0]
            mi = 0
            while self.is_r_slanted(i):
                self.r_staircase_move(i, check=False)
                mi += 1
            representatives.append(i)
            multiplicities.append(mi)
        self.symmetry()
        return representatives, multiplicities

    def r_staircase_move(self, i, check=True):
        r"""
        Apply the right staircase move on the right staircase containing ``i``.
        """
        if check and not self.is_r_slanted(i):
            raise ValueError('the right staircase containing {} is not right slanted'.format(i))

        n = self._n
        pr = self._pr
        pl = self._pl
        zr = self._zr
        zl = self._zl

        j = i
        while True:
            zl[j] += zr[pl[j]]
            j = pr[j]
            if j == i:
                break
        pli = pl[i]
        j = i
        jj = pr[j]
        while jj != i:
            pl[j] = pl[jj]
            j = jj
            jj = pr[jj]
        pl[j] = pli
        if check:
            self._check()

    def symmetry(self):
        r"""
        Exchange left and right (non-orientable symmetry with respect to the vertical axis).
        """
        self._pr, self._pl = self._pl, self._pr
        self._zr, self._zl = self._zl, self._zr
        for i in range(self._n):
            self._zl[i][0] *= -1
            self._zr[i][0] *= -1

    def apply_diagonal_matrix(self, a, b):
        self._zl = [self._V((a*x, b*y)) for x,y in self._zl]
        self._zr = [self._V((a*x, b*y)) for x,y in self._zr]

    def normalize(self):
        r"""
        Apply a diagonal matrix so that sum of x-coordinates of right vectors is n
        and sum of y-coordinates of left vectors is n where n is the number of quadrilaterals
        """
        n = self._n
        sx = sum(self._zr[i][0] for i in range(n)) / n
        sy = sum(self._zl[i][1] for i in range(n)) / n
        for i in range(n):
            self._zr[i][0] /= sx
            self._zl[i][0] /= sx
            self._zr[i][1] /= sy
            self._zl[i][1] /= sy

    def plot(self):
        from sage.plot.graphics import Graphics
        from sage.plot.line import line2d
        from sage.plot.polygon import polygon2d
        from sage.plot.plot import graphics_array
        output = []
        n = self._n
        for i in range(n):
            G = Graphics()
            o = (0,0)
            r = self._zr[i]
            l = self._zl[i]
            d = self._zr[i] + self._zl[self._pr[i]]
            G += polygon2d([o, r, d, l], fill=True, color='gray', alpha=0.4)
            G += line2d([o,r], color='red')
            G += line2d([r,d], color='blue')
            G += line2d([o,l], color='blue')
            G += line2d([l,d], color='red')
            G.set_aspect_ratio(1)
            G.axes(True)
            output.append(G)
        return graphics_array(output, 1, self._n)

@functools.total_ordering
class VeeringQuadrangulationFlipSequence:
    r"""
    Flip sequence of veering quadrangulations.

    Contrarily to the case of triangulation in :mod:`flip_sequence`, the flips
    in the same direction, left or right, are grouped together in blocks.
    Hence, a flip sequence corresponds to a succession of right moves and
    axial symmetries along `x = y`.

    A block in a quadrangulation `Q` is determined by a multiplicity for each
    right staircase in `Q`. In the flip sequence, a staircase is encoded by any
    square that belongs in it.
    """
    def __init__(self, start, sequence=None, relabelling=None):
        if not isinstance(start, VeeringQuadrangulation):
            raise TypeError
        self._start = start.copy()
        self._end = start.copy()
        self._flips = []
        self._relabelling = perm_id(start._n)
        if sequence is not None:
            for representatives, multiplicities in sequence:
                self.append_flip(representatives, multiplicities)
        if relabelling is not None:
            self.append_relabelling(relabelling)

        # the sequence is normalized if:
        # - no 0 multiplicity
        # - no go through blocks to the left
        # - representatives are the mins
        self._normalized = False

    def __iter__(self):
        r"""
        Iterator through the list of triples ``(Q, representatives, multiplicities)``
        """
        q = self._start.copy()
        for step in self._flips:
            yield q, step[0], step[1]
            for i,m in zip(*step):
                q.r_staircase_move(i, m)
            q.symmetry()

    def append_flip(self, representatives, multiplicities):
        r"""
        Add a flip to this flip sequence.

        INPUT:

        - ``representatives`` -- list of quadrilateral number that identify the
          staircases in the final veering quadrangulation of this flip
          sequence
        - ``multiplicities`` -- list of positive integers of twist multiplicity
          for each staircase
        """
        if not any(multiplicities):
            raise ValueError('at least one multiplicity must be nonzero')
        for i, mult in zip(representatives, multiplicities):
            self._end.r_staircase_move(i, mult)
        self._end.symmetry()
        representatives = [perm_preimage(self._relabelling, i) for i in representatives]
        self._flips.append((tuple(representatives), tuple(multiplicities)))

    def append_relabelling(self, r):
        r"""
        Add a relabelling to this flip sequence.
        """
        end = self._end
        if not perm_check(r, end._n):
            r = perm_init(r, end._n)
            if not perm_check(r, end._n):
                raise ValueError('invalid relabelling permutation')
        end.relabel(r)
        self._relabelling = perm_compose(self._relabelling, r)

    def automorphism_group(self):
        r"""
        EXAMPLES::

            sage: from veerer.veering_quadrangulation import *

            sage: V = VeeringQuadrangulation("(0,1)(2,3)", "(0)(1,2)(3)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0, 2], [1, 1])
            sage: fp.append_flip([0], [1])
            sage: fp.append_flip([1,2],[1,1])
            sage: fp.append_relabelling("(0,3)(1,2)")
            sage: fp.automorphism_group()  # optional - surface_dynamics
            Subgroup generated by [(1,4)(2,3)] of (Symmetric group of order 4! as a permutation group)

            sage: V = VeeringQuadrangulation("(0,1)(2,3)", "(0)(1,2)(3)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0, 2], [1, 1])
            sage: fp.append_flip([0], [1])
            sage: fp.append_flip([1,2],[2,1])
            sage: fp.append_relabelling("(0,3)(1,2)")
            sage: fp.automorphism_group()  # optional - surface_dynamics
            Subgroup generated by [()] of (Symmetric group of order 4! as a permutation group)
        """
        A = self._start.to_origami().automorphism_group()
        S = A.ambient_group()
        A = A.gap()
        if A.Size().is_one():
            return S.subgroup(gap_group=A)

        B = self._end.to_origami().automorphism_group()
        A = A.Intersection(B)
        if A.Size() == 1:
            return S.subgroup(gap_group=A)

        # now check for symmetries of the steps
        q = self._start.copy()
        for representatives, multiplicities in self._flips:
            d = {}
            for i,m in zip(representatives, multiplicities):
                if m not in d:
                    d[m] = []
                d[m].extend([j+1 for j in perm_orbit(q._pr, i)])
                if m:
                    q.r_staircase_move(i, m)
            partition = list(d.values())
            for x in partition: x.sort()
            A = A.Stabilizer(partition, libgap.OnTuplesSets)
            if A.Size() == 1:
                return S.subgroup(gap_group=A)
            q.symmetry()

        return S.subgroup(gap_group=A)

    def relabel(self, r):
        r"""
        Relabel the whole flip sequence.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence
            sage: V = VeeringQuadrangulation("(0,1,2,3,4)", "(0)(1,4)(2,3)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0], [1])
            sage: fp.append_flip([2], [1])
            sage: fp.append_relabelling("(0,3,1,4,2)")
            sage: assert fp.is_loop() and fp.is_complete()

            sage: V = VeeringQuadrangulation("(0)(1,2)", "(0,1,2)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0],[1])
            sage: fp.append_flip([0],[1])
            sage: fp.append_relabelling(fp.find_closure())  # optional - surface_dynamics
            sage: assert fp.is_loop() and fp.is_complete()  # optional - surface_dynamics
            sage: fp.num_flips()  # optional - surface_dynamics
            4
        """
        self._normalized = False
        n = self._start._n
        if not perm_check(r, n):
            r = perm_init(r, n)
            if not perm_check(r, n):
                raise ValueError('invalid relabeling permutation')
        for i in range(len(self._flips)):
            representatives, multiplicities = self._flips[i]
            representatives = [r[x] for x in representatives]
            self._flips[i] = (tuple(representatives), multiplicities)
        self._start.relabel(r)
        self._end.relabel(r)
        self._relabelling = perm_compose(perm_compose(perm_invert(r), self._relabelling), r)

    def __eq__(self, other):
        if type(self) is not type(other):
            raise TypeError
        if self._start != other._start or \
           self._end != other._end or \
           self._relabelling != other._relabelling or \
           len(self._flips) != len(other._flips):
            return False

        self.normalize()
        other.normalize()
        return self._flips == other._flips

    def __lt__(self, other):
        if type(self) is not type(other):
            raise TypeError
        if self._start < other._start:
            return True
        elif self._start > other._start:
            return False
        self.normalize()
        other.normalize()
        return self._flips < other._flips

    def normalize(self):
        r"""
        Normalize the flip sequence inplace.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence

            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([1,2],[2,0])
            sage: fp.append_flip([0,1],[0,1])
            sage: fp.append_flip([0,2],[0,1])
            sage: fp.append_flip([2,1],[2,0])
            sage: fp
            VeeringQuadrangulationFlipSequence(VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)"), [((1, 2), (2, 0)), ((0, 1), (0, 1)), ((0, 2), (0, 1)), ((2, 1), (2, 0))], "(0)(1)(2)")

            sage: fp.matrix()
            [2 1 1|0 1 1]
            [0 2 1|1 0 0]
            [1 1 2|0 1 1]
            [-----+-----]
            [0 1 1|1 0 0]
            [0 1 1|0 1 0]
            [1 0 0|0 0 1]

            sage: fp.normalize()
            sage: fp
            VeeringQuadrangulationFlipSequence(VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)"), [((0,), (2,)), ((1,), (1,)), ((2,), (1,)), ((0,), (2,))], "(0)(1)(2)")
            sage: fp.matrix()
            [2 1 1|0 1 1]
            [0 2 1|1 0 0]
            [1 1 2|0 1 1]
            [-----+-----]
            [0 1 1|1 0 0]
            [0 1 1|0 1 0]
            [1 0 0|0 0 1]
        """
        if self._normalized:
            return
        q = self._start.copy()
        for i, (representatives, multiplicities) in enumerate(self._flips):
            representatives_and_multiplicities = [(min(perm_orbit(q._pr, i)), m) for i,m in zip(representatives, multiplicities) if m]
            representatives_and_multiplicities.sort()
            for j in range(len(representatives_and_multiplicities) - 1):
                i0, m0 = representatives_and_multiplicities[j]
                i1, m1 = representatives_and_multiplicities[j+1]
                if i0 == i1:
                    raise ValueError('invalid flip sequence: same staircase appears twice')
            representatives = tuple(i for i, m in representatives_and_multiplicities)
            multiplicities = tuple(m for i, m in representatives_and_multiplicities)
            self._flips[i] = (representatives, multiplicities)
            for i,m in zip(representatives, multiplicities):
                q.r_staircase_move(i, m)
            q.symmetry()
        self._normalized = True

    def is_relabelling_equivalent(self, other, certificate=False):
        r"""
        Test whether ``self`` and ``other`` are equivalent up to a relabelling.
        """
        if self.num_blocks() != other.num_blocks() or \
           self.num_flips() != other.num_flips():
            return (False, None) if certificate else False
        ans, p = self._start.is_isomorphic(other._start, certificate=True)
        if not ans:
            return (False, None) if certificate else False
        self.normalize()
        rself = self.copy()
        rself.relabel(p)

        q = rself._start
        n = q._n
        A = q.to_origami().automorphism_group()
        A = [perm_init([p(i)-1 for i in range(1,n+1)], n) for p in A]
        for a in A:
            rother = other.copy()
            rother.relabel(a)
            rother.normalize()
            assert rself._start == rother._start
            if rself == rother:
                return (True, perm_compose(p, a)) if certificate else True
        return (False, None) if certificate else False

    def is_conjugate(self, other, certificate=False):
        r"""
        Test whether ``self`` and ``other`` are equivalent up to cyclic rotation and relabelling.
        """
        if not self.is_loop():
            raise ValueError('self not a loop')
        if not other.is_loop():
            raise ValueError('other not a loop')

        if self.num_blocks() != other.num_blocks() or \
           self.num_flips() != other.num_flips():
            return False

        # From now on we assume that the flip sequence are reduced

        # 0. Compare relabelling invariant quantities
        l = len(self._flips)
        infpb = self.num_flips_per_blocks()
        imin, iperiod = least_rotation(infpb)
        jnfpb = other.num_flips_per_blocks()
        jmin, jperiod = least_rotation(jnfpb)
        assert self.is_loop()

        infpb = infpb[imin:] + infpb[:imin]
        jnfpb = jnfpb[jmin:] + jnfpb[:jmin]
        if infpb != jnfpb:
            return False
        period = iperiod

        # 1. Do compare
        rself = self.copy()
        assert rself.is_loop()
        rself.rotate(imin)
        rother = other.copy()
        rother.rotate(jmin)
        for k in range(l // period):
            ans, p = rself.is_relabelling_equivalent(rother, certificate=True)
            if ans:
                return (True, (jmin - imin + k * period, p)) if certificate else True
            rother.rotate(period)

        return (False, None) if certificate else False


    def copy(self):
        res = VeeringQuadrangulationFlipSequence.__new__(VeeringQuadrangulationFlipSequence)
        res._start = self._start.copy()
        res._end = self._end.copy()
        res._flips = self._flips[:]
        res._relabelling = self._relabelling[:]
        res._normalized = self._normalized
        return res

    def flips_string(self, latex_color=False, shift=False):
        q = self._start.copy()
        n = q._n
        l = []
        color = 0
        for istep,step in enumerate(self._flips):
            if latex_color:
                sstep = '\\textcolor{{{}}}'.format('red' if color%2 == 0 else 'blue')
                sstep += '{'
            else:
                sstep = ''
            cycles = q.r_cycles()
            index_to_cycle = [None] * n
            for c in cycles:
                for j in c:
                    index_to_cycle[j] = c
            RL = 'RL'
            for i,m in zip(*step):
                if m:
                    q.r_staircase_move(i, m)
                    if shift:
                        sstep += RL[istep%2] + '_{(' + ','.join(str(x+1) for x in index_to_cycle[i]) + ')}' + '^%d' %m
                    else:
                        sstep += RL[istep%2] + '_{(' + ','.join(map(str,index_to_cycle[i])) + ')}' + '^%d' %m

            if latex_color:
                sstep += '}'
            q.symmetry()
            l.append(sstep)
            color = 1 - color
        if latex_color:
            return '\\,'.join(l)
        else:
            return '|'.join(l)

    def str(self):
        if len(self._flips) % 2 == 1:
            self = self * self
        p = self._start
        n = p._n
        return perm_cycle_string(p._pr, True, n) + '; ' + \
               perm_cycle_string(p._pl, True, n) + '; ' + \
               self.flips_string() + '; ' + \
               perm_cycle_string(self._relabelling, False, n)

    def latex(self, right_first=True, shift=False):
        if len(self._flips) % 2 == 1:
            self = self * self
        p = self._start
        n = p._n
        csr = perm_cycle_string(p._pr, True, n, shift=shift)
        csl = perm_cycle_string(p._pl, True, n, shift=shift)
        if right_first:
            perm = csr + ' & ' + csl
        else:
            perm = csl + ' & ' + csr
        return perm + ' & ' +\
               self.flips_string(latex_color=True, shift=shift) + ' & ' + \
               perm_cycle_string(self._relabelling, False, n, shift=shift)

    def __repr__(self):
        return "VeeringQuadrangulationFlipSequence({}, {}, \"{}\")".format(
                self._start,
                self._flips,
                perm_cycle_string(self._relabelling, True, self._start._n))

    def _check(self):
        fp = VeeringQuadrangulationFlipSequence(self._start)
        for representatives, multiplicities in self._flips:
            fp.append_flip(representatives, multiplicities)
        fp.append_relabelling(self._relabelling)
        if fp._end != self._end:
            raise ValueError('incompatible stuff')

    def __mul__(self, other):
        if self._end != other._start:
            raise ValueError('left hand side should end with the start of right hand side')
        res = self.copy()
        for representatives, multiplicities in other._flips:
            res.append_flip(representatives, multiplicities)
        res.append_relabelling(other._relabelling)
        return res

    def __pow__(self, n):
        r"""
        TESTS::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence
            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0,2],[2,1])
            sage: fp.append_flip([0,1],[2,1])
            sage: fp**0
            VeeringQuadrangulationFlipSequence(VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)"), [], "(0)(1)(2)")
            sage: fp**1 == fp
            True
            sage: fp**2 == fp * fp
            True
        """
        if self._end != self._start:
            raise ValueError('can not take power of a non-loop')
        n = int(n)
        if n < 0:
            raise ValueError
        if n == 0:
            return VeeringQuadrangulationFlipSequence(self._start)
        res = self.copy()
        for _ in range(n - 1):
            res = res * self
        return res

    def num_blocks(self):
        r"""
        Return the number of blocks.
        """
        return len(self._flips)

    def num_flips(self):
        r"""
        Return the number of flips.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence

            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0,2],[2,1])
            sage: fp.append_flip([0,1],[2,1])
            sage: fp.num_flips()
            10

        The pseudo-Anosov with minimal dilatation in `H(2)` is also the one
        with the least possible number of flips::

            sage: V = VeeringQuadrangulation("(0)(1,2)", "(0,1,2)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0],[1])
            sage: fp.append_flip([0],[1])
            sage: fp.append_relabelling(fp.find_closure())  # optional - surface_dynamics
            sage: assert fp.is_loop() and fp.is_complete()  # optional - surface_dynamics
            sage: fp.num_flips()  # optional - surface_dynamics
            4
        """
        ans = 0
        for q, representatives, multiplicities in self:
            for i,m in zip(representatives, multiplicities):
                ans += m * len(perm_orbit(q._pr, i))
        return ans

    def num_staircase_moves(self):
        r"""
        Return the number of staircase moves.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence

            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0,2],[2,1])
            sage: fp.append_flip([0,1],[7,3])
            sage: fp.num_staircase_moves()
            13

        The pseudo-Anosov with minimal dilatation in `H(2)` is also the one
        with the least possible number of staircase moves::

            sage: V = VeeringQuadrangulation("(0)(1,2)", "(0,1,2)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0],[1])
            sage: fp.append_flip([0],[1])
            sage: fp.append_relabelling(fp.find_closure())  # optional - surface_dynamics
            sage: assert fp.is_loop() and fp.is_complete()  # optional - surface_dynamics
            sage: fp.num_staircase_moves()  # optional - surface_dynamics
            2
        """
        ans = 0
        for q, representatives, multiplicities in self:
            ans += sum(multiplicities)
        return ans

    def num_flips_per_blocks(self):
        r"""
        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence
            sage: V = VeeringQuadrangulation("(0)(1,2)", "(0,1,2)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0],[1])
            sage: fp.append_flip([0],[1])
            sage: fp.append_relabelling(fp.find_closure())  # optional - surface_dynamics
            sage: assert fp.is_loop() and fp.is_complete()  # optional - surface_dynamics
            sage: fp.num_flips_per_blocks()  # optional - surface_dynamics
            [((1, (3,)),), ((1, (1, 2, 2)),)]
        """
        q = self._start.copy()
        ans = []
        for representatives, multiplicities in self._flips:
            flip = []
            for i,m in zip(representatives, multiplicities):
                if m:
                    sizes = [len(perm_orbit(q._pl, j)) for j in perm_orbit(q._pr, i)]
                    jmin, p = least_rotation(sizes)
                    sizes = sizes[jmin:] + sizes[:jmin]
                    flip.append((m, tuple(sizes)))
                    q.r_staircase_move(i, m)
            flip.sort()
            ans.append(tuple(flip))
            q.symmetry()
        return ans

    def rotate(self, n=1):
        r"""
        Remove the first ``n`` flips at the beginning of the
        sequence and move it at the end.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence
            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0,2],[2,1])
            sage: fp.append_flip([0,1],[2,1])
            sage: assert fp.is_loop() and fp.is_complete()
            sage: sorted(fp.matrix().eigenvalues())
            [0.2277771042343813?, 0.5441132198971334?, 1, 1, 1.837852791352972?, 4.390256884515514?]
            sage: fp.rotate()
            sage: sorted(fp.matrix().eigenvalues())
            [0.2277771042343813?, 0.5441132198971334?, 1, 1, 1.837852791352972?, 4.390256884515514?]

            sage: V = VeeringQuadrangulation("(0,1,2)", "(0,1)(2)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0],[1])
            sage: fp.append_flip([0],[1])
            sage: fp.append_relabelling([2, 0, 1])
            sage: assert fp.is_loop() and fp.is_complete()
            sage: fp.matrix().charpoly()
            x^6 - x^4 - 3*x^3 - x^2 + 1
            sage: fp.rotate()
            sage: fp.matrix().charpoly()
            x^6 - x^4 - 3*x^3 - x^2 + 1
            sage: fp.rotate()
            sage: fp.matrix().charpoly()
            x^6 - x^4 - 3*x^3 - x^2 + 1

            sage: V = VeeringQuadrangulation("(0,1,2,3,4)", "(0)(1,4)(2,3)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0], [1])
            sage: fp.append_flip([2], [1])
            sage: fp.append_relabelling("(0,3,1,4,2)")
            sage: assert fp.is_loop() and fp.is_complete()
            sage: poly = fp.matrix().charpoly()
            sage: fp.rotate()
            sage: assert fp.matrix().charpoly() == poly
        """
        if not self.is_loop():
            raise ValueError
        if not self._flips:
            raise ValueError('empty flip sequence')
        for _ in range(n):
            step = self._flips.pop(0)
            for i, mult in zip(*step):
                self._start.r_staircase_move(i, mult)
            self._start.symmetry()
            self.append_flip(*step)

    def is_loop(self):
        r"""
        A flip sequence is a loop if the initial and final quadrangulations are the same.

        A flip sequence that is a loop gives rise to a mapping class of the
        underlying surface.
        """
        return self._start == self._end

    def flipped_edges(self):
        r"""
        Return the set of edges that are flipped along this flip sequence
        assuming it is a loop.
        """
        if self._start != self._end:
            raise ValueError('not a loop')
        n = self._start._n

        right_flipped = set()
        left_flipped = set()
        V = self._start.copy()
        for nstep, (representatives, multiplicities) in enumerate(self._flips):
            for i,m in zip(representatives, multiplicities):
                if not m:
                    continue
                for j in perm_orbit(V._pr, i):
                    orbit = perm_orbit(self._relabelling, j)
                    if len(self._flips) % 2:
                        if len(orbit) % 2:
                            right_flipped.update(orbit)
                            left_flipped.update(orbit)
                        else:
                            right_flipped.update(orbit[0::2])
                            left_flipped.update(orbit[1::2])
                    else:
                        right_flipped.update(orbit)
                V.r_staircase_move(i, m)
            right_flipped, left_flipped = left_flipped, right_flipped
            V.symmetry()
            if len(right_flipped) == n and len(left_flipped) == n:
                return right_flipped, left_flipped
        return right_flipped, left_flipped

    def is_complete(self):
        r"""
        Return whether this flip sequence is complete.

        A flip sequence is complete if it is a loop and if, by repeating this loop,
        all edges are flipped. Equivalently, if the associated mapping class is of
        pseudo-Anosov type.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence

            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0],[2])
            sage: fp.append_flip([0],[2])
            sage: assert fp.is_loop()
            sage: fp.is_complete()
            False
            sage: fp.matrix().is_primitive()
            False
            sage: fp.append_flip([2],[1])
            sage: fp.append_flip([1],[1])
            sage: assert fp.is_loop()
            sage: fp.is_complete()
            True
            sage: fp.matrix().is_primitive()
            True

            sage: V = VeeringQuadrangulation("(0,1,2)", "(0,1)(2)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0], [2])
            sage: fp.append_flip([1], [1])
            sage: fp.append_flip([1], [1])
            sage: fp.append_flip([2], [1])
            sage: assert fp.is_loop()
            sage: fp.is_complete()
            False
            sage: fp.matrix().is_primitive()
            False

            sage: V = VeeringQuadrangulation("(0,1,2)", "(0,1)(2)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0], [1])
            sage: fp.append_flip([0], [1])
            sage: fp.append_flip([0], [2])
            sage: fp.append_flip([2], [1])
            sage: assert fp.is_loop()
            sage: fp.is_complete()
            False
            sage: fp.matrix().is_primitive()
            False

            sage: V = VeeringQuadrangulation("(0,1,2)", "(0,1)(2)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0], [1])
            sage: fp.append_flip([0,1], [1,0])
            sage: fp.append_flip([0], [2])
            sage: fp.append_flip([0,2], [0,1])
            sage: assert fp.is_loop()
            sage: fp.is_complete()
            False
            sage: fp.matrix().is_primitive()
            False
        """
        right_flipped, left_flipped = self.flipped_edges()
        return len(right_flipped) == len(left_flipped) == self._start._n

    def reduce(self):
        r"""
        Put self in standard form.
        """
        if not self.is_loop():
            raise ValueError("must be a loop")
        positions = [[0] * n, [1] * n]
        new_representatives = []
        new_multiplicities = []
        parity = 0
        for q, representatives, multiplicities in self:
            for i,m in zip(representatives, multiplicities):
                if m:
                    pos = max(positions[parity][j] for j in perm_orbit(q._pr, i))
                    while len(flips) <= pos:
                        flips.append([])
                    new_representatives[pos].append(i)
                    new_multiplicities[pos].append(m)
                    for j in perm_orbit(q._pr, i):
                        positions[1-parity][j] = max(positions[1-parity][j], pos + 1)
            parity = 1 - parity
        # We need to figure out one starting point...
        # at this stage, all levels could be wrong because of the starting point!!!!
        return new_representatives, new_multiplicities

    def is_reduced_at_end(self):
        r"""
        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence

            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0,2],[2,1])
            sage: fp.append_flip([0,1],[2,1])
            sage: fp.append_flip([0,2],[0,1])
            sage: fp.append_flip([0,1],[0,1])
            sage: fp.is_reduced_at_end()
            False

            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0,2],[2,1])
            sage: fp.append_flip([0,1],[2,1])
            sage: fp.is_reduced_at_end()
            True

            sage: V = VeeringQuadrangulation("(0)(1,2)", "(0,1,2)")
            sage: flips = [((0, 1), (0, 1)), ((1, 2), (1, 0)), ((1,), (1,))]
            sage: relabelling = "(0,1)(2)"
            sage: fp = VeeringQuadrangulationFlipSequence(V, flips, relabelling)
            sage: x = polygen(QQ)
            sage: p = (x - 1) * (x + 1) * (x^4 - 3*x^3 - x^2 + 3*x + 1)
            sage: for _ in range(6):
            ....:     assert fp.is_loop() and fp.is_complete() and fp.is_reduced_at_end()
            ....:     assert fp.matrix().charpoly() == p
            ....:     fp.rotate()

            sage: V = VeeringQuadrangulation("(0,1,2)", "(0)(1,2)")
            sage: flips = [((0,), (1,)), ((0, 1), (1, 0)), ((0, 2), (0, 1))]
            sage: relabelling = "(0)(1,2)"
            sage: fp = VeeringQuadrangulationFlipSequence(V, flips, relabelling)
            sage: fp.is_reduced_at_end()
            True

            sage: V = VeeringQuadrangulation("(0)(1,2)", "(0,1)(2)")
            sage: flips = [((0, 1), (1, 0)), ((0,), (2,)), ((0, 1), (1, 0))]
            sage: relabelling = "(0,2)(1)"
            sage: fp = VeeringQuadrangulationFlipSequence(V, flips, relabelling)
            sage: fp.is_reduced_at_end()
            False

            sage: V = VeeringQuadrangulation("(0)(1,2)", "(0,1)(2)")
            sage: flips = [((0, 1), (2, 0)), ((0,), (1,)), ((0,), (1,)), ((0, 1), (1, 0))]
            sage: relabelling = "(0,1,2)"
            sage: fp = VeeringQuadrangulationFlipSequence(V, flips, relabelling)
            sage: fp.is_reduced_at_end()
            True
        """
        if not self.is_loop():
            raise ValueError('the path must be a loop')

        # check compatibility of first move and completeness
        allowed_last = set()
        for i,m in zip(*self._flips[-1]):
            if m:
                i = self._relabelling[i]
                allowed_last.update(perm_orbit(self._end._pl, i))
        for i,m in zip(*self._flips[0]):
            if m and not any(j in allowed_last for j in perm_orbit(self._start._pr, i)):
                return False
        return True

    def find_closure(self):
        r"""
        Return an isomorphism identifying the end with start.
        """
        ans, p = self._end.is_isomorphic(self._start, True)
        return p if ans else None

    # dilatations should be at the same time
    # PF(Q1 + Q0) = PF(P1 + P0)
    def matrix_by_blocks(self):
        n = self._start._n
        P0 = identity_matrix(n)
        P1 = zero_matrix(n)
        Q0 = zero_matrix(n)
        Q1 = identity_matrix(n)

        # p_{i+1} = A p_i + p_{i-1}
        # q_{i+1} = A q_i + q_{i-1}
        V = self._start.copy()
        for representatives, multiplicities in self._flips:
            A = V.r_matrix(representatives, multiplicities)
            P0, P1 = P1, A*P1 + P0
            Q0, Q1 = Q1, A*Q1 + Q0
            for i, mult in zip(representatives, multiplicities):
                V.r_staircase_move(i, mult)
            V.symmetry()

        swap = sage.matrix.matrix0.Matrix.swap_rows
        perm_on_list(self._relabelling, P0, n, swap)
        perm_on_list(self._relabelling, P1, n, swap)
        perm_on_list(self._relabelling, Q0, n, swap)
        perm_on_list(self._relabelling, Q1, n, swap)
        return Q1, P1, Q0, P0

    def matrix(self):
        r"""
        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence

        A square-tiled example with coordinates in `Q[\sqrt{2}]`::

            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0,2],[2,2])
            sage: fp.append_flip([0,1],[2,2])
            sage: fp.is_loop() and fp.is_complete()
            True
            sage: fp.matrix()
            [3 1 1|0 1 1]
            [0 3 2|2 0 0]
            [2 1 2|0 1 1]
            [-----+-----]
            [0 1 1|1 0 0]
            [0 1 1|0 1 0]
            [2 0 0|0 0 1]

        Trace field `Q[\sqrt{5}]`::

            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0,2],[2,1])
            sage: fp.append_flip([0,1],[2,1])
            sage: fp.is_loop() and fp.is_complete()
            True
        """
        Q1, P1, Q0, P0 = self.matrix_by_blocks()
        return block_matrix(2, [Q1,P1,Q0,P0])

    def _self_similar_data(self):
        if len(self._flips) % 2:
            path = self * self
        else:
            path = self
        M = path.matrix()
        poly = M.charpoly()
        K, a, _ = number_field_elements_from_algebraics(max(poly.roots(AA, False)), minimal=True, embedded=True)
        assert a > 1, poly
        assert poly(a).is_zero(), poly
        assert poly(1/a).is_zero(), poly

        heights = (M-a).right_kernel_matrix()[0]
        if heights[0] < 0:
            heights = -heights
        widths = (M-1/a).right_kernel_matrix()[0]
        if widths[0] < 0:
            widths -= widths

        return (a, widths, heights)

    def self_similar_surface(self):
        r"""
        Return a pair ``(flat surface, dilatation)``.

        EXAMPLES::

            sage: from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence

            sage: V = VeeringQuadrangulation("(0,1)(2)", "(0,2)(1)")
            sage: fp = VeeringQuadrangulationFlipSequence(V)
            sage: fp.append_flip([0,2],[2,2])
            sage: fp.append_flip([0,1],[2,2])
            sage: F, a = fp.self_similar_surface()
            sage: FF = F.copy()
            sage: FF.flip()
            ([0, 2], [2, 2])
            sage: FF.flip()
            ([0, 1], [2, 2])
            sage: FF.apply_diagonal_matrix(a, 1/a)
            sage: F == FF
            True
        """
        n = self._start._n
        a, widths, heights = self._self_similar_data()
        K = a.parent()
        z = [vector(K, (x,y)) for x,y in zip(widths, heights)]
        zr = z[:n]
        zl = z[n:]
        return FlatVeeringQuadrangulation(self._start._pr, self._start._pl, zr, zl, base_ring=K), a

def weighted_composition_iterator(
        n,
        weights=None,
        max_weight=None,
        max_sum=None,
        max_part=None):
    r"""
    Iterator through the triple ``(c, w, s)`` where ``c`` is a weighted composition
    of maximum weight at most ``max_weight``, maximum sum at most ``max_sum`` and
    maximum part at most ``max_part``. The integers ``w`` and ``s`` are
    respectively the remaining weight and sum.

    Be careful that the compositions ``c`` is modified inplace so that the
    output has to be copied before usage.

    TESTS::

        sage: from veerer.veering_quadrangulation import weighted_composition_iterator

        sage: [c[:] for c,_,_ in weighted_composition_iterator(1, [2], max_weight=6)] == [[1], [2], [3]]
        True
        sage: [c[:] for c,_,_ in weighted_composition_iterator(1, max_sum=3)] == [[1], [2], [3]]
        True
        sage: [c[:] for c,_,_ in weighted_composition_iterator(2, [1,3], max_weight=7, max_sum=3)] == [[1, 0], [2, 0], [3, 0], [0, 1], [1, 1], [2, 1], [0, 2], [1, 2]]
        True
    """
    if (weights is None) + (max_weight is None) == 1:
        raise ValueError("weights and max_weight should be specified together")

    if weights is None:
        if max_sum is None:
            raise ValueError("either weights or max_sum must be set")
        weights = [1] * n
        max_weight = max_sum
    elif len(weights) != n or any(w < 0 for w in weights):
        raise ValueError('invalid weights specification')

    if max_sum is None:
        if any(w == 0 for w in weights):
            raise ValueError('when zero weights are used, max_sum must be specified')
        max_sum = max_weight

    if max_part is None:
        max_part = max_sum

    c = [0] * n      # current composition
    w = max_weight   # available weight
    s = max_sum      # available sum

    while True:
        # find the first element available to increase, and set
        # all terms in between to zero
        i = 0
        while i < n and (s < 1 or w < weights[i] or c[i] == max_part):
            w += c[i] * weights[i]
            s += c[i]
            c[i] = 0
            i += 1

        if i == n:
            # can not increase any
            return
        c[i] += 1
        w -= weights[i]
        s -= 1

        yield c, w, s

def path_to_string(q, path):
    qq = q.copy()
    spath = []
    for step in path:
        sstep = ''
        for i,m in zip(*step):
            if m:
                sstep += '(' + ','.join(map(str,qq.r_orbit(i))) + ')^{}'.format(m)
                qq.r_staircase_move(i, m)
        qq.symmetry()
        spath.append(sstep)
    return ' | '.join(spath)

def all_paths(q,
        max_staircase_moves,     # upper bound on the number of staircase moves
        max_flips,               # upper bound on the number of diagonal changes
        max_blocks,              # upper bound on the number of blocks
        last_move=None):
    r"""
    Iterator through all reduced flip sequences starting from a given veering
    quadrangulation ``q`` and at most ``limit`` flips.

    INPUT::

    - ``q`` - initial quadrangulation

    - ``max_staircase_moves`` - bound on the number of staircase moves (or
      equivalently the number of flips)
      (sum)

    - ``max_flips`` - bound on the number of diagonal changes

    - ``max_blocks`` - bound on the number of blocks
      (alternation between left and right)

    - ``last_move`` - (optional) internal parameter used for the recursive calls
    """
    n = q._n
    if max_staircase_moves <= 0 or max_flips <= 0 or max_blocks <= 0:
        yield []
        return
    q0 = q.copy() # NOTE: just for extra check at the end of the loop
                  # could be removed in the future
    representatives, sizes = q.r_cycles_representatives(last_move)
    for multiplicities, new_max_flips, new_max_staircase_moves in weighted_composition_iterator(len(sizes), sizes, max_flips, max_staircase_moves):
        flipped = []

        for i,m in zip(representatives, multiplicities):
            if m:
                q.r_staircase_move(i, m)
                flipped.append(i)
        q.symmetry()

        prefix = [(representatives, multiplicities)]
        yield prefix
        for path in all_paths(q,
            max_staircase_moves=new_max_staircase_moves,
            max_flips=new_max_flips,
            max_blocks=max_blocks-1,
            last_move=flipped):
            yield prefix + path

        q.symmetry()
        for i,m in zip(representatives, multiplicities):
            if m:
                q.r_staircase_move(i, -m)
        assert q == q0

def all_paths_up_to_symmetry(q, limit, A, last_move=None):
    r"""
    Generate all paths starting from q up to the automorphism group A
    """
    # TODO

# TODO: in what form should we return the result?
# TODO: be smarter for the bound on eigenvalues
def hyperelliptic_length_spectrum(k, max_staircase_moves=None, max_flips=None, max_blocks=None, primitive=True, verbose=False):
    r"""
    Return the primitive hyperelliptic length spectrum in C^hyp(k).
    """
    k = int(k)

    if max_staircase_moves is None and max_flips is None:
        raise ValueError('at least one of max_staircase_moves or max_flips must be specified')
    elif max_staircase_moves is None:
        max_staircase_moves = max_flips
    elif max_flips is None:
        max_flips = k * max_staircase_moves

    if max_blocks is None:
        max_blocks = max_flips

    ans = {}
    for fp in all_reduced_and_complete_loops(k, max_staircase_moves, max_flips, max_blocks, verbose):
        cp = fp.matrix().charpoly()
        if cp not in ans:
            ans[cp] = []
        assert fp.num_blocks() <= fp.num_staircase_moves() <= fp.num_flips()
        ans[cp].append(fp)

    if not primitive:
        return ans

    # remove powers
    Rxy = QQ['x,y']
    x, y = Rxy.gens()
    Rx = QQ['x']
    Ry = QQ['y']
    for cp in sorted(ans, key = lambda x: min(fp.num_flips() for fp in ans[cp])):
        if not ans[cp]:
            # non-primitive data in the spectrum
            continue

        minblocks = min(fp.num_blocks() for fp in ans[cp])
        minstair = min(fp.num_staircase_moves() for fp in ans[cp])
        minflip = min(fp.num_flips() for fp in ans[cp])

        # ans[cp] should only be primitive elements, check that
        # each appear with the correct multiplicity
        invariants = {}
        for fp in ans[cp]:
            key = (fp.num_blocks(), fp.num_staircase_moves(), fp.num_flips())
            if key not in invariants:
                invariants[key] = 0
            invariants[key] += 1
        for (nblocks, nstairs, nflips), multiplicity in invariants.items():
            if multiplicity % nblocks != 0:
                print("WARNING: wrong multiplicity of (num_blocks, num_staircase_moves, num_flips) for charpoly={}".format(cp))
            elif multiplicity % nblocks != 0:
                print("We have an example where (char poly, num blocks, num staircase moves, num flips) does not determine the conjugacy class")

        # remove multiples
        k = 2
        while k * minflip < max_flips and k * minstair < max_staircase_moves and k * minblocks < max_blocks:
            # minimal polynomial of rho^k
            q = Rxy(cp).resultant(y - x**k, Rxy(x))
            cpk = Rx(Ry(q))
            assert cpk in ans
            for fp in ans[cp]:
                nblocks = fp.num_blocks()
                nstairs = fp.num_staircase_moves()
                nflips = fp.num_flips()
                if k * nblocks < max_blocks and k * nstairs < max_staircase_moves and k * nflips < max_flips:
                    ans[cpk].remove(fp**k)
            k += 1

    return ans

# TODO: marked loops, ie take into account possible automorphisms of the starting point
def all_reduced_and_complete_loops(start, max_staircase_moves, max_flips, max_blocks, verbose=False):
    r"""
    Each loop must appear as many times as the number of blocks.

    Currently only works for H(2g-2)^hyp (symmetries in the H(g-1,g-1)^hyp
    cases are not taken care of).
    """
    if isinstance(verbose, str):
        output = open(verbose, 'w')
    elif verbose:
        import sys
        output = sys.stdout

    # it is better to pick a Lyndon order that consider first the combinatorics and then
    # the multiplicities
    if isinstance(start, numbers.Integral):
        start = ferenczi_zamboni_class(start)
    else:
        raise TypeError('invalid input for q')


    for q in start:
        n = q._n
        qit = q.copy()
        A = qit.to_origami().automorphism_group()
        assert A.cardinality() == 1

        for path in all_paths(qit, max_staircase_moves, max_flips, max_blocks):
            ans, close = qit.is_isomorphic(q, True)
            if ans:
                fp = VeeringQuadrangulationFlipSequence.__new__(VeeringQuadrangulationFlipSequence)
                fp._start = q.copy()
                fp._end = q.copy()
                fp._flips = [(tuple(representatives), tuple(multiplicities)) for representatives, multiplicities in path]
                fp._relabelling = close
                fp._normalized = False
                fp._check()
                assert fp.is_loop()

                if verbose:
                    print("fp = {}".format(fp), file=output)
                    print("num blocks = {}".format(fp.num_blocks()), file=output)
                    print("num flips = {}".format(fp.num_flips()), file=output)
                    print("charpoly = {} = {}".format(fp.matrix().charpoly(), fp.matrix().charpoly().factor()), file=output)
                    print("closure = {}".format(fp._relabelling), file=output)
                    print("reduced at end = {}".format(fp.is_reduced_at_end()), file=output)
                    print("complete = {}".format(fp.is_complete()), file=output)

                if fp.is_reduced_at_end() and fp.is_complete():
                    yield fp

    if isinstance(verbose, str):
        output.close()

def one_quadrangulation(k):
    pr = perm_init([k-1] + list(range(k-1)))
    pl = perm_init(list(range(k-1,-1,-1)))
    return VeeringQuadrangulation(pr, pl)

def ferenczi_zamboni_class(k, labelled=False):
    r"""
    Ferenczi-Zamboni quadrangulations in hyperelliptic strata.

    EXAMPLES::

        sage: from veerer.veering_quadrangulation import ferenczi_zamboni_class
        sage: len(ferenczi_zamboni_class(2, labelled=True))  # optional - surface_dynamics
        3
        sage: len(ferenczi_zamboni_class(2, labelled=False))  # optional - surface_dynamics
        3

        sage: len(ferenczi_zamboni_class(3, labelled=True))  # optional - surface_dynamics
        9
        sage: len(ferenczi_zamboni_class(3, labelled=False))  # optional - surface_dynamics
        3

        sage: len(ferenczi_zamboni_class(4, labelled=True))  # optional - surface_dynamics
        28
        sage: len(ferenczi_zamboni_class(4, labelled=False))  # optional - surface_dynamics
        10
    """
    from surface_dynamics import AbelianStratum
    o = one_quadrangulation(k).to_origami()
    if k % 2:
        assert o.stratum_component() == AbelianStratum(k-1).hyperelliptic_component()
    elif k != 2:
        assert o.stratum_component() == AbelianStratum((k-2)//2, (k-2)//2).hyperelliptic_component()

    if not labelled:
        ocan = o.relabel(inplace=False)
        orbit = {ocan: o} # we keep the original as values of the dictionary
        todo = [ocan]
        while todo:
            ocan = todo.pop(0)
            o = orbit[ocan]
            for i in o.horizontal_cycle_representatives():
                oo = o.horizontal_twist(1, i)
                oocan = oo.relabel(inplace=False)
                if oocan not in orbit:
                    orbit[oocan] = oo
                    todo.append(oocan)
            for i in o.vertical_cycle_representatives():
                oo = o.vertical_twist(1, i)
                oocan = oo.relabel(inplace=False)
                if oocan not in orbit:
                    orbit[oocan] = oo
                    todo.append(oocan)
        orbit = orbit.values()

    else:
        todo = [o]
        orbit = set([o])
        while todo:
            o = todo.pop(0)
            for i in o.horizontal_cycle_representatives():
                oo = o.horizontal_twist(1, i)
                if oo not in orbit:
                    orbit.add(oo)
                    todo.append(oo)
            for i in o.vertical_cycle_representatives():
                oo = o.vertical_twist(1, i)
                if oo not in orbit:
                    orbit.add(oo)
                    todo.append(oo)

    return [VeeringQuadrangulation(o.r_tuple(), o.u_tuple()) for o in orbit]

# NOTE: assume no automorphism
def _uniq_primitive(paths):
    r"""
    Remove duplicates for paths with same dilatation.

    The implementation only discriminate paths by their characteristic
    polynomial, their number of blocks and number of flips.
    """
    nb_nf = {}
    for p in paths:
        k = (p.num_blocks(), p.num_flips())
        if k in nb_nf:
            nb_nf[k].append(p)
        else:
            nb_nf[k] = [p]

    # TODO: when the path is not symmetric but the start is, we get
    # multiple times the guy... (with different relabelling)

    for k, paths in nb_nf.items():
        nb, nf = k

        # check for single orbit
        if nb == len(paths):
            # what about non-primitive orbits!?
            yield min(paths)
        else:
            assert len(paths) % nb == 0
            classes = []
            for path in paths:
                is_new = True
                for c in classes:
                    if c[0].is_conjugate(path):
                        c.append(path)
                        is_new = False
                        break
                if is_new:
                    classes.append([path])
            for c in classes:
                assert len(c) == nb, (nb, [len(c) for c in classes])
                yield min(c)

def latex_array_small_dilatations(filename, k, limit, right_first=True, paper='a4paper', landscape=False, shift=False):
    from datetime import datetime
    dilatations = {}
    paths = {}
    for path in all_reduced_and_complete_loops(k, limit):
        cp = path.matrix().charpoly()
        if cp in dilatations:
            paths[cp].append(path)
        else:
            rho = max(cp.roots(RDF, False))
            if rho > 4:
                continue
            rho = max(cp.roots(AA, False))
            dilatations[cp] = rho
            paths[cp] = [path]

    lines = []
    rhomin = min(dilatations.values())
    for cp in sorted(paths, key=lambda x: dilatations[x]):
        rho = dilatations[cp]
        if dilatations[cp] == rhomin**2:
            break

        for path in _uniq_primitive(paths[cp]):
            nf = (2*path.num_flips()) if path.num_blocks() % 2 else path.num_flips()
            lines.append(path.latex(right_first=right_first, shift=shift) + ' & ' + str(rho.n(digits=6)) + ' & ' + latex(rho.minpoly()) + ' & ' + str(nf))

    with open(filename, 'w') as output:
        date = datetime.now()
        output.write('%' * 60)
        output.write('\n')
        output.write('%% Paths with small dilatations in Chyp(%d)\n' % k)
        output.write('%% limit = %d\n' % limit)
        output.write('%% computed on %s/%s/%s' % (date.year, date.month, date.day))
        output.write('%' * 60)
        output.write('\n')

        documentclass_options = [paper]
        if landscape:
            documentclass_options.append('landscape')
        output.write('\\documentclass[{}]{{article}}\n'.format(','.join(documentclass_options)))

        output.write('\\usepackage{tikz}\n')
        output.write('\\usepackage{amsmath}\n')
        output.write('\\usepackage[margin=2cm]{geometry}\n')

        output.write('\\begin{document}\n')
        output.write('\\[\n')
        output.write('\\def\\arraystretch{1.4}\n')
        output.write('\\begin{array}{lllllll}\n')
        if right_first:
            perm = '\\pi_r & \\pi_\\ell'
        else:
            perm = '\\pi_\\ell & \\pi_r'
        output.write('{} & \\gamma & p & \\text{{dilatation}} & \\text{{minimal polynomial}} & \\text{{diag. changes}} \\\\\\hline\n'.format(perm))
        output.write(' \\\\\n'.join(lines))
        output.write('\n')
        output.write('\\end{array}\n')
        output.write('\\]\n')
        output.write('\\end{document}')

    return len(lines)
