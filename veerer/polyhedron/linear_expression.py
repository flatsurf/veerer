r"""
Linear expression and linear constraints

This module provides a common interface to ppl, PyNormaliz and sage to build
cones and polyhedra.
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2023 Vincent Delecroix
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

import ppl

from sage.categories.modules import Modules
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import ModuleElement, Vector, parent, get_coercion_model
from sage.structure.parent import Parent
from sage.structure.richcmp import op_LE, op_LT, op_EQ, op_NE, op_GT, op_GE
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.functions import lcm
from sage.arith.functions import LCM_list


cm = get_coercion_model()


# TODO: allow base_ring=int or Fraction so that we can handle pure
# Python things
class LinearExpression(ModuleElement):
    r"""
    EXAMPLES::

        sage: from veerer.polyhedron import LinearExpressions  # random output due to deprecation warnings in realalg
        sage: L = LinearExpressions(QQ)
        sage: L()
        0
        sage: L({0: 5})
        5*x0
        sage: L({0: 5, 1: -2})
        5*x0 - 2*x1
    """
    def __init__(self, parent, f=None, inhomogeneous_term=None):
        ModuleElement.__init__(self, parent)
        self._f = {} if not f else f
        base_ring = parent.base_ring()
        self._inhomogeneous_term = base_ring.zero() if not inhomogeneous_term else inhomogeneous_term

    def change_ring(self, base_ring):
        r"""
        Return the same linear expression on a different base ring.

        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions
            sage: L = LinearExpressions(QQ)
            sage: L({0: 5}, 3)
            5*x0 + 3
        """
        if base_ring is self.base_ring():
            return self
        parent = LinearExpressions(base_ring)
        f = {i: base_ring(coeff) for i, coeff in self._f.items()}
        inhomogeneous_term = base_ring(self._inhomogeneous_term)
        return parent.element_class(parent, f, inhomogeneous_term)

    def is_homogeneous(self):
        return self._inhomogeneous_term.is_zero()

    def _repr_(self):
        if not self._f:
            return str(self._inhomogeneous_term)

        def term_string(i, s):
            if s == '1':
                return 'x%d' % i
            else:
                return '%s*x%d' % (s, i)

        # TODO:nicer string representation with + and -
        ind_coeffs = [(i, str(self._f[i])) for i in sorted(self._f)]
        lin_part = term_string(*ind_coeffs[0])
        for i, s in ind_coeffs[1:]:
            if s[0] == '-':
                lin_part += ' - '
                s = s[1:]
            else:
                lin_part += ' + '
            lin_part += term_string(i, s)
        if self._inhomogeneous_term:
            s = str(self._inhomogeneous_term)
            if s[0] == '-':
                return lin_part + ' - ' + s[1:]
            else:
                return lin_part + ' + ' + s
        else:
            return lin_part

    def __call__(self, x=None):
        r"""
        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions
            sage: L = LinearExpressions(QQ)
            sage: l = L({0: 5, 7: -1}, 3)
            sage: l
            5*x0 - x7 + 3
            sage: l([])
            3
            sage: l([1, 2, 3, 4, 5])
            8
            sage: l([1, 2, 3, 4, 5, 6, 7, 8, 9])
            0

            sage: l({0: 1, 1: 1, 7: 5})
            3
        """
        if not x or not self._f:
            return self._inhomogeneous_term
        if isinstance(x, dict):
            lin = sum(coeff * x[i] for i, coeff in self._f.items() if i in x)
            return lin + self._inhomogeneous_term
        elif isinstance(x, (tuple, list, Vector)):
            lin = sum(coeff * x[i] for i, coeff in self._f.items() if i < len(x))
            return lin + self._inhomogeneous_term
        else:
            raise TypeError('invalid input')

    def _lmul_(self, other):
        f = {i: other * coeff for i, coeff in self._f.items()}
        inhomogeneous_term = other * self._inhomogeneous_term
        return self.parent().element_class(self.parent(), f, inhomogeneous_term)

    def _rmul_(self, other):
        f = {i: coeff * other for i, coeff in self._f.items()}
        inhomogeneous_term = self._inhomogeneous_term * other
        return self.parent().element_class(self.parent(), f, inhomogeneous_term)

    def _add_(self, other):
        data = self._f.copy()
        for i, coeff in other._f.items():
            if i in data:
                data[i] += coeff
                if not data[i]:
                    del data[i]
            else:
                data[i] = coeff
        inhomogeneous_term = self._inhomogeneous_term + other._inhomogeneous_term
        return self.parent().element_class(self.parent(), data, inhomogeneous_term)

    def _sub_(self, other):
        data = self._f.copy()
        for i, coeff in other._f.items():
            if i in data:
                data[i] -= coeff
                if not data[i]:
                    del data[i]
                else:
                    data[i] = -coeff
            else:
                data[i] = -coeff
        inhomogeneous_term = self._inhomogeneous_term - other._inhomogeneous_term
        return self.parent().element_class(self.parent(), data, inhomogeneous_term)

    def __neg__(self):
        return LinearExpression(self.parent(), {i: -c for (i, c) in self._f.items()}, -self._inhomogeneous_term)

    def __le__(self, other):
        return LinearConstraint(op_LE, self, other)

    def __lt__(self, other):
        return LinearConstraint(op_LT, self, other)

    def __eq__(self, other):
        return LinearConstraint(op_EQ, self, other)

    def __ne__(self, other):
        return LinearConstraint(op_NE, self, other)

    def __gt__(self, other):
        return LinearConstraint(op_GT, self, other)

    def __ge__(self, other):
        return LinearConstraint(op_GE, self, other)

    def _richcmp_(self, other, op):
        return LinearExpression(op, self, other)

    def denominator(self):
        r"""
        Return the common denominator of the coefficients of this linear expression.
        """
        return LCM_list([coeff.denominator() for coeff in self._f.values()] + [self._inhomogeneous_term.denominator()])

    # TODO: try as much as possible to avoid calls to integral
    def integral(self):
        if self.base_ring() is ZZ:
            return self
        if self.base_ring() is QQ:
            return (self.denominator() * self).change_ring(ZZ)
        raise ValueError('invalid base ring')

    def ppl(self):
        r"""
        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions
            sage: L = LinearExpressions(QQ)
            sage: (3 * L.variable(2) - 7/2 * L.variable(5) + 1/3).ppl()
            18*x2-21*x5+2
        """
        from gmpy2 import mpz
        lin = self.integral()
        # TODO: the line below is too costly : it accounts for 80% of the
        # geometric automaton computation
        return sum(mpz(coeff) * ppl.Variable(i) for i, coeff in lin._f.items()) + mpz(lin._inhomogeneous_term)

    def coefficients(self, dim=None, homogeneous=False):
        r"""
        Return the coefficients of this linear form as a list.

        If ``homogeneous`` is ``False`` (default) the first coefficient is the homogeneous term.

        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions
            sage: L = LinearExpressions(QQ)
            sage: lin = 3 * L.variable(2) - 7/2 * L.variable(5) + 1/3
            sage: lin.coefficients()
            [-1/3, 0, 0, 3, 0, 0, -7/2]
            sage: lin.coefficients(homogeneous=True)
            [0, 0, 3, 0, 0, -7/2]

            sage: lin.coefficients(10)
            [-1/3, 0, 0, 3, 0, 0, -7/2, 0, 0, 0, 0]
        """
        if dim is None:
            if not self._f:
                dim = 0
            else:
                dim = max(self._f) + 1
        zero = self.base_ring().zero()
        if homogeneous:
            return [self._f.get(i, zero) for i in range(dim)]
        else:
            return [-self._inhomogeneous_term] + [self._f.get(i, zero) for i in range(dim)]


class LinearConstraint:
    r"""
    EXAMPLES::

        sage: from veerer.polyhedron import LinearExpressions
        sage: L = LinearExpressions(QQ)
        sage: 3 * L.variable(0) - 5 * L.variable(2) == 7
        3*x0 - 5*x2 - 7 == 0
    """
    def __init__(self, op, left, right):
        self._expression = left - right
        self._op = op

    @staticmethod
    def _from_ppl(cst):
        f = cst.coefficients()
        inhomogeneous_term = cst.inhomogeneous_term()
        L = LinearExpressions(ZZ)
        lin = L.element_class(L, {i: c for i, c in enumerate(f) if c}, inhomogeneous_term)
        if cst.is_equality():
            return LinearConstraint(op_EQ, lin, ZZ.zero())
        elif cst.is_inequality():
            return LinearConstraint(op_GE, lin, ZZ.zero())
        raise ValueError('invalid constraint cst={}'.format(cst))

    def __repr__(self):
        if self._op == op_LE:
            op = '<='
        elif self._op == op_LT:
            op = '<'
        elif self._op == op_EQ:
            op = '=='
        elif self._op == op_NE:
            op = '!='
        elif self._op == op_GT:
            op = '>'
        elif self._op == op_GE:
            op = '>='
        else:
            raise RuntimeError
        return '{} {} 0'.format(self._expression, op)

    def __call__(self, x=None):
        r"""
        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions
            sage: L = LinearExpressions(QQ)
            sage: x0 = L.variable(0)
            sage: x1 = L.variable(1)
            sage: x2 = L.variable(2)
            sage: constraint = 3*x0 - 2*x1 + 7*x2 >= 5
            sage: constraint([-4, 7, 3])
            False
            sage: constraint([0, 1, 1])
            True

            sage: constraint({0:1, 1:1, 2:1})
            True
        """
        ev = self._expression(x)
        if self._op == op_LE:
            return ev <= 0
        elif self._op == op_LT:
            return ev < 0
        elif self._op == op_EQ:
            return ev == 0
        elif self._op == op_NE:
            return ev != 0
        elif self._op == op_GT:
            return ev > 0
        elif self._op == op_GE:
            return ev >= 0
        else:
            raise RuntimeError

    def is_homogeneous(self):
        return self._expression.is_homogeneous()

    def coefficients(self, dim=None, homogeneous=False):
        return self._expression.coefficients(dim, homogeneous)

    def ppl(self):
        r"""
        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions
            sage: L = LinearExpressions(QQ)
            sage: (3 * L.variable(0) >= 1).ppl()
            3*x0-1>=0
            sage: (3 * L.variable(0) > 1).ppl()
            3*x0-1>0
            sage: (3 * L.variable(0) == 1).ppl()
            3*x0-1==0
            sage: (3 * L.variable(0) < 1).ppl()
            -3*x0+1>0
            sage: (3 * L.variable(0) <= 1).ppl()
            -3*x0+1>=0
        """
        if self._op == op_LE:
            return self._expression.ppl() <= 0
        elif self._op == op_LT:
            return self._expression.ppl() < 0
        elif self._op == op_EQ:
            return self._expression.ppl() == 0
        elif self._op == op_NE:
            return self._expression.ppl() != 0
        elif self._op == op_GT:
            return self._expression.ppl() > 0
        elif self._op == op_GE:
            return self._expression.ppl() >= 0

    def integral(self):
        return LinearConstraint(self._op, self._expression.integral(), ZZ.zero())


class ConstraintSystem:
    r"""
    EXAMPLES::

        sage: from veerer.polyhedron import *
        sage: L = LinearExpressions(QQ)
        sage: cs = ConstraintSystem()
        sage: cs.insert(L.variable(0) >= 2)
        sage: cs.insert(2 * L.variable(1) - L.variable(3) <= 18)
        sage: cs
        {x0 - 2 >= 0, 2*x1 - x3 - 18 <= 0}
    """
    def __init__(self):
        self._data = []
        self._dim = None

    def base_ring(self):
        r"""
        EXAMPLES::

            sage: from veerer.polyhedron import *
            sage: L = LinearExpressions(QQ)
            sage: cs = ConstraintSystem()
            sage: cs.insert(L.variable(0) >= 2)
            sage: cs.insert(2 * L.variable(1) - L.variable(3) <= 18)
            sage: cs.base_ring()
            Rational Field
        """
        if not self._data:
            raise ValueError
        return cm.common_parent(*[constraint._expression.parent().base_ring() for constraint in self._data])

    def __repr__(self):
        return '{' + ', '.join(map(str, self)) + '}'

    def __call__(self, x=None):
        return [constraint(x) for constraint in self._data]

    def copy(self):
        cs = ConstraintSystem.__new__(ConstraintSystem)
        cs._data = self._data[:]
        cs._dim = self._dim
        return cs

    def insert(self, constraint):
        if not isinstance(constraint, LinearConstraint):
            raise TypeError('invalid type; expected LinearConstraint but got {}'.format(type(constraint).__name__))
        if self._dim is None:
            self._dim = max(constraint._expression._f) + 1
        else:
            self._dim = max(self._dim, max(constraint._expression._f) + 1)
        self._data.append(constraint)

    def __iter__(self):
        return iter(self._data)

    def ppl(self):
        r"""
        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions, ConstraintSystem
            sage: L = LinearExpressions(QQ)
            sage: cs = ConstraintSystem()
            sage: cs.insert(2 * L.variable(0) - 3/5 * L.variable(2) >= 1)
            sage: cs.insert(L.variable(0) + L.variable(1) + L.variable(2) == 3)
            sage: cs.ppl()
            Constraint_System {10*x0-3*x2-5>=0, x0+x1+x2-3==0}
        """
        cs = ppl.Constraint_System()
        for constraint in self._data:
            cs.insert(constraint.ppl())
        return cs

    def ieqs_eqns(self, dim=None, homogeneous=False):
        r"""
        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions, ConstraintSystem
            sage: L = LinearExpressions(QQ)
            sage: cs = ConstraintSystem()
            sage: cs.insert(2 * L.variable(0) - 3/5 * L.variable(2) >= 1)
            sage: cs.insert(L.variable(3) <= 1)
            sage: cs.insert(L.variable(0) + L.variable(1) + L.variable(2) == 3)
            sage: cs.ieqs_eqns()
            ([[1, 2, 0, -3/5, 0], [-1, 0, 0, 0, -1]], [[3, 1, 1, 1, 0]])
            sage: cs.ieqs_eqns(homogeneous=True)
            ([[2, 0, -3/5, 0], [0, 0, 0, -1]], [[1, 1, 1, 0]])

            sage: ieqs, eqns = cs.ieqs_eqns()
            sage: Polyhedron(ieqs=ieqs, eqns=eqns)
            A 3-dimensional polyhedron in QQ^4 defined as the convex hull of 1 vertex, 2 rays, 1 line
        """
        if dim is None:
            dim = self._dim
        ieqs = []
        eqns = []
        for constraint in self._data:
            if constraint._op == op_GE:
                ieqs.append(constraint._expression.coefficients(dim, homogeneous))
            elif constraint._op == op_LE:
                ieqs.append((-constraint._expression).coefficients(dim, homogeneous))
            elif constraint._op == op_EQ:
                eqns.append(constraint._expression.coefficients(dim, homogeneous))
            else:
                raise ValueError('invalid constraint {}'.format(constraint))
        return ieqs, eqns

    def is_linear_subspace(self):
        r"""
        Return whether all constraints are equalities with vanishing inhomogeneous term.
        """
        return all(constraint._op == op_EQ and constraint._expression._inhomogeneous_term.is_zero() for constraint in self)

    def linear_generators_matrix(self, dim=None):
        from sage.matrix.constructor import matrix
        if dim is None:
            dim = self._dim
        if not self.is_linear_subspace():
            raise ValueError('not a linear subspace')
        mat = matrix([constraint._expression.coefficients(dim)[1:] for constraint in self])
        return mat.right_kernel_matrix()

    def integral(self):
        cs = ConstraintSystem()
        for constraint in self:
            cs.insert(constraint.integral())
        return cs

    def cone(self, backend='sage'):
        r"""
        Return the cone defined from these constraints.

        INPUT:

        - ``backend`` (optional): either ``'ppl'``, ``'sage'``, ``'normaliz-QQ'`` or ``'noramliz-NF'``
        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions, ConstraintSystem

            sage: L = LinearExpressions(QQ)
            sage: cs = ConstraintSystem()
            sage: cs.insert(L.variable(0) >= 0)
            sage: cs.insert(L.variable(1) >= 0)
            sage: cs.insert(L.variable(2) >= 0)
            sage: cs.insert(L.variable(3) >= 0)
            sage: cs.insert(L.variable(0) + 1/5 * L.variable(1) - 2/3 * L.variable(2) + 6 * L.variable(3) == 0)

            sage: cs.cone('sage')
            Cone of dimension 3 in ambient dimension 4 made of 3 facets (backend=sage)
            sage: cs.cone('ppl')
            Cone of dimension 3 in ambient dimension 4 made of 3 facets (backend=ppl)
            sage: cs.cone('normaliz-QQ')  # optional - pynormaliz
            Cone of dimension 3 in ambient dimension 4 made of 3 facets (backend=normaliz-QQ)

        An example over a number field::

            sage: K = NumberField(x^2 - x - 1, 'phi', embedding=(AA(5).sqrt() + 1)/2)
            sage: phi = K.gen()
            sage: L = LinearExpressions(K)
            sage: cs = ConstraintSystem()
            sage: cs.insert(L.variable(0) - phi * L.variable(1) >= 0)
            sage: cs.insert(L.variable(0) + L.variable(1) + L.variable(2) == 0)
            sage: cs.cone('sage')
            Cone of dimension 2 in ambient dimension 3 made of 1 facets (backend=sage)
            sage: cs.cone('normaliz-NF')  # optional - pynormaliz
            Cone of dimension 2 in ambient dimension 3 made of 1 facets (backend=normaliz-NF)
        """
        # homogeneous case
        if not all(constraint.is_homogeneous() for constraint in self):
            raise ValueError

        if backend is None:
            from .cone import default_backend
            backend = default_backend(self.base_ring())

        if backend == 'ppl':
            from .cone import Cone_ppl
            return Cone_ppl(QQ, ppl.C_Polyhedron(self.ppl()))
        elif backend == 'sage':
            from .cone import Cone_sage
            ieqs, eqns = self.ieqs_eqns()
            return Cone_sage._new(ieqs, eqns)
        elif backend == 'normaliz-QQ':
            from ..features import pynormaliz_feature
            pynormaliz_feature.require()

            from .cone import Cone_normalizQQ
            ieqs, eqns = self.ieqs_eqns(homogeneous=True)
            return Cone_normalizQQ._new(ieqs, eqns)
        elif backend == 'normaliz-NF':
            from ..features import pynormaliz_feature
            pynormaliz_feature.require()

            from .cone import Cone_normalizNF
            ieqs, eqns = self.ieqs_eqns(homogeneous=True)
            return Cone_normalizNF._new(ieqs, eqns)
        else:
            raise ValueError('invalid backend')


class LinearExpressions(UniqueRepresentation, Parent):
    r"""
    Linear expressions over a given base ring.

    EXAMPLES::

        sage: from veerer.polyhedron import LinearExpressions
        sage: LinearExpressions(QQ)
        Linear expressions over Rational Field
        sage: LinearExpressions(QuadraticField(5))
        Linear expressions over Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
    """
    Element = LinearExpression

    def __init__(self, base_ring):
        Parent.__init__(self, base=base_ring, category=Modules(base_ring))
        self._populate_coercion_lists_(coerce_list=[base_ring])

    def _repr_(self):
        return 'Linear expressions over {}'.format(self.base_ring())

    def variable(self, i):
        r"""
        EXAMPLES::

        sage: from veerer.polyhedron import LinearExpressions
        sage: L = LinearExpressions(QQ)
        sage: L.variable(0)
        x0
        sage: L.variable(1)
        x1
        sage: 5 * L.variable(2) - 3 * L.variable(7)
        5*x2 - 3*x7
        """
        return LinearExpression(self, {int(i): self.base_ring().one()})

    def _element_constructor_(self, *args):
        if not args:
            return LinearExpression(self)
        elif len(args) == 1:
            base_ring = self.base_ring()
            data = args[0]
            if isinstance(data, (tuple, list, Vector)):
                # TODO: should we consider vector as homogeneous or inhomogeneous?
                f = {i: base_ring(coeff) for i, coeff in enumerate(data) if coeff}
                inhomogeneous_term = base_ring.zero()
            elif isinstance(data, dict):
                f = {int(i): base_ring(coeff) for i, coeff in data.items() if coeff}
                inhomogeneous_term = base_ring.zero()
            elif data in base_ring:
                f = {}
                inhomogeneous_term = base_ring(data)
            else:
                raise ValueError('can not construct linear expression from {}'.format(data))
            return LinearExpression(self, f, inhomogeneous_term)

        elif len(args) == 2:
            data0 = args[0]
            data1 = args[1]
            base_ring = self.base_ring()

            if isinstance(data0, (tuple, list, Vector)):
                data0 = {i: base_ring(coeff) for i, coeff in enumerate(data) if coeff}
            elif isinstance(data0, dict):
                data0 = {int(i): base_ring(coeff) for i, coeff in data0.items() if coeff}
            else:
                raise ValueError('invalid first argument {}'.format(data0))
            data1 = base_ring(data1)
            return LinearExpression(self, data0, data1)
        else:
            raise ValueError('can not construct linear expression from {}'.format(args))
