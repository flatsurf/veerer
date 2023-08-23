r"""
Linear expression and linear constraints

This module provides a common interface to ppl, PyNormaliz and sage to build
polyhedra.
"""
from sage.categories.modules import Modules

from sage.structure.element import ModuleElement, parent
from sage.structure.parent import Parent
from sage.structure.richcmp import op_LE, op_LT, op_EQ, op_NE, op_GT, op_GE

from sage.arith.functions import lcm


class LinearExpression(ModuleElement):
    r"""
    EXAMPLES::

        sage: from veerer.linear_expression import LinearExpressions
        sage: L = LinearExpressions(QQ)
        sage: L()
        0
        sage: L({0: 5})
        5 * x[0]
        sage: L({0: 5, 1: -2})
        5 * x[0] + -2 * x[1]
    """
    def __init__(self, parent, f=None, inhomogeneous_term=None):
        ModuleElement.__init__(self, parent)
        self._f = {} if not f else f
        base_ring = parent.base_ring()
        self._inhomogeneous_term = base_ring.zero() if not inhomogeneous_term else inhomogeneous_term

    def is_homogeneous(self):
        return self._inhomogeneous_term.is_zero()

    def _repr_(self):
        if not self._f:
            return str(self._inhomogeneous_term)

        # TODO:nicer string representation with + and -
        lin_part = ' + '.join('%s * x[%d]' % (self._f[i], i) for i in sorted(self._f))
        if self._inhomogeneous_term:
            return lin_part + ' + ' + str(self._inhomogeneous_term)
        else:
            return lin_part

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
        return lcm(list(coeff.denominator() for coeff in self._f.values()) + [self._inhomogeneous_term.denominator()])

    def ppl(self):
        r"""
        EXAMPLES::

            sage: from veerer.linear_expression import LinearExpressions
            sage: L = LinearExpressions(QQ)
            sage: (3 * L.variable(2) - 7/2 * L.variable(5) + 1/3).ppl()
            18*x2-21*x5+2
        """
        from gmpy2 import mpz
        import ppl
        # need to chase denominator first
        den = self.denominator()
        if den.is_one():
            lin = self
        else:
            lin = den * self
        return sum(mpz(coeff) * ppl.Variable(i) for i, coeff in lin._f.items()) + mpz(lin._inhomogeneous_term)

    def coefficients(self, dim=None):
        r"""
        EXAMPLES::

            sage: from veerer.linear_expression import LinearExpressions
            sage: L = LinearExpressions(QQ)
            sage: (3 * L.variable(2) - 7/2 * L.variable(5) + 1/3).vector()
            sage: (3 * L.variable(2) - 7/2 * L.variable(5) + 1/3).vector(10)
        """
        if dim is None:
            if not self._f:
                dim = 0
            else:
                dim = max(self._f) + 1
        zero = self.base_ring().zero()
        return [-self._inhomogeneous_term] + [self._f.get(i, zero) for i in range(dim)]


class LinearConstraint:
    r"""
    EXAMPLES::

        sage: from veerer.linear_expression import LinearExpressions
        sage: L = LinearExpressions(QQ)
        sage: 3 * L.variable(0) - 5 * L.variable(2) == 7
    """
    def __init__(self, op, left, right):
        self._expression = left - right
        self._op = op

    def __repr__(self):
        if self._op == op_LE:
            op = '<='
        elif self._op == op_LT:
            op = '<'
        elif self._op == op_EQ:
            op = '='
        elif self._op == op_NE:
            op = '!='
        elif self._op == op_GT:
            op = '>'
        elif self._op == op_GE:
            op = '>='
        else:
            raise RuntimeError
        return '{} {} 0'.format(self._expression, op)

    def ppl(self):
        r"""
        EXAMPLES::

            sage: from veerer.linear_expression import LinearExpressions
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


class ConstraintSystem:
    r"""
    EXAMPLES::

        sage: from veerer.linear_expression import *
        sage: L = LinearExpressions(QQ)
        sage: cs = ConstraintSystem()
        sage: cs.insert(L.variable(0) >= 2)
        sage: cs.insert(2 * L.variable(1) - L.variable(3) <= 18)
        sage: cs
        {1 * x[0] + -2 >= 0, 2 * x[1] + -1 * x[3] + -18 <= 0}
    """
    def __init__(self):
        self._data = []
        self._dim = None

    def __repr__(self):
        return '{' + ', '.join(map(str, self)) + '}'

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

            sage: from veerer.linear_expression import LinearExpressions, ConstraintSystem
            sage: L = LinearExpressions(QQ)
            sage: cs = ConstraintSystem()
            sage: cs.insert(2 * L.variable(0) - 3/5 * L.variable(2) >= 1)
            sage: cs.insert(L.variable(0) + L.variable(1) + L.variable(2) == 3)
            sage: cs.ppl()
            Constraint_System {10*x0-3*x2-5>=0, x0+x1+x2-3==0}
        """
        import ppl
        cs = ppl.Constraint_System()
        for constraint in self._data:
            cs.insert(constraint.ppl())
        return cs

    def ieqs_eqns(self, dim=None):
        r"""
        EXAMPLES::

            sage: from veerer.linear_expression import LinearExpressions, ConstraintSystem
            sage: L = LinearExpressions(QQ)
            sage: cs = ConstraintSystem()
            sage: cs.insert(2 * L.variable(0) - 3/5 * L.variable(2) >= 1)
            sage: cs.insert(L.variable(3) <= 1)
            sage: cs.insert(L.variable(0) + L.variable(1) + L.variable(2) == 3)
            sage: cs.ieqs_eqns()
            ([[1, 2, 0, -3/5, 0], [-1, 0, 0, 0, -1]], [[3, 1, 1, 1, 0]])

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
                ieqs.append(constraint._expression.coefficients(dim))
            elif constraint._op == op_LE:
                ieqs.append((-constraint._expression).coefficients(dim))
            elif constraint._op == op_EQ:
                eqns.append(constraint._expression.coefficients(dim))
            else:
                raise ValueError('invalid constraint {}'.format(constraint))
        return ieqs, eqns

    def is_linear_subspace(self):
        return all(constraint._op == op_EQ and constraint._expression._inhomogeneous_term.is_zero() for constraint in self)

    def linear_generators_matrix(self, dim=None):
        from sage.matrix.constructor import matrix
        if dim is None:
            dim = self._dim
        if not self.is_linear_subspace():
            raise ValueError('not a linear subspace')
        mat = matrix([constraint._expression.coefficients(dim)[1:] for constraint in self])
        return mat.right_kernel_matrix()


class LinearExpressions(Parent):
    r"""
    EXAMPLES::

        sage: from veerer.linear_expression import LinearExpressions
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

        sage: from veerer.linear_expression import *
        sage: L = LinearExpressions(QQ)
        sage: L.variable(0)
        1 * x[0]
        sage: L.variable(1)
        1 * x[1]
        sage: 5 * L.variable(2) - 3.variable(7)
        5 * x[2] + -3 * x[7]
        """
        return LinearExpression(self, {i: self.base_ring().one()})

    def _element_constructor_(self, *args):
        if not args:
            return LinearExpression(self)
        elif len(args) == 1:
            data = args[0]
            if isinstance(data, dict):
                base_ring = self.base_ring()
                data = {int(i): base_ring(coeff) for i, coeff in data.items()}
            elif data in self.base_ring():
                return LinearExpression(self, inhomogeneous_term=data)
            else:
                raise ValueError('can not construct linear expression from {}'.format(data))
        elif len(args) == 2:
            data0 = args[0]
            data1 = args[1]
            if not isinstance(data0, dict) or not data1 in self.base_ring():
                raise ValueError('can not construct linear expression from {}'.format(args))
                base_ring = self.base_ring()
            data0 = {int(i): base_ring(coeff) for i, coeff in data0.items()}
            return LinearExpression(self, data0, data1)
        else:
            raise ValueError('can not construct linear expression from {}'.format(args))


def polyhedron(cs, backend):
    if backend == 'ppl':
        import ppl
        return ppl.C_Polyhedron(cs.ppl())
    elif backend == 'sage':
        from sage.geometry.polyhedron.constructor import Polyhedron
        ieqs, eqns = cs.ieqs_eqns()
        return Polyhedron(ieqs=ieqs, eqns=eqns)
    else:
        raise RuntimeError

def polyhedron_add_constraints(P, cs, backend):
    if backend == 'ppl':
        import ppl
        P = ppl.C_Polyhedron(P)
        if isinstance(cs, LinearConstraint):
            P.add_constraint(cs.ppl())
        elif isinstance(cs, ConstraintSystem):
            for constraint in cs:
                P.add_constraint(constraint.ppl())
        else:
            raise TypeError
        return P
    elif backend == 'sage':
        from sage.geometry.polyhedron.constructor import Polyhedron
        if isinstance(cs, LinearConstraint):
            constraint = cs
            cs = ConstraintSystem()
            cs.insert(constraint)
        elif isinstance(cs, ConstraintSystem):
            pass
        else:
            raise TypeError
        ieqs, eqns = cs.ieqs_eqns(P.ambient_dimension())
        return P.intersection(Polyhedron(ieqs=ieqs, eqns=eqns))

def polyhedron_dimension(P, backend):
    if backend == 'ppl':
        return P.affine_dimension()
    elif backend == 'sage':
        return P.dimension()
    else:
        raise TypeError

def polyhedron_to_hashable(P, backend):
    r"""
    EXAMPLES::

        sage: from veerer.veering_triangulation import VeeringTriangulation, ppl_cone_to_hashable, ppl_cone_from_hashable
        sage: from surface_dynamics import *  # optional - surface_dynamics

        sage: T = VeeringTriangulation.from_stratum(AbelianStratum(2))   # optional - surface_dynamics
        sage: P = T.geometric_polytope()      # optional - surface_dynamics
        sage: h = ppl_cone_to_hashable(P)     # optional - surface_dynamics
        sage: P == ppl_cone_from_hashable(h)  # optional - surface_dynamics
        True
    """
    if backend == 'ppl':
        eqns = []
        ieqs = []
        for constraint in P.minimized_constraints():
            if constraint.inhomogeneous_term():
                raise ValueError('not a cone')
            if constraint.is_equality():
                eqns.append(tuple(constraint.coefficients()))
            elif constraint.is_inequality():
                ieqs.append(tuple(constraint.coefficients()))
            else:
                raise RuntimeError
        eqns.sort()
        ieqs.sort()
        return (P.space_dimension(), tuple(eqns), tuple(ieqs))
    elif backend == 'sage':
        if P.is_mutable():
            P = P.__copy__()
            P.set_immutable()
            return P
        else:
            return P
    else:
        raise TypeError


def polyhedron_from_hashable(args, backend):
    if backend == 'ppl':
        d, eqns, ieqs = args
        P = ppl.C_Polyhedron(d)
        for constraint in eqns:
            P.add_constraint(sum(coeff * ppl.Variable(i) for i,coeff in enumerate(constraint)) == 0)
        for constraint in ieqs:
            P.add_constraint(sum(coeff * ppl.Variable(i) for i,coeff in enumerate(constraint)) >= 0)
        return P
    else:
        return args



