r"""
Cone

A cone in R^d is a domain delimited by finitely linear homogeneous inequalities
`a_1 x_1 + a_2 x_2 + \ldots + a_d x_d \geq 0`.

EXAMPLES::

    sage: from veerer.polyhedron.linear_expression import LinearExpressions, ConstraintSystem  # random output due to deprecation warnings in realalg

    sage: L = LinearExpressions(QQ)
    sage: cs = ConstraintSystem()
    sage: cs.insert(L.variable(0) >= 0)
    sage: cs.insert(L.variable(1) >= 0)
    sage: cs.insert(L.variable(2) >= 0)
    sage: cs.insert(L.variable(3) >= 0)
    sage: cs.insert(L.variable(0) + 1/5 * L.variable(1) - 2/3 * L.variable(2) + 6 * L.variable(3) == 0)

    sage: cones = []
    sage: cone_sage = cs.cone('sage')
    sage: cones.append(cone_sage)
    sage: cone_ppl = cs.cone('ppl')
    sage: cones.append(cone_ppl)
    sage: cone_nmz = cs.cone('normaliz-QQ')  # optional - pynormaliz
    sage: cones.append(cone_nmz)  # optional - pynormaliz
    sage: assert all(cone.space_dimension() == 4 for cone in cones)
    sage: assert all(cone.affine_dimension() == 3 for cone in cones)
    sage: assert all(sorted(cone.ieqs()) == [[-15, -3, 10, 0], [0, 1, 0, 0], [1, 0, 0, 0]] for cone in cones)
    sage: assert all(sorted(cone.eqns()) == [[15, 3, -10, 90]] for cone in cones)
    sage: assert all(sorted(cone.rays()) == [[0, 0, 9, 1], [0, 10, 3, 0], [2, 0, 3, 0]] for cone in cones)

    sage: cone_ppl.add_constraint(L.variable(0) == L.variable(3))
    Cone of dimension 2 in ambient dimension 4 made of 2 facets (backend=ppl)
    sage: cone_sage.add_constraint(L.variable(0) == L.variable(3))
    Cone of dimension 2 in ambient dimension 4 made of 2 facets (backend=sage)
    sage: cone_nmz.add_constraint(L.variable(0) == L.variable(3)) # optional - pynormaliz
    Cone of dimension 2 in ambient dimension 4 made of 2 facets (backend=normaliz-QQ)

An example over a number field::

    sage: K = NumberField(x^2 - x - 1, 'phi', embedding=(AA(5).sqrt() + 1)/2)
    sage: phi = K.gen()
    sage: L = LinearExpressions(K)
    sage: cs = ConstraintSystem()
    sage: cs.insert(L.variable(0) - phi * L.variable(1) >= 0)
    sage: cs.insert(L.variable(0) + L.variable(1) + L.variable(2) == 0)
    sage: cone_sage = cs.cone('sage')
    sage: cone_nmz = cs.cone('normaliz-NF') # optional - pynormaliz

Note that sage and noramliz use different noramlization::

    sage: cone_sage.ieqs()
    [[1, -phi, 0]]
    sage: cone_nmz.ieqs() # optional - pynormaliz
    [[phi - 1, -1, 0]]

    sage: cone_sage.eqns()
    [[1, 1, 1]]
    sage: cone_nmz.eqns() # optional - pynormaliz
    [[1, 1, 1]]

Testing adding constraints::

    sage: cone_sage.add_constraint(phi * L.variable(0) - 3 * L.variable(2) <= 0)
    Cone of dimension 2 in ambient dimension 3 made of 2 facets (backend=sage)
    sage: cone_nmz.add_constraint(phi * L.variable(0) - 3 * L.variable(2) <= 0) # optional - pynoramliz
    Cone of dimension 2 in ambient dimension 3 made of 2 facets (backend=normaliz-NF)
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

from fractions import Fraction

from sage.categories.number_fields import NumberFields
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_arb import RealBallField
from sage.rings.real_double import RDF
from sage.rings.qqbar import AA, number_field_elements_from_algebraics

from ..features import pynormaliz_feature


RBF = RealBallField(64)
_NumberFields = NumberFields()

_backends = {}

_backends[int] = _backends[Fraction] = 'ppl'
_backends[ZZ] = _backends[QQ] = 'ppl'
def default_backend(base_ring):
    try:
        return _backends[base_ring]
    except KeyError:
        pass
    if base_ring in _NumberFields:
        return 'normaliz-NF' if pynormaliz_feature.is_present() else 'sage'
    raise ValueError('no suitable polytope backend for base ring {}'.format(base_ring))


class Cone:
    _name = 'none'

    def __init__(self, base_ring, cone):
        self._base_ring = base_ring
        self._cone = cone

    def __repr__(self):
        return 'Cone of dimension {} in ambient dimension {} made of {} facets (backend={})'.format(
                self.affine_dimension(), self.space_dimension(), len(self.ieqs()), self._name)

    def space_dimension(self):
        raise NotImplementedError

    def affine_dimension(self):
        raise NotImplementedError

    def ieqs(self):
        raise NotImplementedError

    def eqns(self):
        raise NotImplementedError

    def rays(self):
        raise NotImplementedError

    def lines(self):
        raise NotImplementedError

    def facets(self):
        raise NotImplementedError

    def add_constraint(self, constraint):
        from .linear_expression import ConstraintSystem
        cs = ConstraintSystem()
        cs.insert(constraint)
        return self.add_constraints(cs)

    def add_constraints(self, constraints):
        raise NotImplementedError


class Cone_ppl(Cone):
    r"""
    PPL cone wrapper.
    """
    _name = 'ppl'

    @staticmethod
    def _new(base_ring, ieqs, eqns):
        raise TypeError('do not use Cone_ppl._new')

    def __hash__(self):
        ieqs, eqns = self.ieqs_eqns()
        ieqs.sort()
        eqns.sort()
        return hash(tuple(tuple(ieq) for ieq in ieqs), tuple(tuple(eqn) for eqn in eqns))

    def copy(self):
        r"""
        EXAMPLES::

            sage:
        """
        import ppl
        return Cone_ppl(ZZ, ppl.C_Polyhedron(self._cone))

    def space_dimension(self):
        return self._cone.space_dimension()

    def affine_dimension(self):
        return self._cone.affine_dimension()

    def ieqs(self):
        return [[ZZ(x) for x in c.coefficients()] for c in self._cone.minimized_constraints() if c.is_inequality()]

    def eqns(self):
        return [[ZZ(x) for x in c.coefficients()] for c in self._cone.minimized_constraints() if c.is_equality()]

    def rays(self):
        return [[ZZ(x) for x in g.coefficients()] for g in self._cone.minimized_generators() if g.is_ray()]

    def lines(self):
        return [[ZZ(x) for x in g.coefficients()] for g in self._cone.minimized_generators() if g.is_line()]

    def add_constraint(self, constraint):
        r"""
        Add a single constraint to the cone.

        EXAMPLES::

            sage: from veerer.polyhedron import LinearExpressions, ConstraintSystem
            sage: L = LinearExpressions(ZZ)
            sage: x0, x1, x2, x3 = (L.variable(i) for i in range(4))
            sage: cs = ConstraintSystem()
            sage: cs.insert(x0 >= 0)
            sage: cs.insert(x1 >= 0)
            sage: cs.insert(x2 >= 0)
            sage: cs.insert(x3 >= 0)
            sage: cs.insert(x0 + x1 - 2*x2 + 3*x3 == 0)

            sage: cone = cs.cone('ppl')
            sage: cone.add_constraint(x0 == x1)
            Cone of dimension 2 in ambient dimension 4 made of 2 facets (backend=ppl)
        """
        import ppl
        cone = ppl.C_Polyhedron(self._cone)
        cone.add_constraint(constraint.ppl())
        return Cone_ppl(ZZ, cone)

    def add_constraints(self, constraints):
        r"""
        Add a system of constraints to the cone.
        """
        import ppl
        cone = ppl.C_Polyhedron(self._cone)
        for constraint in constraints:
            cone.add_constraint(constraint.ppl())
        return Cone_ppl(ZZ, cone)


class Cone_sage(Cone):
    r"""
    Sage cone wrapper.
    """
    _name = 'sage'

    @staticmethod
    def _new(ieqs, eqns):
        from sage.geometry.polyhedron.constructor import Polyhedron
        P = Polyhedron(ieqs=ieqs, eqns=eqns)
        return Cone_sage(P.base_ring(), P)

    def __hash__(self):
        return hash(self._cone)

    def space_dimension(self):
        return self._cone.ambient_dim()

    def affine_dimension(self):
        return self._cone.dim()

    def ieqs(self):
        return [ieq[1:] for ieq in self._cone.inequalities_list()]

    def eqns(self):
        return [eqn[1:] for eqn in self._cone.equations_list()]

    def rays(self):
        return self._cone.rays_list()

    def lines(self):
        return self._cone.lines_list()

    def add_constraints(self, cs):
        from sage.geometry.polyhedron.constructor import Polyhedron
        ieqs, eqns = cs.ieqs_eqns(self._cone.ambient_dim())
        new_cone = self._cone.intersection(Polyhedron(ieqs=ieqs, eqns=eqns))
        return Cone_sage(new_cone.base_ring(), new_cone)


class NFElementHandler:
    def __init__(self, nf):
        self._nf = nf

    def __call__(self, l):
        nf = self._nf
        l = list(l) + [0] * (nf.degree() - len(l))
        l = nf(l)
        return l


class Cone_normaliz(Cone):
    @property
    def _rational_handler(self):
        return lambda l: QQ((l[0], l[1]))

    @property
    def _nfelem_handler(self):
        return NFElementHandler(self._base_ring)

    def _nmz_result_raw(self, prop):
        from PyNormaliz import NmzResult
        return NmzResult(self._cone, prop)

    def _nmz_result(self, prop):
        from PyNormaliz import NmzResult
        return NmzResult(self._cone, prop,
                         RationalHandler=self._rational_handler,
                         NumberfieldElementHandler=self._nfelem_handler)

    def __hash__(self):
        ieqs, eqns = self.ieqs_eqns()
        ieqs.sort()
        eqns.sort()
        return hash(tuple(tuple(ieq) for ieq in ieqs), tuple(tuple(eqn) for eqn in eqns))

    def space_dimension(self):
        return self._nmz_result("EmbeddingDim")

    def affine_dimension(self):
        return self._nmz_result("Rank")

    def ieqs(self):
        return self._nmz_result("SupportHyperplanes")

    def eqns(self):
        return self._nmz_result("Equations")

    def ieqs_eqns(self):
        return self.ieqs(), self.eqns()

    def rays(self):
        return self._nmz_result("ExtremeRays")

    def lines(self):
        raise NotImplementedError


class Cone_normalizQQ(Cone_normaliz):
    r"""
    Normaliz cone wrapper.
    """
    _name = 'normaliz-QQ'

    @staticmethod
    def _new(ieqs, eqns):
        from PyNormaliz import NmzCone
        return Cone_normalizQQ(QQ, NmzCone(inequalities=ieqs, equations=eqns))

    def add_constraints(self, cs):
        ieqs = self._nmz_result_raw("SupportHyperplanes")
        eqns = self._nmz_result_raw("Equations")
        ieqs2, eqns2 = cs.ieqs_eqns(self.space_dimension(), homogeneous=True)
        ieqs2 = [[[int(x.numerator()), int(x.denominator())] for x in ieq] for ieq in ieqs2]
        eqns2 = [[[int(x.numerator()), int(x.denominator())] for x in eqn] for eqn in eqns2]
        return self._new(ieqs + ieqs2, eqns + eqns2)


class Cone_normalizNF(Cone_normaliz):
    _name='normaliz-NF'

    @staticmethod
    def _new(ieqs, eqns):
        from PyNormaliz import NmzCone
        base_ring = embedded_nf([x for l in (ieqs + eqns) for x in l])
        ieqs = [[str(base_ring(x)) for x in ieq] for ieq in ieqs]
        eqns = [[str(base_ring(x)) for x in eqn] for eqn in eqns]
        nf_data = nmz_number_field_data(base_ring)
        ans = Cone_normalizNF(base_ring, NmzCone(number_field=nf_data, inequalities=ieqs, equations=eqns))
        ans._nf_data = nf_data
        return ans

    def add_constraints(self, cs):
        from PyNormaliz import NmzCone
        ieqs = self._nmz_result("SupportHyperplanes")
        eqns = self._nmz_result("Equations")
        ieqs2, eqns2 = cs.ieqs_eqns(self.space_dimension(), homogeneous=True)
        base_ring = self._base_ring
        ieqs.extend(ieqs2)
        eqns.extend(eqns2)
        ieqs = [[str(base_ring(x)) for x in ieq] for ieq in ieqs]
        eqns = [[str(base_ring(x)) for x in eqn] for eqn in eqns]
        ans = Cone_normalizNF(self._base_ring, NmzCone(number_field=self._nf_data, inequalities=ieqs, equations=eqns))
        ans._nf_data = self._nf_data
        return ans


def embedded_nf(l):
    from sage.structure.element import get_coercion_model
    cm = get_coercion_model()
    K = cm.common_parent(*l)
    if K in _NumberFields:
        if not RDF.has_coerce_map_from(K):
            raise ValueError("invalid base ring {} (no embedding)".format(K))
        return K
    elif K == AA:
        K, ll, hom = number_field_elements_from_algebraics(l, embedded=True, minimal=True)
        if K == QQ:
            raise ValueError('rational base ring')
        return K


def nmz_number_field_data(base_ring):
    gen = base_ring.gen()
    s_gen = str(gen)
    emb = RBF(gen)
    return [str(base_ring.polynomial()).replace("x", s_gen), s_gen, str(emb)]
