r"""
Utility functions.
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2018-2023 Vincent Delecroix
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


def det2(u, v):
    return u[0]*v[1] - u[1]*v[0]


def flipper_edge(T, e):
    r"""
    EXAMPLES::

        sage: from veerer.layout import flipper_edge  # random output due to deprecation warnings from realalg
        sage: import flipper                                                # optional - flipper
        sage: T = flipper.create_triangulation([(0r,1r,2r),(-1r,-2r,-3r)])  # optional - flipper
        sage: sorted([flipper_edge(T, e) for e in T.edges])                 # optional - flipper
        [0, 1, 2, 3, 4, 5]
    """
    n = (3 * T.num_triangles)
    e = e.label
    return n * (e < 0) + e


def flipper_edge_perm(n):
    from array import array
    return array('i', range(n-1,-1,-1))


def flipper_face_edge_perms(T):
    r"""
    Return a pair ``(face permutation, edge permutation)`` from a flipper triangulation.
    """
    from .permutation import perm_init
    n = 3 * T.num_triangles # number of half edges

    # extract triangulation
    triangles = [(flipper_edge(T, e), flipper_edge(T, f), flipper_edge(T, g)) for e,f,g in T]
    return perm_init(triangles), flipper_edge_perm(n)


def flipper_isometry_to_perm(isom, ep, inv=False):
    r"""
    Question: how do we create an isometry in flipper?
    """
    from array import array
    n = isom.zeta
    p = array('i', [-1]*(2*n))
    if inv:
        dic = isom.inverse_label_map
    else:
        dic = isom.label_map

    for i,j in dic.items():
        if i < 0:
            i = ep[~i]
        if j < 0:
            j = ep[~j]
        p[i] = j
    return p


def flipper_nf_to_sage(K, name='a'):
    r"""
    Convert a flipper number field to Sage.

    .. NOTE::

        Currently, the code is not careful at all with root isolation.
    """
    from sage.rings.all import ZZ, QQ, AA, RIF
    from sage.rings.number_field.number_field import NumberField
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    QQx = PolynomialRing(QQ, 'x')

    r = K.lmbda.interval()
    l = r.lower * ZZ(10)**(-r.precision)
    u = r.upper * ZZ(10)**(-r.precision)

    p = QQx(K.coefficients)
    s = AA.polynomial_root(p, RIF(l,u))
    return NumberField(p, name, embedding=s)


def flipper_nf_element_to_sage(x, K=None):
    r"""
    Convert a flipper nf element to Sage.

    EXAMPLES::

        sage: from veerer.misc import flipper_nf_element_to_sage
        sage: import flipper                            # optional - flipper
        sage: S = flipper.load('S_2_1')                 # optional - flipper
        sage: h = S.mapping_class('acBD')               # optional - flipper
        sage: F = h.flat_structure()                    # optional - flipper
        sage: for v in F.edge_vectors.values():         # optional - flipper
        ....:     x = flipper_nf_element_to_sage(v.x)   # optional - flipper
        ....:     y = flipper_nf_element_to_sage(v.y)   # optional - flipper
    """
    from sage.rings.rational_field import QQ
    if K is None:
        K = flipper_nf_to_sage(x.field)
    coeffs = list(map(QQ, x.coefficients))
    if K.degree() != len(coeffs):
        coeffs.extend([0] * (K.degree() - len(coeffs)))
    return K(list(map(QQ, coeffs)))


def rays_to_ppl_cone(rays):
    gs = ppl.Generator_System()
    gs.insert(ppl.point())
    for r in rays:
        gs.insert(ppl.ray(sum(coeff * ppl.Variable(i) for i,coeff in enumerate(r) if coeff)))
    return ppl.C_Polyhedron(gs)
