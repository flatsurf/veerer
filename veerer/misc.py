r"""
Utility functions.
"""

from __future__ import print_function, absolute_import


def det2(u, v):
    return u[0]*v[1] - u[1]*v[0]

def flipper_edge(T, e):
    r"""
    EXAMPLES::

        sage: import flipper
        sage: from veerer.layout import flipper_edge

        sage: T = flipper.create_triangulation([(0r,1r,2r),(-1r,-2r,-3r)])
        sage: sorted([flipper_edge(T, e) for e in T.edges])
        [0, 1, 2, 3, 4, 5]
    """
    n = (3 * T.num_triangles)
    e = e.label
    return n * (e < 0) + e

def flipper_edge_perm(n):
    from array import array
    return array('l', range(n-1,-1,-1))

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
    p = array('l', [-1]*(2*n))
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

    r = K.lmbda.interval_approximation()
    l = r.lower * ZZ(10)**(-r.precision)
    u = r.upper * ZZ(10)**(-r.precision)

    p = QQx(K.polynomial.coefficients)
    s = AA.polynomial_root(p, RIF(l,u))
    return NumberField(p, name, embedding=s)

def flipper_nf_element_to_sage(x, K=None):
    r"""
    Convert a flipper nf element to Sage.
    
    EXAMPLES::

        sage: import flipper
        sage: from veerer.misc import flipper_nf_element_to_sage

        sage: S = flipper.load('S_2_1')
        sage: h = S.mapping_class('acBD')
        sage: F = h.flat_structure()
        sage: x = F.edge_vectors.values()[0].x
        sage: flipper_nf_element_to_sage(x)
        2*a^3 - 14*a^2 + 26*a - 14
    """
    if K is None:
        K = flipper_nf_to_sage(x.number_field)
    return K(x.linear_combination)

