r"""
Utility functions.
"""


def det2(u, v):
    return u[0]*v[1] - u[1]*v[0]

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

