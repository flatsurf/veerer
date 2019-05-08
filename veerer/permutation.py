r"""
Partial permutation on `\{0, 1, ..., n-1\}`.

TODO:

- We want much faster datastructure (ie C array)

- In many situations, we need a bitarray of the size of
  the permutation (conjugation, composition, etc). But
  such array should not be allocated each time the function
  is called.
"""
from __future__ import absolute_import, print_function
from six.moves import range, map, filter, zip

from array import array
from math import log

try:
    import sage.all
except ImportError:
    from random import shuffle, randint
else:
    from sage.misc.prandom import shuffle, randint

def argmin(l):
    r"""
    Return the position of the minimal element in the list ``l``.

    EXAMPLES::

        sage: from veerer.permutation import argmin
        sage: argmin([3,0,1,2])
        1
        sage: argmin([-1,3,5,-2,50])
        3
    """
    if not(l):
        raise ValueError('empty list')
    imin = 0
    jmin = l[0]
    for i,j in enumerate(l):
        if j < jmin:
            jmin = j
            imin = i
    return imin

#####################################################################
# Initialization and conversion
#####################################################################

def perm_check(l, n=None):
    r"""
    Checks that ``l`` is a partial permutation of `\{0, 1, ..., n-1\}`.

    INPUT:

    - ``n`` - integer (optional)

    EXAMPLES::

        sage: from veerer.permutation import perm_check
        sage: from array import array

    Good permutations::

        sage: perm_check(array('l', [1, 0, 3, 2]), 4)
        True
        sage: perm_check(array('l', [-1]), 1)
        True
        sage: perm_check(array('l', [-1, 3, -1, 1]), 4)
        True

        sage: perm_check(array('l', [1,0,-1,-1,-1]), 2)
        True

    Bad permutations::

        sage: perm_check(array('l', [1, 0, 3, 2]), 3)
        False
        sage: perm_check(array('l', [2, 0]))
        False
        sage: perm_check(array('l', [1, 0, 1]))
        False
    """
    if not isinstance(l, array):
        return False
    if n is None:
        n = len(l)
    else:
        n = int(n)

    seen = [False] * n
    for i in range(n):
        if l[i] == -1:
            continue
        if l[i] < 0 or l[i] >= n or seen[l[i]]:
            return False
        seen[l[i]] = True
    return True

def perm_id(n):
    r"""
    Return the identity permutation.

    EXAMPLES::

        sage: from veerer.permutation import perm_id

        sage: perm_id(4)
        array('l', [0, 1, 2, 3])
    """
    return array('l', range(n))

def perm_init(data, n=None, involution=None):
    """
    Returns a permutation from the given data.

    If data is a list of integers, then they are considered to be
    the images of the permutation. If ``data`` is a list of list
    then each list in ``data`` is considered as a cycle. Finally,
    string input in cycle notation is allowed.

    EXAMPLES::

        sage: from veerer.permutation import perm_init

    As a list of integers::

        sage: perm_init([3,2,1,4])
        array('l', [3, 2, 1, 4])
        sage: perm_init([3,1,None,0])
        array('l', [3, 1, -1, 0])

    Cycle notation (not mentioned elements are considered to be fixed
    point)::

        sage: perm_init(([2,1],[3,4,0]))
        array('l', [3, 2, 1, 4, 0])
        sage: perm_init([[2,1],[3,4,0]])
        array('l', [3, 2, 1, 4, 0])

    As a string::

        sage: perm_init('(0,1)(3,2)')
        array('l', [1, 0, 3, 2])

    Zerology::

        sage: perm_init([])
        array('l')
        sage: perm_init([[]])
        array('l')
        sage: perm_init('')
        array('l')
        sage: perm_init('()')
        array('l')
    """
    if isinstance(data, (array, tuple, list)):
        if not data:
            if n is not None:
                return array('l', range(n))
            else:
                return array('l', [])
        elif isinstance(data[0], (tuple, list)):
            return perm_from_cycles(data, n=n, involution=involution)
        else:
            return array('l', (-1 if x is None else x for x in data))

    if isinstance(data, str):
        c = str_to_cycles(data)
        return perm_from_cycles(c, n=n, involution=involution)

    if data.__module__.startswith('flipper'):
        if involution is None:
            raise ValueError("involution must be provided")
        from .misc import flipper_isometry_to_perm
        return flipper_isometry_to_perm(data, involution)

    raise TypeError("The input must be list, tuple or string")

def perm_from_cycles(t, n=None, involution=None):
    r"""
    Returns a permutation on `[0, n-1]` from a list of cycles on `[0, n-1]`

    INPUT:

    - ``t`` - cycles

    - ``n`` - optional domain size

    - ``involution`` (optional) - if provided use it to convert minus
      signs

    EXAMPLES::

        sage: from veerer.permutation import perm_from_cycles

        sage: perm_from_cycles([[1,3,5],[0,2,4],[6]])
        array('l', [2, 3, 4, 5, 0, 1, 6])

        sage: perm_from_cycles([])
        array('l')
        sage: perm_from_cycles([[],[]])
        array('l')

        sage: perm_from_cycles([[1,-2,0,3]], n=6, involution=[0,4,5,3,1,2])
        array('l', [3, 4, 2, 1, 0, 5])
    """
    if not any(tt for tt in t):
        return array('l', [])

    if n is None:
        n = max(map(max, t)) + 1
    else:
        n = int(n)

    res = array('l', range(n))

    for c in t:
        a = int(c[0])
        if a < 0:
            a = n+a if involution is None else involution[~a]
        for j in range(1,len(c)):
            b = int(c[j])
            if b < 0:
                b = n+b if involution is None else involution[~b]
            res[a] = b
            a = b
        b = int(c[0])
        if b < 0:
            b = n+b if involution is None else involution[~b]
        res[a] = b

    return res

def str_to_cycles(s):
    """
    Returns a list of cycles from a string.

    EXAMPLES::

        sage: from veerer.permutation import str_to_cycles
        sage: str_to_cycles('(0,1)')
        [[0, 1]]
        sage: str_to_cycles('(0,1)(3,2)')
        [[0, 1], [3, 2]]

        sage: str_to_cycles('()(0,1)()(2,3)')
        [[0, 1], [2, 3]]

        sage: str_to_cycles('(0,1,2)(~0,~1,~2)')
        [[0, 1, 2], [-1, -2, -3]]
    """
    r = []
    for c_str in s[1:-1].split(')('):
        if not c_str:
            continue
        r.append([~int(c[1:]) if c[0] == '~' else int(c) for c in c_str.replace(' ', '').split(',')])
    return r

def perm_random(n):
    r"""
    Return a random permutation.

    EXAMPLES::

        sage: from veerer.permutation import perm_random, perm_check
        sage: perm_check(perm_random(13), 13)
        True
    """
    r = list(range(n))
    shuffle(r)
    return array('l', r)

def perm_random_centralizer(p):
    r"""
    Return a random permutation in the centralizer of ``p``.

    EXAMPLES::

        sage: from veerer.permutation import *
        sage: p = perm_random(10)
        sage: q = perm_random_centralizer(p)
        sage: perm_compose(p, q) == perm_compose(q, p)
        True
    """
    if not p:
        return p

    cyc = perm_cycles(p)
    cyc.sort(key = lambda c: len(c))
    i = 0
    ans = array('l', [-1] * len(p))
    while i < len(cyc):
        j = i + 1
        s = len(cyc[i])
        while j < len(cyc) and len(cyc[j]) == s:
            j += 1

        # permutation of the cycles i, i+1, ..., j-1
        m = perm_random(j - i)

        for ii in range(i, j):
            jj = i + m[ii - i]
            shift = randint(0, s - 1)  # random shift
            for k in range(len(cyc[i])):
                ans[cyc[ii][k]] = cyc[jj][(k + shift) % s]

        # next round
        i = j

    return ans

def perm_random_conjugacy_class(c):
    r"""
    Return a random permutation with given conjugacy class ``c``.

    EXAMPLES::

        sage: from veerer.permutation import *

        sage: p = perm_random_conjugacy_class([5,2])
        sage: perm_cycle_type(p)
        [5, 2]

        sage: p = perm_random_conjugacy_class([7,3,3,1])
        sage: perm_cycle_type(p)
        [7, 3, 3, 1]
    """
    n = sum(c)
    r = list(range(n))
    shuffle(r)
    p = array('l', [-1]*n)
    i = 0
    for k in c:
        # add a k-cycle following the list r
        for j in range(i, i+k-1):
            p[r[j]] = r[j+1]
        p[r[i+k-1]] = r[i]
        i += k
    return p

#####################################################################
# Serialization
#####################################################################
# TODO: this is called often and would better be cythonized

CHARS = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+-'
CHARS_INV = {j:i for i,j in enumerate(CHARS)}

def uint_base64_str(n, l=None):
    r"""
    EXAMPLES::

        sage: from veerer.permutation import uint_base64_str

        sage: uint_base64_str(15)
        'f'
        sage: uint_base64_str(15, 3)
        '00f'
    """
    n = int(n)
    s = ''
    while n:
        s = CHARS[n % 64] + s
        n //= 64
    if l is not None:
        if len(s) > l:
            raise ValueError
        else:
            s = CHARS[0] * (l - len(s)) + s
    return s

def uint_from_base64_str(s):
    r"""
    EXAMPLES::

        sage: from veerer.permutation import uint_from_base64_str, uint_base64_str

        sage: uint_from_base64_str('mqb')
        91787
        sage: uint_base64_str(91787)
        'mqb'

        sage: uint_from_base64_str('00mqb')
        91787
    """
    n = 0
    d = 1
    for c in reversed(s):
        n += CHARS_INV[c] * d
        d *= 64
    return n

def perm_base64_str(p):
    r"""
    Make a canonical ASCII string out of ``p``.

    EXAMPLES::

        sage: from veerer.permutation import perm_base64_str, perm_from_base64_str
        sage: from array import array

        sage: perm_base64_str([3,1,0,2])
        '3102'
        sage: s = perm_base64_str(range(2048))
        sage: s
        '00010203...v+v-'
        sage: perm_from_base64_str(s, 2048) == array('l', range(2048))
        True
    """
    n = len(p)
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    return ''.join(uint_base64_str(i, l) for i in p)

def perm_from_base64_str(s, n):
    r"""
    EXAMPLES::

        sage: from veerer.permutation import perm_from_base64_str, perm_base64_str
        sage: from array import array

        sage: p = array('l', [3,0,2,1])
        sage: s = perm_base64_str(p)
        sage: perm_from_base64_str(s, 4) == p
        True

        sage: r = list(range(3000))
        sage: shuffle(r)
        sage: p = array('l', r)
        sage: perm_from_base64_str(perm_base64_str(p), 3000) == p
        True
    """
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    if len(s) != n * l:
        raise ValueError('wrong length')
    return array('l', (uint_from_base64_str(s[i:i+l]) for i in range(0,len(s),l)))

#####################################################################
# Cycles and action
#####################################################################

def perm_dense_cycles(p, n=None):
    r"""
    EXAMPLES::

        sage: from veerer.permutation import perm_dense_cycles

        sage: perm_dense_cycles([1,2,0])
        ([0, 0, 0], [3])

        sage: perm_dense_cycles([0,2,1])
        ([0, 1, 1], [1, 2])

        sage: perm_dense_cycles([2,1,0])
        ([0, 1, 0], [2, 1])
    """
    if n is None:
        n = len(p)
    deg = []
    res = [-1] * n
    k = 0
    for i in range(n):
        if res[i] != -1:
            continue
        d = 0
        while res[i] == -1:
            res[i] = k
            i = p[i]
            d += 1
        k += 1
        deg.append(d)
    return res, deg

def perm_cycles(p, singletons=True, n=None):
    r"""
    Return the cycle decomposition of ``p`` as a list of lists.

    INPUT:

    - ``p`` -- the permutation

    - ``singletons`` -- bool (default: ``True``) - return or not the singletons

    - ``n`` -- (optional) only use the first ``n`` elements of the permutation ``p``

    EXAMPLES::

        sage: from veerer.permutation import perm_cycles

        sage: perm_cycles([0,2,1])
        [[0], [1, 2]]
        sage: perm_cycles([0,2,1], False)
        [[1, 2]]

        sage: perm_cycles([2,-1,0])
        [[0, 2]]

        sage: perm_cycles([2,0,1,None,None], n=3)
        [[0, 2, 1]]
    """
    if n is None:
        n = len(p)
    elif n < 0 or n > len(p):
        raise ValueError

    seen = [False] * n
    res = []

    for i in range(n):
        if seen[i] or p[i] == -1:
            continue
        if p[i] == i and not singletons:
            continue
        cycle = []
        j = i
        while not seen[j]:
            seen[j] = True
            cycle.append(j)
            j = p[j]
        res.append(cycle)

    return res

def perm_num_cycles(p, n=None):
    r"""
    Return the number of cycles of the permutation ``p``.

    EXAMPLES::

        sage: from veerer.permutation import perm_num_cycles
        sage: perm_num_cycles([1,2,3,0])
        1
        sage: perm_num_cycles([0,2,3,1])
        2
        sage: perm_num_cycles([3,2,1,0])
        2
        sage: perm_num_cycles([3,1,2,0])
        3
        sage: perm_num_cycles([0,1,2,3])
        4
    """
    if n is None:
        n = len(p)
    seen = [False] * n
    ans = 0
    for i in range(n):
        if seen[i] or p[i] == -1:
            continue
        ans += 1
        j = i
        while not seen[j]:
            seen[j] = True
            j = p[j]
    return ans

def perm_cycle_type(p, n=None):
    r"""
    Return the lengths of the cycles of the permutation ``p`` in size of
    decreasing order.

    EXAMPLES::

        sage: from veerer.permutation import perm_cycle_type
        sage: perm_cycle_type([1,2,3,0])
        [4]
        sage: perm_cycle_type([0,2,3,1])
        [3, 1]
        sage: perm_cycle_type([3,2,1,0])
        [2, 2]
        sage: perm_cycle_type([3,1,2,0])
        [2, 1, 1]
        sage: perm_cycle_type([0,1,2,3])
        [1, 1, 1, 1]
    """
    if n is None:
        n = len(p)
    seen = [False] * n
    c = []
    for i in range(n):
        if seen[i] or p[i] == -1:
            continue
        k = 0
        j = i
        while not seen[j]:
            seen[j] = True
            k += 1
            j = p[j]
        c.append(k)
    c.sort(reverse=True)
    return c

def perm_cycle_string(p, singletons=True, n=None):
    r"""
    Returns a string representing the cycle decomposition of `p`

    EXAMPLES::

        sage: from veerer.permutation import perm_cycle_string
        sage: perm_cycle_string([0,2,1])
        '(0)(1,2)'
        sage: perm_cycle_string([0,2,1],False)
        '(1,2)'
    """
    return ''.join(map(lambda x: '('+','.join(map(str, x))+')',
                       perm_cycles(p, singletons, n)))

def perm_orbit(p, i):
    r"""
    Returns the orbit of an integer `i` under the permutation `p`

    EXAMPLES::

        sage: from veerer.permutation import perm_orbit
        sage: perm_orbit([0,3,1,2],2)
        [2, 1, 3]
    """
    res = [i]
    j = p[i]
    while j != i:
        res.append(j)
        j = p[j]
    return res

def perm_on_list(p, a, n=None, swap=None):
    r"""
    Inplace action of permutation on list-like objects.

    INPUT:

    - ``p`` - permutation

    - ``a`` - list, array

    - ``n`` - (optional) size of permutation

    - ``swap`` - (optional) a swap function

    EXAMPLES::

        sage: from veerer.permutation import *
        sage: l = [0,1,2,3,4]
        sage: p = [4,2,3,0,1]
        sage: perm_on_list(p,l)
        sage: l
        [3, 4, 1, 2, 0]

    Permutation action on matrix rows::

        sage: m1 = matrix(ZZ, 5, 5, 1)
        sage: m2 = matrix(ZZ, 5, 5, 1)
        sage: m = matrix(ZZ, 5, 5, 1)
        sage: p1 = perm_init([4,1,3,2,0])
        sage: p2 = perm_init([1,0,3,4,2])
        sage: perm_on_list(p1, m1, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: perm_on_list(p2, m2, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: perm_on_list(perm_compose(p1, p2), m, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: m == m2 * m1
        True
    """
    if n is None:
        n = len(p)
    seen = [False] * n
    for i in range(n):
        if seen[i]:
            continue
        seen[i] = True
        j = p[i]
        while seen[j] is False:
            if swap:
                swap(a, i, j)
            else:
                tmp = a[i]
                a[i] = a[j]
                a[j] = tmp
            seen[j] = True
            j = p[j]

# WARNING: this is NOT inplace
def perm_on_cyclic_list(p, t):
    r"""
    Action of the permutation ``p`` on the list ``t`` up to cyclic order.

    EXAMPLES::

        sage: from veerer.permutation import perm_on_cyclic_list

        sage: perm_on_cyclic_list([0,1,2], [2,1,2])
        [1, 2, 2]
        sage: perm_on_cyclic_list([0,1], [0,1,0,0,1,0,0,0,1,1])
        [0, 0, 0, 1, 1, 0, 1, 0, 0, 1]

        sage: a = [1, 0, 3, 2, 5, 4]
        sage: perm_on_cyclic_list(a, [0, 5, 3])
        [1, 4, 2]
        sage: perm_on_cyclic_list(a, [1, 4, 2])
        [0, 5, 3]

        sage: a1 = [0, 1, 4, 2, 3, 5]
        sage: a2 = [1, 5, 3, 4, 2, 0]
        sage: a3 = [2, 3, 1, 5, 0, 4]
        sage: a4 = [5, 0, 2, 3, 4, 1]
        sage: t1 = [0, 5, 1]
        sage: t2 = [2, 4, 3]
        sage: perm_on_cyclic_list(a1, t1) == perm_on_cyclic_list(a2, t1) == perm_on_cyclic_list(a4, t1) == t1
        True
        sage: perm_on_cyclic_list(a3, t1) == t2
        True
        sage: perm_on_cyclic_list(a3, t2) == t1
        True
    """
    res = r = [p[i] for i in t]
    # the thing below is very stupid!
    for i in range(1, len(r)):
        rr = r[i:] + r[:i]
        if rr < res:
            res = rr
    return res

#####################################################################
# Group operations
#####################################################################

def perm_invert(l, n=None):
    r"""
    Returns the inverse of the permutation ``l``.

    EXAMPLES::

        sage: from veerer.permutation import perm_invert

        sage: perm_invert([0, 3, 1, 2])
        array('l', [0, 2, 3, 1])

        sage: perm_invert([2, -1, 5, 0, -1, 3])
        array('l', [3, -1, 0, 5, -1, 2])

    TESTS::

        sage: from itertools import permutations
        sage: from veerer.permutation import perm_invert, perm_compose, perm_id

        sage: q = perm_id(3)
        sage: all(perm_compose(perm_invert(p),p) == q for p in permutations(range(3)))
        True
        sage: all(perm_compose(p,perm_invert(p)) == q for p in permutations(range(3)))
        True

    """
    if n is None:
        n = len(l)
    res = array('l', [0]*n)
    for i in range(n):
        if l[i] == -1:
            res[i] = -1
        else:
            res[l[i]] = i
    return res

def perm_compose(p1, p2, n=None):
    r"""
    Returns the product ``p1 p2``.

    In the notation ``p1 p2`` we use the right action, in other words
    ``p1`` is applied first.

    EXAMPLES::

        sage: from veerer.permutation import perm_compose

        sage: perm_compose([0,2,1], [0,2,1])
        array('l', [0, 1, 2])
        sage: perm_compose([-1,2,3,1],[-1,2,1,3])
        array('l', [-1, 1, 3, 2])

        sage: perm_compose([1,0,2,None,None], [2,1,0,None], 3)
        array('l', [1, 2, 0])
    """
    if n is None:
        n = len(p1)
    r = array('l', [-1] * n)
    for i in range(n):
        if p1[i] != -1:
            r[i] = p2[p1[i]]
    return r

perm_compose_00 = perm_compose

def perm_compose_10(p1, p2, n=None):
    r"""
    Return the product ``p1^(-1) p2``

    EXAMPLES::

        sage: from veerer.permutation import perm_compose, perm_invert, perm_compose_10
        sage: p1 = [0,5,2,1,3,4]
        sage: p2 = [3,1,5,4,2,0]
        sage: perm_compose_10(p1, p2) == perm_compose(perm_invert(p1), p2)
        True
        sage: shuffle(p1)
        sage: shuffle(p2)
        sage: perm_compose_10(p1, p2) == perm_compose(perm_invert(p1), p2)
        True
    """
    if n is None:
        n = len(p1)
    r = array('l', [-1] * n)
    for i in range(n):
        r[p1[i]] = p2[i]
    return r

def perm_compose_01(p1, p2, n=None):
    r"""
    Return the product ``p1 p2^(-1)``

    EXAMPLES::

        sage: from veerer.permutation import perm_compose, perm_invert, perm_compose_01
        sage: p1 = [0,5,2,1,3,4]
        sage: p2 = [3,1,5,4,2,0]
        sage: perm_compose_01(p1, p2) == perm_compose(p1, perm_invert(p2)) # not tested
        True
        sage: shuffle(p1)
        sage: shuffle(p2)
        sage: perm_compose_01(p1, p2) == perm_compose(p1, perm_invert(p2)) # not tested
        True
    """
    raise NotImplementedError

def perm_compose_11(p1, p2, n=None):
    r"""
    Return the product ``p1^(-1) p2^(-1)

    EXAMPLES::

        sage: from veerer.permutation import perm_compose, perm_invert, perm_compose_11
        sage: p1 = [0,5,2,1,3,4]
        sage: p2 = [3,1,5,4,2,0]
        sage: perm_compose_11(p1, p2) == perm_compose(perm_invert(p1), perm_invert(p2))
        True
        sage: shuffle(p1)
        sage: shuffle(p2)
        sage: perm_compose_11(p1, p2) == perm_compose(perm_invert(p1), perm_invert(p2))
        True

    TESTS::

        sage: from veerer.permutation import perm_invert, perm_compose
        sage: from itertools import permutations

        sage: for p1 in permutations(range(4)):
        ....:     for p2 in permutations(range(4)):
        ....:         assert perm_compose_11(p1, p2) == perm_compose(perm_invert(p1), perm_invert(p2))
    """
    if n is None:
        n = len(p1)
    r = array('l', [-1] * n)
    for i in range(n):
        r[p1[p2[i]]] = i
    return r

def perm_conjugate(p1, p2, n=None):
    r"""
    Conjugate ``p1`` by ``p2``.

    Let call ``res`` the output of this function. If ``p1`` was
    mapping ``a`` to ``b`` then ``res`` will map ``p2[a]``
    to ``p2[b]``.

    EXAMPLES::

        sage: from veerer.permutation import perm_conjugate, perm_random

        sage: p1 = perm_random(23)
        sage: p2 = perm_random(23)
        sage: res = perm_conjugate(p1, p2)
        sage: res[p2[14]] == p2[p1[14]]
        True
        sage: res[p2[19]] == p2[p1[19]]
        True
    """
    if n is None:
        n = len(p1)
    res = array('l', [-1] * n)
    for i in range(n):
        res[p2[i]] = p2[p1[i]]
    return res

#####################################################################
# Transitivity test
#####################################################################

def perms_transitive_components(p, n=None):
    r"""
    Return the list of transitive components of the subgroup generated by the
    permutations ``p``.

    INPUT:

    - ``p`` -- a list of permutations given as lists

    EXAMPLES::

        sage: from veerer.permutation import perms_transitive_components

        sage: perms_transitive_components([[1,0,2,3],[0,1,3,2]])
        [(0, 1), (2, 3)]

        sage: perms_transitive_components([[2,3,0,1]])
        [(0, 2), (1, 3)]

        sage: perms_transitive_components([[3,1,2,0], [0,3,2,1], [0,1,3,2]])
        [(0, 1, 2, 3)]
    """
    if n is None:
        n = len(p[0])
    seen = [-1] * n
    cc_num = 0
    for i in range(n):
        if seen[i] != -1:
            continue

        todo = [i]
        seen[i] = cc_num
        while todo:
            j = todo.pop()
            for pp in p:
                k = pp[j]
                if seen[k] == -1:
                    todo.append(k)
                    seen[k] = cc_num
        cc_num += 1

    return [tuple(i for i in range(n) if seen[i] == j) for j in range(cc_num)]


def perms_are_transitive(p, n=None):
    """
    Test whether the group generated by the permutations in ``p`` is transitive.

    We assume that the list of partial permutations act on
    the same domain (ie the -1 occur at the same positions).

    INPUT:

    - ``p`` - a list of permutations of `[0, n-1]`

    EXAMPLES::

        sage: from veerer.permutation import perms_are_transitive
        sage: perms_are_transitive([[0,1,2],[0,2,1]])
        False
        sage: perms_are_transitive([[0,1,2],[1,2,0]])
        True

        sage: p0 = [0,2,3,1,7,5,6,4]
        sage: p1 = [7,1,2,3,4,5,6,0]
        sage: p2 = [6,1,2,3,4,5,0,7]
        sage: p3 = [1,0,2,3,4,5,6,7]
        sage: p4 = [0,1,4,5,2,3,6,7]
        sage: perms_are_transitive([p0,p1,p2,p3,p4])
        True
        sage: perms_are_transitive([p0,p1,p2,p3])
        False
    """
    if not p:
        raise ValueError("empty list")

    p0 = p[0]
    if n is None:
        n = len(p0)

    # compute the connected component of 0
    cc0 = [True if j == -1 else False for j in p0]
    todo = [0]
    cc0[0] = True
    while todo:
        j = todo.pop()
        for pp in p:
            k = pp[j]
            if cc0[k] is False:
                todo.append(k)
                cc0[k] = True

    return all(cc0)


def perms_relabel(p, m):
    """
    Relabel the list of permutations ``p`` according to ``m``.

    INPUT:

    - ``p`` - a list of permutations

    - ``m`` - the relabeling permutation

    EXAMPLES::

        sage: from surface_dynamics.misc.constellation import perms_relabel
        sage: perms_relabel([[0,1,2],[0,2,1]],[2,1,0])
        [[0, 1, 2], [1, 0, 2]]
    """
    q = [k[:] for k in p]
    for i in range(len(m)):
        for j in range(len(p)):
            q[j][m[i]] = m[p[j][i]]
    return q


def perms_canonical_labels_from(x, y, j0):
    r"""
    Return canonical labels for ``x``, ``y`` that starts at ``j0``.

    .. WARNING:

        The group generated by ``x`` and the elements of ``y`` should be
        transitive.

    INPUT:

    - ``x`` -- list - a permutation of `[0, ..., n]` as a list

    - ``y`` -- list of permutations of `[0, ..., n]` as a list of lists

    - ``j0`` -- an index in [0, ..., n]

    OUTPUT:

    mapping: a permutation that specify the new labels

    EXAMPLES::

        sage: from surface_dynamics.misc.constellation import perms_canonical_labels_from
        sage: perms_canonical_labels_from([0,1,2],[[1,2,0]],0)
        [0, 1, 2]
        sage: perms_canonical_labels_from([1,0,2],[[2,0,1]],0)
        [0, 1, 2]
        sage: perms_canonical_labels_from([1,0,2],[[2,0,1]],1)
        [1, 0, 2]
        sage: perms_canonical_labels_from([1,0,2],[[2,0,1]],2)
        [2, 1, 0]
    """
    n = len(x)

    k = 0
    mapping = [None] * n
    waiting = [[] for i in range(len(y))]

    while k < n:
        # initialize at j0
        mapping[j0] = k
        waiting[0].append(j0)
        k += 1
        # complete x cycle from j0
        j = x[j0]
        while j != j0:
            mapping[j] = k
            waiting[0].append(j)
            k += 1
            j = x[j]

        # find another guy
        l = 0
        while l < len(waiting):
            i = 0
            while i < len(waiting[l]):
                j1 = waiting[l][i]
                if mapping[y[l][j1]] is None:
                    break
                i += 1

            if i == len(waiting[l]):  # not found: go further in waiting
                if l < len(waiting)-1:
                    waiting[l+1].extend(waiting[l])
                waiting[l] = []
                l += 1
                i = 0

            else:  # found: complete cycle from new guy
                j0 = y[l][j1]
                if l < len(waiting)-1:
                    waiting[l+1].extend(waiting[l][:i+1])
                del waiting[l][:i+1]
                break

    return mapping


def perms_canonical_labels(p, e=None):
    assert(len(p) > 1)
    n = len(p[0])

    c_win = None
    m_win = range(n)

    x = p[0]
    y = p[1:]

    if e is None:
        e = range(n)

    # get canonical label from i in to_test and compare
    while e:
        i = e.pop()
        m_test = perms_canonical_labels_from(x, y, i)
        c_test = perms_relabel(p, m_test)
        if c_win is None or c_test < c_win:
            c_win = c_test
            m_win = m_test

    return c_win, m_win
