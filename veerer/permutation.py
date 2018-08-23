r"""
Partial permutation on `\{0, 1, ..., n-1\}`.
"""
from __future__ import absolute_import, print_function

from array import array
from math import log
from six.moves import range

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

def permutation_to_perm(p):
    r"""
    Return a list on `[0, n-1]` from a permutation on `[1, n]`

    EXAMPLES::

        sage: from veerer.permutation import permutation_to_perm
        sage: permutation_to_perm(PermutationGroupElement([3,1,2]))
        [2, 0, 1]
    """
    return map(lambda x: x-1, p.domain())

def perm_to_permutation(l):
    r"""
    Returns a permutation on `[1, n]` from a list on `[0, n-1]`

    EXAMPLES::

        sage: from veerer.permutation import perm_to_permutation
        sage: perm_to_permutation([2,1,0])
        (1,3)
    """
    from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
    return PermutationGroupElement(map(lambda x: x+1, l))

def perm_id(n):
    r"""
    Return the identity permutation.

    EXAMPLES::

        sage: from veerer.permutation import perm_id

        sage: perm_id(4)
        array('l', [0, 1, 2, 3])
    """
    return array('l', range(n))

def perm_init(data):
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
            return array('l', [])
        elif isinstance(data[0], (tuple, list)):
            return perm_from_cycles(data)
        else:
            return array('l', (-1 if x is None else x for x in data))

    if isinstance(data, str):
        c = str_to_cycles(data)
        return perm_from_cycles(c)

    raise TypeError("The input must be list, tuple or string")

def perm_from_cycles(t):
    r"""
    Returns a permutation on `[0, n-1]` from a list of cycles on `[0, n-1]`

    EXAMPLES::

        sage: from veerer.permutation import perm_from_cycles

        sage: perm_from_cycles([[1,3,5],[0,2,4],[6]])
        array('l', [2, 3, 4, 5, 0, 1, 6])

        sage: perm_from_cycles([])
        array('l')
        sage: perm_from_cycles([[],[]])
        array('l')
    """
    if not any(tt for tt in t):
        return array('l', [])

    res = array('l', range(max(map(max, t))+1))

    for c in t:
        for j in range(len(c)-1):
            res[c[j]] = int(c[j+1])
        res[c[-1]] = int(c[0])

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

#####################################################################
# Serialization
#####################################################################

CHARS = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+-'
CHARS_INV = {j:i for i,j in enumerate(CHARS)}

def uint_to_base64_str(n, l=None):
    r"""
    EXAMPLES::

        sage: from veerer.permutation import uint_to_base64_str

        sage: uint_to_base64_str(15)
        'f'
        sage: uint_to_base64_str(15, 3)
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

        sage: from veerer.permutation import uint_from_base64_str, uint_to_base64_str

        sage: uint_from_base64_str('mqb')
        91787
        sage: uint_to_base64_str(91787)
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
        '4_3102'
        sage: s = perm_base64_str(range(2048))
        sage: s
        'w0_00010203...v+v-'
        sage: perm_from_base64_str(s) == array('l', range(2048))
        True
    """
    n = len(p)
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    return uint_to_base64_str(n) + '_' + ''.join(uint_to_base64_str(i, l) for i in p)

def perm_from_base64_str(s):
    r"""
    EXAMPLES::

        sage: from veerer.permutation import perm_from_base64_str, perm_base64_str
        sage: from array import array

        sage: p = array('l', [3,0,2,1])
        sage: s = perm_base64_str(p)
        sage: perm_from_base64_str(s) == p
        True

        sage: r = list(range(3000))
        sage: shuffle(r)
        sage: p = array('l', r)
        sage: perm_from_base64_str(perm_base64_str(p)) == p
        True
    """
    n, p = s.split('_')
    n = uint_from_base64_str(n)
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    if len(p) != n * l:
        raise ValueError('wrong length')
    return array('l', (uint_from_base64_str(p[i:i+l]) for i in range(0,len(p),l)))

#####################################################################
# Cycles and action
#####################################################################

def perm_cycles(p, singletons=True):
    r"""
    Return the cycle decomposition of ``p`` as a list of lists.

    INPUT:

    - ``p`` -- the permutation

    - ``singletons`` -- bool (default: ``True``) - return or not the singletons

    EXAMPLES::

        sage: from veerer.permutation import perm_cycles

        sage: perm_cycles([0,2,1])
        [[0], [1, 2]]
        sage: perm_cycles([0,2,1], False)
        [[1, 2]]

        sage: perm_cycles([2,-1,0])
        [[0, 2]]
    """
    seen = [False] * len(p)
    res = []

    for i in range(len(p)):
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

def perm_cycle_string(p, singletons=True):
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
                       perm_cycles(p, singletons)))

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


def perm_on_list(p, t):
    r"""
    Action of the permutation ``p`` on the list ``t``.

    EXAMPLES::

        sage: from veerer.permutation import perm_on_list
        sage: perm_on_list([2,1,3,0], [2,1,2,0])
        [3, 1, 3, 2]
    """
    return [p[i] for i in t]


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

def perm_invert(l):
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
    res = array('l', [0]*len(l))
    for i in range(len(l)):
        if l[i] == -1:
            res[i] = -1
        else:
            res[l[i]] = i
    return res

def perm_compose(p1, p2):
    r"""
    Returns the product ``p1 p2``.

    In the notation ``p1 p2`` we use the left action, in other words
    ``p1`` is applied first.

    EXAMPLES::

        sage: from veerer.permutation import perm_compose

        sage: perm_compose([0,2,1], [0,2,1])
        array('l', [0, 1, 2])
        sage: perm_compose([-1,2,3,1],[-1,2,1,3])
        array('l', [-1, 1, 3, 2])
    """
    r = array('l', [-1] * len(p1))
    for i in range(len(p1)):
        if p1[i] != -1 and p1[i] < len(p2):
            r[i] = p2[p1[i]]
    return r


def perm_compose_i(p1, p2):
    r"""
    Returns the product ``p1^{-1} p2^{-1}``.

    EXAMPLES::

        sage: from veerer.permutation import perm_compose_i

        sage: perm_compose_i([0,1,2],[1,2,0])
        array('l', [2, 0, 1])

    TESTS::

        sage: from veerer.permutation import perm_invert, perm_compose
        sage: from itertools import permutations

        sage: for p1 in permutations(range(4)):
        ....:     for p2 in permutations(range(4)):
        ....:         assert perm_compose_i(p1, p2) == perm_compose(perm_invert(p1), perm_invert(p2))
    """
    res = array('l', [-1]*len(p1))
    for i in range(len(p1)):
        res[p1[p2[i]]] = i

    return res

#####################################################################
# Transitivity test
#####################################################################

def perms_transitive_components(p):
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


def perms_are_transitive(p):
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
    n = len(p[0])

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

