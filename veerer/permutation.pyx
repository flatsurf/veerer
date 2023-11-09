r"""
Partial permutations on `\{0, 1, ..., n-1\}`.

TODO:

- In many situations, we need a bitarray of the size of
  the permutation (conjugation, composition, etc). But
  such array would better not be allocated each time the
  function is called.
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2018 Mark Bell
#                     2018-2023 Vincent Delecroix
#                     2018 Saul Schleimer
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

from cpython cimport array
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

        sage: from veerer.permutation import argmin  # random output due to deprecation warnings in realalg
        sage: argmin([3, 0, 1, 2])
        1
        sage: argmin([-1, 3, 5, -2, 50])
        3
    """
    if not l:
        raise ValueError('empty list')
    imin = 0
    jmin = l[0]
    for i, j in enumerate(l):
        if j < jmin:
            jmin = j
            imin = i
    return imin


def least_rotation(S):
    """
    Return the pair ``(index of smallest rotation, period)`` of the list ``S``.

    Implementation of Booth's algorithm.

    EXAMPLES::

        sage: from veerer.permutation import least_rotation

        sage: least_rotation([1,0])
        (1, 2)
        sage: least_rotation([0,1,0])
        (2, 3)
        sage: least_rotation([0,1,1,0,1])
        (3, 5)
        sage: least_rotation([0,1,0,1,1])
        (0, 5)

    Adding some period to the above examples::

        sage: least_rotation([1,0]*4)
        (1, 2)
        sage: least_rotation([0,1,0]*5)
        (2, 3)
        sage: least_rotation([0,1,1,0,1]*3)
        (3, 5)
        sage: least_rotation([0,1,0,1,1]*2)
        (0, 5)
    """
    l = len(S)
    S = S + S
    f = [-1] * len(S)  # failure function
    k = 0              # least rotation of string found so far
    period = 0
    for j in range(1, len(S)):
        sj = S[j]
        i = f[j - k - 1]
        while i != -1 and sj != S[k + i + 1]:
            if sj < S[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if sj != S[k + i + 1]:  # if sj != S[k+i+1], then i == -1
            if sj < S[k]:  # better index
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    # NOTE: the loop below could probably be included in the above loop
    for period in range(l):
        if f[period+l] == l:
            return (k, period)
    return (k, l)

#####################################################################
# Initialization and conversion
#####################################################################

def perm_check(l, int n=-1, involution=None):
    r"""
    Checks that ``l`` is a partial permutation of `\{0, 1, ..., n-1\}`.

    INPUT:

    - ``n`` - integer (optional)

    EXAMPLES::

        sage: from veerer.permutation import perm_check
        sage: from array import array

    Good permutations::

        sage: perm_check(array('i', [1, 0, 3, 2]), 4)
        True
        sage: perm_check(array('i', [-1]), 1)
        True
        sage: perm_check(array('i', [-1, 3, -1, 1]), 4)
        True

        sage: perm_check(array('i', [1,0,-1,-1,-1]), 2)
        True

    Bad permutations::

        sage: perm_check(array('i', [1, 0, 3, 2]), 3)
        False
        sage: perm_check(array('i', [2, 0]))
        False
        sage: perm_check(array('i', [1, 0, 1]))
        False

    With involution::

        sage: perm_check(array('i', [2,1,0]), involution=array('i', [2,1,0]))
        True
        sage: perm_check(array('i', [1,0,2]), involution=array('i', [2,1,0]))
        False
    """
    if not isinstance(l, array.array):
        return False
    if n == -1:
        n = len(l)
    else:
        if len(l) < n:
            return False

    seen = [False] * n
    for i in range(n):
        if l[i] == -1:
            continue
        if l[i] < 0 or l[i] >= n or seen[l[i]]:
            return False
        seen[l[i]] = True

    if involution is not None:
        for i in range(n):
            if l[involution[i]] != involution[l[i]]:
                return False

    return True


def perm_id(int n):
    r"""
    Return the identity permutation.

    EXAMPLES::

        sage: from veerer.permutation import perm_id

        sage: perm_id(4)
        array('i', [0, 1, 2, 3])
    """
    return array.array('i', range(n))


def perm_init(data, int n=-1, involution=None):
    """
    Return a permutation from the given data.

    If data is a list of integers, then they are considered to be
    the images of the permutation. If ``data`` is a list of list
    then each list in ``data`` is considered as a cycle. Finally,
    string input in cycle notation is allowed.

    EXAMPLES::

        sage: from veerer.permutation import perm_init

    As a list of integers::

        sage: perm_init([3,2,1,4])
        array('i', [3, 2, 1, 4])
        sage: perm_init([3,1,None,0])
        array('i', [3, 1, -1, 0])

    Cycle notation (not mentioned elements are considered to be fixed
    point)::

        sage: perm_init(([2,1],[3,4,0]))
        array('i', [3, 2, 1, 4, 0])
        sage: perm_init([[2,1],[3,4,0]])
        array('i', [3, 2, 1, 4, 0])

    As a string::

        sage: perm_init('(0,1)(3,2)')
        array('i', [1, 0, 3, 2])

    Initialize a permutation in the centralizer of an involution::

        sage: perm_init('(0,2)(1,3)', involution=[0,4,2,5,1,3])
        array('i', [2, 3, 0, 1, 5, 4])
        sage: perm_init('(0,~1)', involution=[0,1])
        array('i', [1, 0])
        sage: perm_init('(0,~1)', involution=[2,3,0,1])
        array('i', [3, 2, 1, 0])
        sage: perm_init('(0,3)', involution=[0,4,2,5,1,3])
        Traceback (most recent call last):
        ...
        ValueError: invalid input

    Zerology::

        sage: perm_init([])
        array('i')
        sage: perm_init([[]])
        array('i')
        sage: perm_init('')
        array('i')
        sage: perm_init('()')
        array('i')
    """
    if n == -1 and involution is not None:
        n = len(involution)
    if isinstance(data, (array.array, tuple, list)):
        if not data:
            if n is not None:
                return array.array('i', range(n))
            else:
                return array.array('i', [])
        elif isinstance(data[0], (tuple, list)):
            return perm_from_cycles(data, n=n, involution=involution)
        else:
            return array.array('i', (-1 if x is None else x for x in data))

    if isinstance(data, str):
        c = str_to_cycles(data)
        return perm_from_cycles(c, n=n, involution=involution)

    # TODO: test flipper conversion
    if data.__module__.startswith('flipper'):
        if involution is None:
            raise ValueError("involution must be provided")
        from .misc import flipper_isometry_to_perm
        return flipper_isometry_to_perm(data, involution)

    raise TypeError("The input must be list, tuple or string")


def perm_from_cycles(t, int n=-1, involution=None):
    r"""
    Return a permutation on `[0, n-1]` from a list of cycles on `[0, n-1]`

    INPUT:

    - ``t`` - cycles

    - ``n`` - optional domain size

    - ``involution`` (optional) - if provided use it to convert minus
      signs

    EXAMPLES::

        sage: from veerer.permutation import perm_from_cycles

        sage: perm_from_cycles([[1,3,5],[0,2,4],[6]])
        array('i', [2, 3, 4, 5, 0, 1, 6])

        sage: perm_from_cycles([])
        array('i')
        sage: perm_from_cycles([[],[]])
        array('i')

        sage: perm_from_cycles([[1,-2],[0,3]], n=6, involution=[0,4,5,3,1,2])
        array('i', [3, 4, 2, 0, 1, 5])
    """
    if not any(tt for tt in t):
        return array.array('i', [])

    if n == -1:
        n = max(map(max, t)) + 1

    res = array.array('i', range(n))

    for c in t:
        a = int(c[0])
        if a < 0:
            a = n+a if involution is None else involution[~a]
        for j in range(1,len(c)):
            b = int(c[j])
            if b < 0:
                b = n+b if involution is None else involution[~b]
            res[a] = b
            if involution is not None:
                if (a == involution[a]) != (b == involution[b]):
                    raise ValueError("invalid input")
                res[involution[a]] = involution[b]
            a = b
        b = int(c[0])
        if b < 0:
            b = n+b if involution is None else involution[~b]

        res[a] = b
        if involution is not None:
            if (a == involution[a]) != (b == involution[b]):
                raise ValueError("invalid input")
            res[involution[a]] = involution[b]

    return res


def str_to_cycles(s):
    """
    Return a list of cycles from a string.

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


def perm_random(int n):
    r"""
    Return a random permutation.

    EXAMPLES::

        sage: from veerer.permutation import perm_random, perm_check
        sage: perm_check(perm_random(13), 13)
        True
    """
    r = list(range(n))
    shuffle(r)
    return array.array('i', r)


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
    cdef array.array ans = array.array('i', [-1] * len(p))
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

        sage: from veerer.permutation import perm_random_conjugacy_class, perm_cycle_type

        sage: p = perm_random_conjugacy_class([5, 2])
        sage: perm_cycle_type(p)
        [5, 2]

        sage: p = perm_random_conjugacy_class([7, 3, 3, 1])
        sage: perm_cycle_type(p)
        [7, 3, 3, 1]
    """
    n = sum(c)
    r = list(range(n))
    shuffle(r)
    p = array.array('i', [-1]*n)
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
        sage: perm_from_base64_str(s, 2048) == array('i', range(2048))
        True

        sage: perm_from_base64_str('vdh0keigmcjfpxtnrwsouyqba987654321zl', 36)
        array('i', [31, 13, 17, 0, 20, 14, 18, 16, 22, 12, 19, 15, 25, 33, 29, 23, 27, 32, 28, 24, 30, 34, 26, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 35, 21])
    """
    n = len(p)
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    return ''.join(uint_base64_str(i, l) for i in p)


def perm_from_base64_str(s, n):
    r"""
    EXAMPLES::

        sage: from veerer.permutation import perm_from_base64_str, perm_base64_str
        sage: from array import array

        sage: p = array('i', [3,0,2,1])
        sage: s = perm_base64_str(p)
        sage: perm_from_base64_str(s, 4) == p
        True

        sage: r = list(range(3000))
        sage: shuffle(r)
        sage: p = array('i', r)
        sage: perm_from_base64_str(perm_base64_str(p), 3000) == p
        True
    """
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    if len(s) != n * l:
        raise ValueError('wrong length')
    return array.array('i', (uint_from_base64_str(s[i:i+l]) for i in range(0,len(s),l)))

#####################################################################
# Boolean properties
#####################################################################

def perm_is_one(array.array p, int n=-1):
    r"""
    Return whether ``p`` is the identity permutation.
    """
    if n == -1:
        n = len(p)
    cdef int i
    for i in range(n):
        if p.data.as_ints[i] != i:
            return False
    return True

#####################################################################
# Cycles and action
#####################################################################

def perm_dense_cycles(array.array p, int n=-1):
    r"""
    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_dense_cycles

        sage: perm_dense_cycles(array('i', [1,2,0]))
        ([0, 0, 0], [3])

        sage: perm_dense_cycles(array('i', [0,2,1]))
        ([0, 1, 1], [1, 2])

        sage: perm_dense_cycles(array('i', [2,1,0]))
        ([0, 1, 0], [2, 1])
    """
    if n == -1:
        n = len(p)
    cdef list deg = []
    cdef list res = [-1] * n
    k = 0
    for i in range(n):
        if res[i] != -1:
            continue
        d = 0
        while res[i] == -1:
            res[i] = k
            i = p.data.as_ints[i]
            d += 1
        k += 1
        deg.append(d)
    return res, deg


def perm_cycles(array.array p, singletons=True, int n=-1):
    r"""
    Return the cycle decomposition of ``p`` as a list of lists.

    INPUT:

    - ``p`` -- the permutation

    - ``singletons`` -- bool (default: ``True``) - return or not the singletons

    - ``n`` -- (optional) only use the first ``n`` elements of the permutation ``p``

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_cycles

        sage: perm_cycles(array('i', [0,2,1]))
        [[0], [1, 2]]
        sage: perm_cycles(array('i', [0,2,1]), False)
        [[1, 2]]

        sage: perm_cycles(array('i', [2,-1,0]))
        [[0, 2]]

        sage: perm_cycles(array('i', [2,0,1,-1,-1]), n=3)
        [[0, 2, 1]]
    """
    if n == -1:
        n = len(p)
    elif n < 0 or n > len(p):
        raise ValueError

    cdef array.array seen = array.clone(p, n, True)
    cdef list res = []
    cdef list cycle
    cdef int i, j

    for i in range(n):
        if seen.data.as_ints[i] or p.data.as_ints[i] == -1:
            continue
        if p.data.as_ints[i] == i and not singletons:
            continue
        cycle = []
        j = i
        while not seen.data.as_ints[j]:
            seen[j] = True
            cycle.append(j)
            j = p.data.as_ints[j]
        res.append(cycle)

    return res


def perm_cycles_lengths(array.array p, int n=-1):
    r"""
    Return the array of orbit sizes.

    INPUT:

    - ``p`` -- the permutation

    - ``n`` -- (optional) only use the first ``n`` elements of the permutation ``p``

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_cycles_lengths

        sage: perm_cycles_lengths(array('i', [0,2,1]))
        array('i', [1, 2, 2])
        sage: perm_cycles_lengths(array('i', [2,-1,0]))
        array('i', [2, -1, 2])
        sage: perm_cycles_lengths(array('i', [2,0,1,-1,-1]), n=3)
        array('i', [3, 3, 3])
    """
    if n == -1:
        n = len(p)
    elif n < 0 or n > len(p):
        raise ValueError

    cdef array.array lengths = array.array('i', [-1] * n)
    cdef int i, j, m

    for i in range(n):
        if lengths.data.as_ints[i] != -1 or p.data.as_ints[i] == -1:
            continue
        j = i
        m = 0
        while lengths.data.as_ints[j] == -1:
            lengths.data.as_ints[j] = 0
            j = p.data.as_ints[j]
            m += 1
        while lengths.data.as_ints[j] == 0:
            lengths.data.as_ints[j] = m
            j = p[j]

    return lengths


def perm_num_cycles(array.array p, n=-1):
    r"""
    Return the number of cycles of the permutation ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_num_cycles

        sage: perm_num_cycles(array('i', [1,2,3,0]))
        1
        sage: perm_num_cycles(array('i', [0,2,3,1]))
        2
        sage: perm_num_cycles(array('i', [3,2,1,0]))
        2
        sage: perm_num_cycles(array('i', [3,1,2,0]))
        3
        sage: perm_num_cycles(array('i', [0,1,2,3]))
        4
    """
    if n == -1:
        n = len(p)
    cdef array.array seen = array.clone(p, n, True)
    ans = 0
    cdef int i, j
    for i in range(n):
        if seen.data.as_ints[i] or p.data.as_ints[i] == -1:
            continue
        ans += 1
        j = i
        while not seen.data.as_ints[j]:
            seen.data.as_ints[j] = 1
            j = p.data.as_ints[j]
    return ans


def perm_cycle_type(array.array p, int n=-1):
    r"""
    Return the lengths of the cycles of the permutation ``p`` in size of
    decreasing order.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_cycle_type
        sage: perm_cycle_type(array('i', [1,2,3,0]))
        [4]
        sage: perm_cycle_type(array('i', [0,2,3,1]))
        [3, 1]
        sage: perm_cycle_type(array('i', [3,2,1,0]))
        [2, 2]
        sage: perm_cycle_type(array('i', [3,1,2,0]))
        [2, 1, 1]
        sage: perm_cycle_type(array('i', [0,1,2,3]))
        [1, 1, 1, 1]
    """
    if n == -1:
        n = len(p)
    cdef array.array seen = array.clone(p, n, True)
    cdef list c = []
    cdef int i, j, k
    for i in range(n):
        if seen.data.as_ints[i] or p.data.as_ints[i] == -1:
            continue
        k = 0
        j = i
        while not seen.data.as_ints[j]:
            seen.data.as_ints[j] = 1
            k += 1
            j = p.data.as_ints[j]
        c.append(k)
    c.sort(reverse=True)
    return c


def perm_cycle_string(array.array p, singletons=True, n=-1, involution=None):
    r"""
    Return a string representing the cycle decomposition of `p`

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_cycle_string

        sage: perm_cycle_string(array('i', [0,2,1]))
        '(0)(1,2)'
        sage: perm_cycle_string(array('i', [0,2,1]), False)
        '(1,2)'
    """
    if involution:
        elt = lambda e: '%d'%e if e <= involution[e] else '~%d'%involution[e]
    else:
        elt = str

    return ''.join(map(lambda x: '('+','.join(map(elt, x))+')',
                       perm_cycles(p, singletons, n)))


def perm_orbit(array.array p, int i):
    r"""
    Return the orbit of ``i`` under the permutation ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_orbit

        sage: perm_orbit(array('i', [0,3,1,2]), 2)
        [2, 1, 3]
    """
    cdef int j
    cdef list res = [i]
    j = p.data.as_ints[i]
    while j != i:
        res.append(j)
        j = p.data.as_ints[j]
    return res


def perm_orbit_size(array.array p, int i):
    r"""
    Return the size of the orbit of ``i`` under the permutation ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_orbit_size

        sage: perm_orbit_size(array('i', [0,3,1,2]), 2)
        3
        sage: perm_orbit_size(array('i', [0,3,1,2]), 0)
        1
    """
    if i < 0 or i >= len(p):
        raise ValueError
    cdef int j, s
    s = 1
    j = p.data.as_ints[i]
    while j != i:
        s += 1
        j = p.data.as_ints[j]
    return s


def perm_preimage(array.array p, int i):
    r"""
    Return the preimage of ``i`` under ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_init, perm_preimage

        sage: p = perm_init("(0,3,1,5)(2,4)")
        sage: perm_preimage(p, 3)
        0
        sage: perm_preimage(p, 2)
        4
    """
    cdef int j = i
    while p.data.as_ints[j] != i:
        j = p.data.as_ints[j]
    return j


def perm_on_list(array.array p, a, int n=-1, swap=None):
    r"""
    Inplace action of permutation on list-like objects.

    INPUT:

    - ``p`` - permutation

    - ``a`` - list, array

    - ``n`` - (optional) size of permutation

    - ``swap`` - (optional) a swap function

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import *

        sage: l = [0,1,2,3,4]
        sage: p = array('i', [4,2,3,0,1])
        sage: perm_on_list(p, l)
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
    if n == -1:
        n = len(p)
    cdef array.array seen = array.clone(p, n, True)
    cdef int i, j
    for i in range(n):
        if seen.data.as_ints[i]:
            continue
        seen.data.as_ints[i] = 1
        j = p.data.as_ints[i]
        while not seen.data.as_ints[j]:
            if swap:
                swap(a, i, j)
            else:
                tmp = a[i]
                a[i] = a[j]
                a[j] = tmp
            seen.data.as_ints[j] = 1
            j = p.data.as_ints[j]

#####################################################################
# Group operations
#####################################################################

def perm_invert(array.array p, int n=-1):
    r"""
    Return the inverse of the permutation ``l``.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_invert

        sage: perm_invert(array('i', [0, 3, 1, 2]))
        array('i', [0, 2, 3, 1])

        sage: perm_invert(array('i', [2, -1, 5, 0, -1, 3]))
        array('i', [3, -1, 0, 5, -1, 2])
    """
    if n == -1:
        n = len(p)
    cdef array.array res = array.clone(p, n, False)
    cdef int i
    for i in range(n):
        if p.data.as_ints[i] == -1:
            res.data.as_ints[i] = -1
        else:
            res.data.as_ints[p.data.as_ints[i]] = i
    return res


def perm_compose(array.array p1, array.array p2, int n=-1):
    r"""
    Return the product ``p1 p2``.

    In the notation ``p1 p2`` we use the right action, in other words
    ``p1`` is applied first.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_compose

        sage: perm_compose(array('i', [0,2,1]), array('i', [0,2,1]))
        array('i', [0, 1, 2])
        sage: perm_compose(array('i', [-1,2,3,1]), array('i', [-1,2,1,3]))
        array('i', [-1, 1, 3, 2])

        sage: perm_compose(array('i', [1,0,2,-1,-1]), array('i', [2,1,0,-1]), 3)
        array('i', [1, 2, 0])
    """
    if n == -1:
        n = len(p1)
    cdef array.array r = array.clone(p1, n, False)
    cdef int i
    for i in range(n):
        if p1.data.as_ints[i] == -1:
            r.data.as_ints[i] = -1
        else:
            r.data.as_ints[i] = p2.data.as_ints[p1.data.as_ints[i]]
    return r


perm_compose_00 = perm_compose


# TODO: do something less stupid
# (do it for each cycle independently, detecting if needed the period)
def perm_pow(array.array p, int k, int n=-1):
    r"""
    Return the power of the permutation ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_pow

        sage: perm_pow(array('i', [3, 0, 1, 2]), 2)
        array('i', [2, 3, 0, 1])
        sage: perm_pow(array('i', [3, 0, 1, 2]), -1)
        array('i', [1, 2, 3, 0])
    """
    if n == -1:
        n = len(p)
    if k == 0:
        return perm_id(n)

    cdef array.array q
    if k < 0:
        p = perm_invert(p, n)
        k = -k

    q = array.copy(p)
    k -= 1
    while k:
        q = perm_compose(q, p)
        k -= 1
    return q


def perm_compose_10(array.array p1, array.array p2, int n=-1):
    r"""
    Return the product ``p1^(-1) p2``

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_compose, perm_invert, perm_compose_10

        sage: p1 = array('i', [0,5,2,1,3,4])
        sage: p2 = array('i', [3,1,5,4,2,0])
        sage: perm_compose_10(p1, p2) == perm_compose(perm_invert(p1), p2)
        True
        sage: shuffle(p1)
        sage: shuffle(p2)
        sage: perm_compose_10(p1, p2) == perm_compose(perm_invert(p1), p2)
        True
    """
    if n == -1:
        n = len(p1)
    cdef int i
    cdef array.array r = array.clone(p1, n, False)
    for i in range(n):
        r.data.as_ints[p1.data.as_ints[i]] = p2.data.as_ints[i]
    return r


def perm_compose_01(array.array p1, array.array p2, int n=-1):
    r"""
    Return the product ``p1 p2^(-1)``

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_compose, perm_invert, perm_compose_01

        sage: p1 = array('i', [0,5,2,1,3,4])
        sage: p2 = array('i', [3,1,5,4,2,0])
        sage: perm_compose_01(p1, p2) == perm_compose(p1, perm_invert(p2)) # not tested
        True
        sage: shuffle(p1)
        sage: shuffle(p2)
        sage: perm_compose_01(p1, p2) == perm_compose(p1, perm_invert(p2)) # not tested
        True
    """
    raise NotImplementedError


def perm_compose_11(array.array p1, array.array p2, int n=-1):
    r"""
    Return the product `p1^(-1) p2^(-1)`.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import perm_compose, perm_invert, perm_compose_11

        sage: p1 = array('i', [0,5,2,1,3,4])
        sage: p2 = array('i', [3,1,5,4,2,0])
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
        ....:         q1 = perm_compose_11(array('i', p1), array('i', p2))
        ....:         q2 = perm_compose(perm_invert(array('i', p1)), perm_invert(array('i', p2)))
        ....:         assert q1 == q2, (p1, p2, q1, q2)
    """
    if n == -1:
        n = len(p1)
    cdef array.array r = array.clone(p1, n, False)
    for i in range(n):
        r.data.as_ints[p1.data.as_ints[p2.data.as_ints[i]]] = i
    return r


def perm_conjugate(array.array p1, array.array p2, int n=-1):
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
        sage: all(res[p2[i]] == p2[p1[i]] for i in range(23))
        True
    """
    if n == -1:
        n = len(p1)
    cdef array.array res = array.clone(p1, n, False)
    cdef int i
    for i in range(n):
        res.data.as_ints[p2.data.as_ints[i]] = p2.data.as_ints[p1.data.as_ints[i]]
    return res

#####################################################################
# Transitivity test
#####################################################################

def perms_transitive_components(p, int n=-1):
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
    if n == -1:
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


def perms_are_transitive(p, int n=-1):
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
    if n == -1:
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

        sage: from veerer.permutation import perms_relabel
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

        sage: from veerer.permutation import perms_canonical_labels_from
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

#####################################################################
# Triangulation relabellings
#####################################################################

def triangulation_relabelling_from(array.array vp, array.array ep, int start_edge):
    r"""
    Return a canonical relabelling where ``start_edge`` is mapped to ``0``.

    EXAMPLES::

        sage: from array import array
        sage: from veerer.permutation import triangulation_relabelling_from
        sage: vp = array('i', [9, 13, 12, 11, 10, 14, 5, 6, 7, 8, 29, 25, 26, 27, 28, 1, 15, 16, 17, 18, 19, 0, 4, 3, 2, 21, 22, 23, 24, 20])
        sage: ep = array('i', [29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0])
        sage: triangulation_relabelling_from(vp, ep, 0)
        array('i', [0, 7, 21, 10, 18, 14, 20, 13, 17, 28, 27, 26, 25, 24, 23, 6, 5, 4, 3, 2, 1, 12, 16, 9, 15, 11, 19, 8, 22, 29])
    """
    # NOTE: the algorithm is as follows
    # 0) we set k=0 and m=n-1 (labelling counter), any time we choose a new
    #    label it is k and k is incremented. If it is not a folded edge, its
    #    twin gets labelled m and m is decremented.
    # 1) the start_edge is relabelled 0 and set as pending
    # 2) while there is a pending half-edge, we look at its vp-cycle and
    #    relabel each edge if needed
    cdef int n = len(vp)

    cdef int k = 0      # current available label at the front.
    cdef int m = n - 1  # current available label at the back.
    cdef array.array relabelling = array.clone(vp, n, False)
    cdef int i
    for i in range(n):
        relabelling.data.as_ints[i] = -1

    relabelling.data.as_ints[start_edge] = 0
    k += 1

    if ep.data.as_ints[start_edge] != start_edge:
        relabelling.data.as_ints[ep.data.as_ints[start_edge]] = m
        m -= 1

    cdef array.array to_process = array.clone(vp, n, False)
    cdef int s = 1
    to_process.data.as_ints[0] = start_edge
    if ep.data.as_ints[start_edge] != start_edge:
        to_process.data.as_ints[1] = ep.data.as_ints[start_edge]
        s = 2

    cdef int e, e0
    while s:
        s -= 1
        e0 = to_process.data.as_ints[s]
        e = vp.data.as_ints[e0]
        while e != e0:
            if relabelling.data.as_ints[e] == -1:
                relabelling.data.as_ints[e] = k
                k += 1
                if ep.data.as_ints[e] != e:
                    relabelling.data.as_ints[ep.data.as_ints[e]] = m
                    m -= 1
                    to_process.data.as_ints[s] = ep.data.as_ints[e]
                    s += 1
            e = vp.data.as_ints[e]

    return relabelling
