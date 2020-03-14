######################################################################
# This file is part of veering.
#
#       Copyright (C) 2020 Vincent Delecroix
#
# flatsurf is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# flatsurf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with veerer. If not, see <https://www.gnu.org/licenses/>.
######################################################################

import sys
import pytest

import itertools
from random import randint
from veerer.permutation import (perm_init, perm_invert, perm_compose, perm_id,
        perm_random, perm_check, perm_cycles, perm_from_cycles)

def test_zero():
    p = perm_init([])
    assert perm_compose(p, p) == p
    assert perm_invert(p) == p

@pytest.mark.parametrize("n", [0, 1, 2, 3, 4])
def test_init_invert_compose(n):
        q = perm_id(n)
        for p in itertools.permutations(range(n)):
            assert perm_compose(perm_invert(p), p) == q

def test_random(repeat=10):
    for i in range(repeat):
        assert perm_check(perm_random(i), i)

def test_cycles(repeat=50):
    for i in range(repeat):
        p = perm_random(randint(0, 10))
        c = perm_cycles(p)
        assert perm_from_cycles(c) == p

if __name__ == '__main__': sys.exit(pytest.main(sys.argv))
