######################################################################
# This file is part of veerer.
#
#       Copyright (C) 2020 Vincent Delecroix
#
# veerer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# veerer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with veerer. If not, see <https://www.gnu.org/licenses/>.
######################################################################

import sys
import pytest

import random
from veerer.permutation import perm_init
from veerer.triangulation import Triangulation

def test_reconstruction():
    T = Triangulation("(0,~5,4)(3,5,6)(1,2,~6)")
    T.flip(0)
    T.flip(3)
    T.relabel("(0,3)")
    assert eval(repr(T)) == T

@pytest.mark.parametrize("fp, repeat",
[("(0,1,2)(~0,~1,~2)", 20),
 ("(0,1,2)(~0,~2,3)", 20),
 ("(0,1,~3)(~0,2,3)", 20),
 ("(1,~3,2)(~2,0,3)", 20),
 ("(0,1,2)", 20),
 ("(0,~5,4)(3,5,6)(1,2,~6)", 30)])
def test_reconstruction_random(fp, repeat):
    T = Triangulation(fp)
    for _ in range(repeat):
        e = random.choice(T.flippable_edges())
        T.flip(e)
        assert eval(repr(T)) == T

if __name__ == '__main__': sys.exit(pytest.main(sys.argv))
