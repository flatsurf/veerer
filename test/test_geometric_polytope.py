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

pytest.importorskip('ppl')

import random
from veerer import VeeringTriangulation

@pytest.mark.parametrize("fp, cols, repeat",
[("(0,8,~7)(1,11,~10)(2,10,~9)(3,9,~8)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 20),
 ("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB", 20),
 ("(0,8,~7)(1,9,~8)(2,11,~10)(3,10,~9)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 20),
 ("(0,8,~7)(1,9,~8)(2,10,~9)(3,~6,7)(4,~2,~11)(5,~3,~4)(6,~0,~5)(11,~10,~1)", "RRRRBBBBBBBB", 20),
 ("(0,10,~9)(1,11,~10)(2,12,~11)(3,13,~12)(4,~7,8)(5,~3,~14)(6,~4,~5)(7,~0,~6)(9,~2,~8)(14,~13,~1)", "RRRRRBBBBBBBBBB", 20),
 ("(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)", "RRRRRRBBBBBBBBBBBB", 20),
 ("(0,1,2)", "RRB", 20),
 ("(0,2,3)(1,4,~0)(5,6,~1)", "BRRBBBB", 20)])
def test_geometric_polytope(fp, cols, repeat):
    V = VeeringTriangulation(fp, cols)

    for _ in range(repeat):
        V.random_forward_flip()
        test1 = V.is_geometric('LP')
        test2 = V.is_geometric('polytope')
        test3 = V.is_geometric('polytope2')
        assert test1 == test2 == test3, V.to_string()

@pytest.mark.parametrize("fp, cols, repeat",
[("(0,8,~7)(1,11,~10)(2,10,~9)(3,9,~8)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 50),
 ("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB", 50),
 ("(0,8,~7)(1,9,~8)(2,11,~10)(3,10,~9)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 50),
 ("(0,8,~7)(1,9,~8)(2,10,~9)(3,~6,7)(4,~2,~11)(5,~3,~4)(6,~0,~5)(11,~10,~1)", "RRRRBBBBBBBB", 50),
 ("(0,10,~9)(1,11,~10)(2,12,~11)(3,13,~12)(4,~7,8)(5,~3,~14)(6,~4,~5)(7,~0,~6)(9,~2,~8)(14,~13,~1)", "RRRRRBBBBBBBBBB", 50),
 ("(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)", "RRRRRRBBBBBBBBBBBB", 50),
 ("(0,1,2)", "RRB", 10),
 ("(0,2,3)(1,4,~0)(5,6,~1)", "BRRBBBB", 50)])
def test_geometric_flips(fp, cols, repeat):
    V = VeeringTriangulation(fp, cols)

    while not V.is_geometric():
        V.random_forward_flip()

    for _ in range(repeat):
        flips = random.choice(V.geometric_flips())
        for e, col in flips:
            V.flip(e, col)
        assert V.is_geometric()

if __name__ == '__main__': sys.exit(pytest.main(sys.argv))
