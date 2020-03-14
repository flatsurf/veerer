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

import random
from veerer import VeeringTriangulation, RED, BLUE, VERTICAL

@pytest.mark.parametrize("fp, cols, repeat",
[("(0,8,~7)(1,11,~10)(2,10,~9)(3,9,~8)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 50),
 ("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB", 50),
 ("(0,8,~7)(1,9,~8)(2,11,~10)(3,10,~9)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 50),
 ("(0,8,~7)(1,9,~8)(2,10,~9)(3,~6,7)(4,~2,~11)(5,~3,~4)(6,~0,~5)(11,~10,~1)", "RRRRBBBBBBBB", 50),
 ("(0,10,~9)(1,11,~10)(2,12,~11)(3,13,~12)(4,~7,8)(5,~3,~14)(6,~4,~5)(7,~0,~6)(9,~2,~8)(14,~13,~1)", "RRRRRBBBBBBBBBB", 50),
 ("(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)", "RRRRRRBBBBBBBBBBBB", 50),
 ("(0,1,2)", "RRB", 10),
 ("(0,2,3)(1,4,~0)(5,6,~1)", "BRRBBBB", 50),
 ])
def test_flip(fp, cols, repeat):
    V = VeeringTriangulation(fp, cols)
    assert V.is_core()

    if V.is_abelian():
        d = 2*V.genus() + V.num_vertices() - 1
    else:
        d = 2*V.genus() + V.num_vertices() - 2

    for _ in range(repeat):
        e = random.choice(V.forward_flippable_edges())
        col = random.choice([RED, BLUE])
        W = V.copy()
        W.flip(e, col)
        test1 = W.edge_has_curve(e)
        test2 = W.is_core()
        assert test1 == test2, (V, W, e, col)

        X = W.copy()
        X.forgot_forward_flippable_color()
        test3 = X.train_track_polytope(VERTICAL).affine_dimension() == X.stratum_dimension()
        assert test1 == test3, (W, X)

        if test1:
            V = W

if __name__ == '__main__': sys.exit(pytest.main(sys.argv))
