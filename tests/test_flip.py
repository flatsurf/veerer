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

import pytest

pytest.importorskip('ppl')

import random
from veerer import VeeringTriangulation, RED, BLUE, PURPLE, VERTICAL

@pytest.mark.parametrize("fp, cols, repeat",
[("(0,8,~7)(1,11,~10)(2,10,~9)(3,9,~8)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 50),
 ("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB", 50),
 ("(0,8,~7)(1,9,~8)(2,11,~10)(3,10,~9)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 50),
 ("(0,8,~7)(1,9,~8)(2,10,~9)(3,~6,7)(4,~2,~11)(5,~3,~4)(6,~0,~5)(11,~10,~1)", "RRRRBBBBBBBB", 50),
 ("(0,10,~9)(1,11,~10)(2,12,~11)(3,13,~12)(4,~7,8)(5,~3,~14)(6,~4,~5)(7,~0,~6)(9,~2,~8)(14,~13,~1)", "RRRRRBBBBBBBBBB", 50),
 ("(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)", "RRRRRRBBBBBBBBBBBB", 50),
 ("(0,1,2)", "RRB", 10),
 ("(0,2,3)(1,4,~0)(5,6,~1)", "BRRBBBB", 50),
 ("(0,3,4)(1,~3,5)(2,6,~4)", "BRBBRRB", 50)
 ])
def test_flip(fp, cols, repeat):
    V = VeeringTriangulation(fp, cols)
    assert V.is_core()

    d = V.stratum_dimension()

    for _ in range(repeat):
        e = random.choice(V.forward_flippable_edges())
        oldcol = V.edge_colour(e)
        col = random.choice([RED, BLUE])

        W = V.copy(mutable=True)
        W.flip(e, col)
        assert W.is_backward_flippable(e)
        test1 = W.edge_has_curve(e)
        test2 = W.is_core()
        assert test1 == test2, (V, W, e, col)

        X = W.copy(mutable=True)
        X.forgot_forward_flippable_colour()
        test3 = X.train_track_polytope(VERTICAL).affine_dimension() == X.stratum_dimension()
        assert test1 == test3, (W, X)

        Y = V.copy(mutable=True)
        Y.forgot_forward_flippable_colour()
        Y.flip(e, col, reduced=False)
        assert Y.is_backward_flippable(e)

        if test1:
            V = W
        else:
            W.flip_back(e, oldcol)
            assert V == W

@pytest.mark.parametrize("fp, cols, repeat",
[("(0,8,~7)(1,11,~10)(2,10,~9)(3,9,~8)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 50),
 ("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB", 50),
 ("(0,8,~7)(1,9,~8)(2,11,~10)(3,10,~9)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 50),
 ("(0,8,~7)(1,9,~8)(2,10,~9)(3,~6,7)(4,~2,~11)(5,~3,~4)(6,~0,~5)(11,~10,~1)", "RRRRBBBBBBBB", 50),
 ("(0,10,~9)(1,11,~10)(2,12,~11)(3,13,~12)(4,~7,8)(5,~3,~14)(6,~4,~5)(7,~0,~6)(9,~2,~8)(14,~13,~1)", "RRRRRBBBBBBBBBB", 50),
 ("(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)", "RRRRRRBBBBBBBBBBBB", 50),
 ("(0,1,2)", "RRB", 10),
 ("(0,2,3)(1,4,~0)(5,6,~1)", "BRRBBBB", 50),
 ("(0,3,4)(1,~3,5)(2,6,~4)", "BRBBRRB", 50)
 ])
def test_flip_reduced(fp, cols, repeat):
    V0 = VeeringTriangulation(fp, cols, mutable=True)
    assert V0.is_core()
    V0.forgot_forward_flippable_colour()
    V0.set_immutable()

    for _ in range(repeat):
        V = V0.copy(mutable=True)
        e = random.choice(V.forward_flippable_edges())
        assert V.edge_colour(e) == PURPLE
        col = random.choice([RED, BLUE])

        # test if this is a valid choice and change if not
        V.flip(e, col, reduced=False)
        if not V.edge_has_curve(e):
            col = BLUE if col == RED else RED
        V.flip_back(e, PURPLE)

        # actually do the flip
        V.flip(e, col, reduced=True)
        assert V.edge_has_curve(e)
        assert V.forward_flippable_edges() == [e[0] for e in V.edges() if V.edge_colour(e[0]) == PURPLE], (V0, e, col)

if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
