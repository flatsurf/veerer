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
 ("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB", 200),
 ("(0,8,~7)(1,9,~8)(2,11,~10)(3,10,~9)(4,~1,~11)(5,~2,~4)(6,~3,~5)(7,~0,~6)", "RRRRBBBBBBBB", 200),
 ("(0,8,~7)(1,9,~8)(2,10,~9)(3,~6,7)(4,~2,~11)(5,~3,~4)(6,~0,~5)(11,~10,~1)", "RRRRBBBBBBBB", 200),
 ("(0,10,~9)(1,11,~10)(2,12,~11)(3,13,~12)(4,~7,8)(5,~3,~14)(6,~4,~5)(7,~0,~6)(9,~2,~8)(14,~13,~1)", "RRRRRBBBBBBBBBB", 200),
 ("(0,12,~11)(1,13,~12)(2,14,~13)(3,15,~14)(4,17,~16)(5,~10,11)(6,~3,~17)(7,~2,~6)(8,~5,~7)(9,~0,~8)(10,~4,~9)(16,~15,~1)", "RRRRRRBBBBBBBBBBBB", 200),
 ("(0,1,2)", "RRB", 10),
 ("(0,2,3)(1,4,~0)(5,6,~1)", "BRRBBBB", 100),
 ("(0,3,4)(1,~3,5)(2,6,~4)", "BRBBRRB", 100)
 ])
def test_flip_reduced(fp, cols, repeat):
    V0 = VeeringTriangulation(fp, cols)
    assert V0.is_core()

    for _ in range(repeat):
        V = V0.copy(mutable=True)
        e = random.choice(V.forward_flippable_edges())
        col = random.choice([RED, BLUE])

        # test if this is a valid choice and change if not
        oldcol = V.edge_colour(e)
        V.flip(e, col)
        if not V.edge_has_curve(e):
            col = BLUE if col == RED else RED
            V.flip_back(e, oldcol)
            V.flip(e, col)

        # square-tiled implies cylindrical
        assert (not V.is_square_tiled(BLUE) or not V.is_square_tiled(RED)) or V.is_cylindrical()

        # Build a purple, blue and red versions of the triangulation
        P = V.copy()
        P.forgot_forward_flippable_colour()
        R = P.copy()
        R.set_colour(RED)
        B = P.copy()
        B.set_colour(BLUE)

        # test that cylinders in the PURPLE surface coincide with
        # the BLUE/RED versions
        assert P.cylinders(RED) == R.cylinders(RED), P
        assert P.cylinders(BLUE) == B.cylinders(BLUE), P

        # test that cylindrical property coicncide
        assert P.is_cylindrical(RED) == R.is_cylindrical(RED), P
        assert P.is_cylindrical(BLUE) == B.is_cylindrical(BLUE), P
        assert P.is_cylindrical(PURPLE) == (B.is_cylindrical(BLUE) and R.is_cylindrical(RED)), P

        # being purple cylindrical is equivalent to square-tiled
        assert P.is_cylindrical(PURPLE) == (B.is_square_tiled(BLUE) and R.is_square_tiled(RED)), P


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
