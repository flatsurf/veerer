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
from veerer.constants import RED, BLUE
from veerer.triangulation import Triangulation
from veerer.veering_triangulation import VeeringTriangulation
from veerer.flip_sequence import VeeringFlipSequence

@pytest.mark.parametrize("fp, repeat",
[("(0,1,2)(~0,~1,~2)", 20),
 ("(0,1,2)(~0,~2,3)", 20),
 ("(0,1,~3)(~0,2,3)", 20),
 ("(1,~3,2)(~2,0,3)", 20),
 ("(0,1,2)", 20),
 ("(0,~5,4)(3,5,6)(1,2,~6)", 30)])
def test_triangulation_reconstruction(fp, repeat):
    T = Triangulation(fp)
    for _ in range(repeat):
        assert eval(repr(T)) == T
        e = random.choice(T.flippable_edges())
        T.flip(e)

@pytest.mark.parametrize("fp, cols, repeat",
[("(0,1,2)(~0,~1,~2)", "BBR", 30),
 ("(0,1,2)(~0,~1,~2)", "PBR", 30),
 ("(0,1,2)", "BBR", 30),
 ("(0,1,2)", "PBR", 30),
 ("(0,~5,4)(3,5,6)(1,2,~6)", "PPBPRBR", 30)])
def test_veering_triangulation_reconstruction(fp, cols, repeat):
    V = VeeringTriangulation(fp, cols)
    for _ in range(repeat):
        assert eval(repr(V)) == V, V
        e = random.choice(V.forward_flippable_edges())
        col = random.choice([RED, BLUE])

        # test flip validity
        oldcol = V.edge_colour(e)
        V.flip(e, col, reduced=False)
        if not V.edge_has_curve(e):
            col = BLUE if col == RED else RED
            V.flip_back(e, oldcol)

        # do flip
        V.flip(e, col)

@pytest.mark.parametrize("fp, cols, repeat",
[("(0,1,2)(~0,~1,~2)", "BBR", 10),
 ("(0,1,2)(~0,~1,~2)", "PBR", 10),
 ("(0,1,2)", "BBR", 10),
 ("(0,1,2)", "PBR", 10),
 ("(0,~5,4)(3,5,6)(1,2,~6)", "PPBPRBR", 10)])
def test_veering_flip_sequence_reconstruction(fp, cols, repeat):
    V = VeeringTriangulation(fp, cols)
    for _ in range(repeat):
        F = V.random_forward_flip_sequence(random.randint(0,10), relabel=True)
        assert eval(repr(F)) == F, F
        V = F.end()

if __name__ == '__main__': sys.exit(pytest.main(sys.argv))
