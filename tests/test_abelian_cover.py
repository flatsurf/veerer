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
from veerer import VeeringTriangulation, RED, BLUE, PURPLE, GREEN

@pytest.mark.parametrize("fp, cols",
[
 ("(0,1,2)", "RRB"),
 ("(0,1,2)", "BBR"),
 ("(0,1,2)", "PBR"),
 ("(0,1,2)", "GRB"),
 ("(0,1,2)", "RPG"),
 ("(0,1,2)", "BGP")
])
def test_triangle(fp, cols):
    V = VeeringTriangulation(fp, cols)
    A = V.abelian_cover()
    assert A.angles() == [2], (V, A)

if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
