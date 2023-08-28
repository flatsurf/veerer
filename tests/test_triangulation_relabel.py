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

from veerer.permutation import perm_init
from veerer.triangulation import Triangulation

def test_relabel():
    T = Triangulation([(0,1,2), (-1,-2,-3)], mutable=True)
    p = perm_init([1,5,0,2,4,3])
    T.relabel(p)
    assert T.faces() == [[0, 1, 5], [2, 3, 4]]
    assert T.edges() == [[0, 2], [1, 3], [4, 5]]
    T._check()

if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
