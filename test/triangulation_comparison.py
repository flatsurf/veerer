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

from veerer import Triangulation, VeeringTriangulation

def test_comparison():
    T1 = Triangulation("(0,1,2)(~0,~1,~2)")
    T2 = Triangulation("(0,1,2)(~0,~1,~2)")
    T3 = Triangulation("(0,1,2)(~0,~2,~1)")
    T4 = Triangulation("(0,1,2)")
    assert T1 == T2
    assert not (T1 == T3)
    assert not (T1 == T4)
    assert not (T3 == T4)

    assert not (T1 != T2)
    assert T1 != T3
    assert T1 != T4
    assert T3 != T4

    V1 = VeeringTriangulation("(0,1,2)", "RRB")
    V2 = VeeringTriangulation("(0,1,2)", "RRB")
    V3 = VeeringTriangulation("(0,1,2)", "RBB")
    V4 = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
    assert V1 == V2
    assert not (V1 == V3)
    assert not (V1 == V4)
    assert not (V3 == V4)

    assert not (V1 != V2)
    assert V1 != V3
    assert V1 != V4
    assert V3 != V4

if __name__ == '__main__': sys.exit(pytest.main(sys.argv))
