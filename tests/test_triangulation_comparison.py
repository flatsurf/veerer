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

from veerer import Triangulation, VeeringTriangulation

def test_triangulation_comparison():
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

    for U in [T1, T2, T3, T4]:
        for V in [T1, T2, T3, T4]:
            assert (U < V) == (V > U)
            assert (U <= V) == (V >= U)
            assert (U < V) == (not (U >= V))
            assert (U > V) == (not (U <= V))
            assert (U <= V) == ((U < V) or (U == V))
            assert (U >= V) == ((U > V) or (U == V))

def test_veering_triangulation_comparison():
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

    for U in [V1, V2, V3, V4]:
        for V in [V1, V2, V3, V4]:
            assert (U < V) == (V > U)
            assert (U <= V) == (V >= U)
            assert (U < V) == (not (U >= V))
            assert (U > V) == (not (U <= V))
            assert (U <= V) == ((U < V) or (U == V))
            assert (U >= V) == ((U > V) or (U == V))


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
