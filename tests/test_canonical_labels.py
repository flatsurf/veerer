######################################################################
# This file is part of veerer.
#
#       Copyright (C) 2024 Vincent Delecroix
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

from veerer.permutation import perm_init, perm_random, perm_random_centralizer
from veerer.triangulation import Triangulation
from veerer.strebel_graph import StrebelGraph
from veerer.veering_triangulation import VeeringTriangulation
from veerer.linear_family import VeeringTriangulationLinearFamily, StrebelGraphLinearFamily

@pytest.mark.parametrize("C", [
    Triangulation("(0,1,2)(~2,~0,~1)"),
    Triangulation("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"),
    Triangulation("(0,10,~9)(1,~8,9)(2,~6,7)(3,~4,5)(4,~3,~11)(6,~2,~5)(8,~1,~7)(11,~10,~0)"),
    Triangulation("(0,1,2)", boundary="(~2:1)(~1:1)(~0:1)"),
    Triangulation("(0,2,1)(3,~1,~0)", boundary="(~3:1,~2:1)"),
    VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB"),
    VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBPBRBPRB"),
    VeeringTriangulation("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)", "BRBBBRRBBBBR"),
    VeeringTriangulation("(0,1,2)(3,4,~0)", boundary="(~4:1)(~3:1)(~2:1)(~1:1)", colouring="RBRBR"),
    StrebelGraph("(0,1,2)(~0,~1:1,~2:2)"),
    StrebelGraph("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"),
    VeeringTriangulationLinearFamily("(0,1,2)(~2,~0,~1)", "RRB", [(1, 0, -1), (0, 1, 1)]),
    VeeringTriangulationLinearFamily("(0,1,2)(3,4,5)(~5,~3,~1)(~4,~2,~0)", "BRRBRR", [(1, 0, 1, 0, 0, 0), (0, 1, 1, 0, 1, 1), (0, 0, 0, 1, 0, 1)]),
    StrebelGraphLinearFamily("(0,1,2)(~2:2,~0,~1:1)", [(1, 0, 1), (0, 1, -1)])
])
def test_canonical_labels(C):
    C = C.copy(mutable=True)
    C.set_canonical_labels()
    C.set_immutable()
    D = C.copy(mutable=True)
    for _ in range(20):
        p = perm_random_centralizer(C.edge_permutation(copy=False))
        D.relabel(p)
        D.set_canonical_labels()
        assert C == D

if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
