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

from veerer import VeeringTriangulation
from veerer.permutation import perm_random_centralizer

@pytest.mark.parametrize("fp, cols, repeat",
    [("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)", "BRBBBRRBBBBR", 10),
     ("(0,~4,3)(1,~9,~6)(2,~3,6)(4,~8,~5)(5,7,11)(8,10,~2)(9,~10,~1)(~11,~0,~7)", "PBBRBPPBRRPR", 10),
     ("(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)", "BRRBRR", 10),
     ("(0,1,2)", "RRB", 10),
     ("(0,1,2)", "RPB", 10)])
def test_isomorphism(fp, cols, repeat):
    V = VeeringTriangulation(fp, cols)
    for _ in range(repeat):
        W = V.copy()
        p = perm_random_centralizer(V.edge_permutation(copy=False))
        W.relabel(p)
        assert V.is_isomorphic(W) is True
        ans, cert = V.is_isomorphic(W, True)
        assert ans is True
        V.relabel(cert)
        assert V == W

if __name__ == '__main__': sys.exit(pytest.main(sys.argv))
