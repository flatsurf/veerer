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
from veerer.permutation import perm_random_centralizer

def test_non_isomorphic():
     T1 = Triangulation("(0,1,2)")
     T2 = Triangulation("(0,~0,1)(~1,2,~2)")
     T3 = Triangulation("(0,~0,1)(~1,2,3)")
     T4 = Triangulation("(0,1,2)(~2,3,4)")
     T5 = Triangulation("(0,1,2)(~0,~1,~2)")
     T6 = Triangulation("(0,1,2)(~0,~2,~1)")
     T = [T1, T2, T3, T4, T5, T6]
     for S1 in T:
         for S2 in T:
             if S1 == S2:
                 continue
             assert S1.is_isomorphic_to(S2) is False
             assert S1.is_isomorphic_to(S2,True) == (False, None)

def test_veering_non_isomorphic():
    V1 = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
    V2 = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "BBR")
    V3 = VeeringTriangulation("(0,1,2)", 'RRB')
    V = [V1, V2, V3]
    for W1 in V:
        for W2 in V:
            if W1 == W2:
                continue
            assert W1.is_isomorphic_to(W2) is False
            assert W1.is_isomorphic_to(W2, True) == (False, None)

@pytest.mark.parametrize("fp, repeat",
    [("(0,5,1)(~0,4,2)(~1,~2,~4)(3,6,~5)", 10),
     ("(0,1,2)", 10),
     ("(0,~0,1)(~1,2,3)", 10)])
def test_isomorphism(fp, repeat):
    T = Triangulation(fp)
    for _ in range(repeat):
        U = T.copy(mutable=True)
        p = perm_random_centralizer(U.edge_permutation(copy=False))
        U.relabel(p)
        assert U.is_isomorphic_to(T) is True
        ans, cert = U.is_isomorphic_to(T, True)
        assert ans is True
        U.relabel(cert)
        assert T == U

@pytest.mark.parametrize("fp, cols, repeat",
    [("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)", "BRBBBRRBBBBR", 10),
     ("(0,~4,3)(1,~9,~6)(2,~3,6)(4,~8,~5)(5,7,11)(8,10,~2)(9,~10,~1)(~11,~0,~7)", "PBBRBPPBRRPR", 10),
     ("(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)", "BRRBRR", 10),
     ("(0,1,2)", "RRB", 10),
     ("(0,1,2)", "RPB", 10)])
def test_veering_isomorphism(fp, cols, repeat):
    V = VeeringTriangulation(fp, cols)
    for _ in range(repeat):
        W = V.copy(mutable=True)
        p = perm_random_centralizer(W.edge_permutation(copy=False))
        W.relabel(p)
        assert W.is_isomorphic_to(V) is True
        ans, cert = W.is_isomorphic_to(V, True)
        assert ans is True
        W.relabel(cert)
        assert V == W

if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
