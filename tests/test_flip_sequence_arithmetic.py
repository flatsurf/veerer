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

from veerer.veering_triangulation import VeeringTriangulation
from veerer.flip_sequence import VeeringFlipSequence

def test_mul_torus():
    T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "PBR")
    B = VeeringFlipSequence(T, "0B", "(1,0,~1,~0)(2,~2)")
    R = VeeringFlipSequence(T, "0R", "(0,2)(1,~1)")

    assert B * (B * R) == (B * B) * R
    assert (B * R) * (R * B) == B * (R * R) * B
    assert (B * R) * (R * B) == B * (R * R * B)

def test_mul_Q_11_m5():
    Vc = VeeringTriangulation("(0,~5,4)(3,5,6)(1,2,~6)", "PPBPRBR")
    Vr = VeeringTriangulation("(0,6,5)(1,2,~6)(3,4,~5)", "BPBBRPR")

    CR5 = VeeringFlipSequence(Vc, "1B", "(1,2)")
    CL5 = VeeringFlipSequence(Vc, "0R", "(0,4)")
    R3 = VeeringFlipSequence(Vc, "0B 3B", "(0,3)")
    R5 = VeeringFlipSequence(Vr, "1B", "(1,2)")
    R1 = VeeringFlipSequence(Vr, "1R 5R", "(0,2,3)(1,4)(5,6)")
    R2 = VeeringFlipSequence(Vr, "5B")

    assert CR5 * CL5 * R3 * R5 * R1 * R5 * R2 == CR5 * (CL5 * (R3 * (R5 * (R1 * (R5 * R2)))))

def test_power_torus():
    T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "PBR")
    B = VeeringFlipSequence(T, "0B", "(1,0,~1,~0)(2,~2)")
    R = VeeringFlipSequence(T, "0R", "(0,2)(1,~1)")

    assert B**0 == R**0 == VeeringFlipSequence(T)

    assert B**1 == B
    assert R**1 == R

    assert B**2 == B*B
    assert R**2 == R*R

    assert B**3 == B*B*B
    assert R**3 == R*R*R

    assert (B*R)**2 == B*R*B*R
    A = B*R*B
    assert (A**2)**3 == (A**3)**2

def test_power_Q_11_m5():
    Vc = VeeringTriangulation("(0,~5,4)(3,5,6)(1,2,~6)", "PPBPRBR")
    Vr = VeeringTriangulation("(0,6,5)(1,2,~6)(3,4,~5)", "BPBBRPR")

    CR5 = VeeringFlipSequence(Vc, "1B", "(1,2)")
    CL5 = VeeringFlipSequence(Vc, "0R", "(0,4)")
    R3 = VeeringFlipSequence(Vc, "0B 3B", "(0,3)")
    R5 = VeeringFlipSequence(Vr, "1B", "(1,2)")
    R1 = VeeringFlipSequence(Vr, "1R 5R", "(0,2,3)(1,4)(5,6)")
    R2 = VeeringFlipSequence(Vr, "5B")

    A = R1 * R2 * R3
    assert A**2 == A*A
    assert A**3 == A*A*A

if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
