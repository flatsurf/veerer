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

import random
from veerer import VeeringTriangulation, VeeringFlipSequence

def test_Q_p4():
    T = VeeringTriangulation("(0,1,2)", "PBR")
    B = VeeringFlipSequence(T, "0B", "(0,1)")
    R = VeeringFlipSequence(T, "0R", "(0,2)")

    for f in [B * R, R * B, B * B * R,
              B * R * B, R * B * B,
              B * R * R, R * B * R,
              R * R * B]:
        assert f.is_pseudo_anosov(), f
        g = f.inverse()
        assert g.is_pseudo_anosov(), (f, g)
        h = g.inverse()
        assert f.is_identical(h), (f, h)

def test_H_01():
    T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "PBR")
    B = VeeringFlipSequence(T, "0B", "(1,0,~1,~0)(2,~2)")
    R = VeeringFlipSequence(T, "0R", "(0,2)(1,~1)")

    for f in [B * R, R * B, B * B * R,
              B * R * B, R * B * B,
              B * R * R, R * B * R,
              R * R * B]:
        assert f.is_pseudo_anosov(), f
        g = f.inverse()
        assert g.is_pseudo_anosov(), (f, g)
        h = g.inverse()
        assert f.is_identical(h), (f, h)

def test_Q_11_p5():
    Vc = VeeringTriangulation("(0,~5,4)(3,5,6)(1,2,~6)", "PPBPRBR")
    Vr = VeeringTriangulation("(0,6,5)(1,2,~6)(3,4,~5)", "BPBBRPR")
    Vl = VeeringTriangulation("(0,~5,4)(1,6,5)(2,3,~6)", "PRBRRBP")

    CR5 = VeeringFlipSequence(Vc, "1B", "(1,2)")
    CL5 = VeeringFlipSequence(Vc, "0R", "(0,4)")
    R3 = VeeringFlipSequence(Vc, "0B 3B", "(0,3)")
    R5 = VeeringFlipSequence(Vr, "1B", "(1,2)")
    R1 = VeeringFlipSequence(Vr, "1R 5R", "(0,2,3)(1,4)(5,6)")
    R2 = VeeringFlipSequence(Vr, "5B")
    L3 = VeeringFlipSequence(Vc, "1R 3R", "(1,3)")
    L5 = VeeringFlipSequence(Vl, "0R", "(0,4)")
    L1 = VeeringFlipSequence(Vl, "0B 6B", "(0,2)(1,4,3)(5,6,~5,~6)")
    L2 = VeeringFlipSequence(Vl, "6R", "(6,~6)")

    for f in [L1 * L5, L5 * L1,
              R1 * R5, R5 * R1,
              L1 * L1 * L5, L1 * L5 * L1, L5 * L1 * L1,
              R1 * R1 * R5, R1 * R5 * R1, R5 * R1 * R1,
              L3 * L5 * L2 * R3 * R5 * R2,
              L3 * L5**2 * L2 * R3 * R5 * R2,
              L3 * L5 * L2 * R3 * R5**2 * R2]:
        assert f.is_pseudo_anosov(), f
        g = f.inverse()
        assert g.is_pseudo_anosov(), (f, g)
        h = g.inverse()
        assert f.is_identical(h), (f, h)

if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
