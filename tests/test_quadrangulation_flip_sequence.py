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
from veerer.permutation import perm_random
from veerer.veering_quadrangulation import VeeringQuadrangulation, VeeringQuadrangulationFlipSequence

def random_flip_sequence(pr, pl, min_length=1, relabel=True, loop=False):
    V = VeeringQuadrangulation(pr, pl)
    n = V._n
    fp = VeeringQuadrangulationFlipSequence(V)
    for _ in range(min_length):
        representatives, sizes = fp._end.r_cycles_representatives()
        mult = [random.randint(0, 2*s) for s in sizes]
        while not any(mult):
            mult = [random.randint(0, 2*s) for s in sizes]
        fp.append_flip(representatives, mult)

    if relabel:
        fp.append_relabelling(perm_random(n))

    if loop:
        p = fp.find_closure()
        while p is None:
            representatives, sizes = fp._end.r_cycles_representatives()
            mult = [random.randint(0, s) for s in sizes]
            while not any(mult):
                mult = [random.randint(0, s) for s in sizes]
            fp.append_flip(representatives, mult)
            p = fp.find_closure()
        fp.append_relabelling(p)
        assert fp.is_loop()

    fp._check()
    return fp

@pytest.mark.parametrize("pr, pl, repeat",
    [("(0,1)(2)", "(0,2)(1)", 20),
     ("(0,1,2,3)(4,5,6,7)", "(0,4,3)", 20),
     ("(0,1,2)(3,4,5)(6,7,8)", "(0,1)(2,3)(4,6)(5,7,8)", 20)])
def test_append_relabelling(pr, pl, repeat):
    V = VeeringQuadrangulation(pr, pl)
    n = V._n
    fp = VeeringQuadrangulationFlipSequence(V)
    for _ in range(repeat):
        representatives, sizes = fp._end.r_cycles_representatives()
        mult = [random.randint(1,s) for s in sizes]
        fp.append_flip(representatives, mult)
        fp.append_relabelling(perm_random(n))
        fp._check()

@pytest.mark.parametrize("pr, pl, repeat",
    [("(0,1)(2)", "(0,2)(1)", 20),
     ("(0,3,2,1)", "(0,3)(1,2)", 20),
     ("(0,1,2,3,4)", "(0)(1,4)(2,3)", 20)])
def test_relabel(pr, pl, repeat):
    pytest.importorskip('surface_dynamics')
    for _ in range(repeat):
        fp = random_flip_sequence(pr, pl, loop=True)
        n = fp._start._n
        fp2 = fp.copy()
        fp2.relabel(perm_random(n))

        fp2._check()
        assert fp2.is_loop()
        assert fp.matrix().charpoly() == fp2.matrix().charpoly()
        assert fp.is_relabelling_equivalent(fp2)

@pytest.mark.parametrize("pr, pl, repeat",
    [("(0,1)(2)", "(0,2)(1)", 20),
     ("(0,3,2,1)", "(0,3)(1,2)", 20),
     ("(0,1,2,3,4)", "(0)(1,4)(2,3)", 20)])
def test_rotate(pr, pl, repeat):
    pytest.importorskip('surface_dynamics')
    for _ in range(repeat):
        fp = random_flip_sequence(pr, pl, loop=True)
        n = fp._start._n
        fp2 = fp.copy()
        fp2.rotate(random.randint(0, fp.num_blocks()))
        fp2._check()
        assert fp2.is_loop()
        assert fp.matrix().charpoly() == fp2.matrix().charpoly()

        fp2.rotate(random.randint(0, fp.num_blocks()))
        fp2._check()
        assert fp2.is_loop()
        assert fp.matrix().charpoly() == fp2.matrix().charpoly()

@pytest.mark.parametrize("pr, pl, repeat",
    [("(0,1)(2)", "(0,2)(1)", 20),
     ("(0,3,2,1)", "(0,3)(1,2)", 20),
     ("(0,1,2,3,4)", "(0)(1,4)(2,3)", 20)])
def test_conjugate(pr, pl, repeat):
    pytest.importorskip('surface_dynamics')
    for _ in range(repeat):
        fp = random_flip_sequence(pr, pl, loop=True)
        n = fp._start._n

        fp2 = fp.copy()
        fp2._check()
        assert fp2.is_loop()
        fp2.rotate(random.randint(0, fp.num_blocks()))
        fp2._check()
        assert fp2.is_loop()
        fp2.relabel(perm_random(n))
        fp2._check()
        assert fp2.is_loop()

        assert fp.matrix().charpoly() == fp2.matrix().charpoly()
        assert fp.is_conjugate(fp2)


def test_failure_conjugate_in_H6():
    data =  [("(0,6,5)(1,4,3,2)", "(0,6)(1,3,5)(2)(4)",
              [((1,), (1,)), ((1,), (1,)), ((0,), (1,)), ((0,), (1,))],
              "(0,4,1,5,2,6,3)"),
             ("(0,6,5)(1,4,3,2)", "(0,6)(1,3,5)(2)(4)",
              [((1,), (1,)), ((1,), (1,)), ((0,), (1,)), ((5,), (1,))],
              "(0,4,1,5,2,6,3)"),
             ("(0,6,5)(1,4,3,2)", "(0,6)(1,2,5)(3,4)",
              [((1,), (1,)), ((4,), (1,)), ((1,), (1,)), ((1,), (1,))],
              "(0,4,1,5,2,6,3)"),
             ("(0,4,2)(1)(3)(5,6)", "(0,5,6)(1,2,3,4)",
              [((1,), (1,)), ((1,), (1,)), ((0,), (1,)), ((0,), (1,))],
              "(0,3,6,2,5,1,4)"),
             ("(0,6,5)(1,4,3,2)", "(0,6)(1,2,5)(3,4)",
              [((1,), (1,)), ((2,), (1,)), ((1,), (1,)), ((1,), (1,))],
              "(0,4,1,5,2,6,3)"),
             ("(0,4,2)(1)(3)(5,6)", "(0,5,6)(1,2,3,4)",
              [((3,), (1,)), ((1,), (1,)), ((0,), (1,)), ((0,), (1,))],
              "(0,3,6,2,5,1,4)"),
             ("(0,4,1)(2,3)(5,6)", "(0,5,6)(1,2,3,4)",
              [((0,), (1,)), ((0,), (1,)), ((0,), (1,)), ((0,), (1,))],
              "(0,3,6,2,5,1,4)"),
             ("(0,4,1)(2,3)(5,6)", "(0,5,6)(1,2,3,4)",
              [((0,), (1,)), ((0,), (1,)), ((5,), (1,)), ((0,), (1,))],
              "(0,3,6,2,5,1,4)")]

    pytest.importorskip('surface_dynamics')

    flip_sequences = [VeeringQuadrangulationFlipSequence(VeeringQuadrangulation(pr, pl), flips, relabelling) for pr, pl, flips, relabelling in data]

    N = len(flip_sequences)

    m = [[0] * N for _ in range(N)]
    for i, fp1 in enumerate(flip_sequences):
        for j, fp2 in enumerate(flip_sequences):
            m[i][j] = fp1.is_conjugate(fp2)

    for i in range(N):
        assert m[i][i], m
        for j in range(i):
            assert m[i][j] == m[j][i]
        assert sum(m[i][j] for j in range(N)) == 4, m
        assert sum(m[j][i] for j in range(N)) == 4, m

#def test_failure_conjugate_in_H1_1():
#    data = [("(0,3,2,1)", "(0,2)(1)(3)", [((0,), (2,)), ((0, 3, 2), (1, 0, 0))], "(0,1,2,3)"),
#            ("(0,3,2,1)", "(0,2)(1)(3)", [((0,), (2,)), ((0, 3, 2), (0, 0, 1))], "(0,1,2,3)"),
#            ("(0,2)(1)(3)", "(0,1,2,3)", [((0, 1, 3), (0, 1, 0)), ((1,), (2,))], "(0,3,2,1)"),
#            ("(0,2)(1)(3)", "(0,1,2,3)", [((0, 1, 3), (0, 0, 1)), ((3,), (2,))], "(0,3,2,1)")]
#
#    flip_sequences = [VeeringQuadrangulationFlipSequence(VeeringQuadrangulation(pr, pl), flips, relabelling) for pr, pl, flips, relabelling in data]
#
#    N = len(flip_sequences)
#
#    m = [[0] * N for _ in range(N)]
#    for i, fp1 in enumerate(flip_sequences):
#        for j, fp2 in enumerate(flip_sequences):
#            m[i][j] = fp1.is_conjugate(fp2)
#
#    for i in range(N):
#        assert m[i][i], m
#        for j in range(i):
#            assert m[i][j] == m[j][i]
#        assert sum(m[i][j] for j in range(N)) == 2, m
#        assert sum(m[j][i] for j in range(N)) == 2, m
#
#def test_multiplicities(k, length):
#    # check that each loop is seen with the correct multiplicities
#    flip_lengths = {}
#    for fp in all_reduced_and_complete_loops(k, length, verbose=False):
#        cp = fp.matrix().charpoly()
#        if cp not in flip_lengths:
#            flip_lengths[cp] = []
#        flip_lengths[cp].append(len(fp._flips))
#
#    # multiplicities must be the sum of lengths
#    for cp in flip_lengths:
#        l = flip_lengths[cp]
#        l.sort()
#        for i in set(l):
#            assert any(l.count(i) % i == 0, (cp, flip_lengths[cp], i)
#
#    # inverses must appear with the same multiplicity
#    for cp in flip_lengths:
#        if cp[0] == -1:
#            reciprocal = -cp.parent()(list(cp)[::-1])(-cp.parent().gen())
#        else:
#            reciprocal = cp.parent()(list(cp)[::-1])
#        assert flip_lengths[cp] == flip_lengths[reciprocal]

if __name__ == '__main__':
    import sys
    sys.exit(pytest.main(sys.argv))
