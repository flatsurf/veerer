r"""
Geometric and balanced veering triangulations in small strata

For each connected component of stratum we compute the number of geometric/balanced triangulations.

genus 2 (holomorphic)
---------------------

H(2)
not geom. not bal. (32)
example: VeeringTriangulation("(0,~8,~3)(1,6,~2)(2,~1,~0)(3,7,~4)(4,8,~5)(5,~7,~6)", "RBBBRBRBB")

geom. not bal. (2)
example: VeeringTriangulation("(0,~2,1)(2,~6,~3)(3,5,~4)(4,8,~5)(6,~8,~7)(7,~1,~0)", "RBRBRBBRB")

geom. bal. (52)
example: VeeringTriangulation("(0,~6,~3)(1,7,~2)(2,~1,~0)(3,5,~4)(4,8,~5)(6,~8,~7)", "RBBBRBBRB")

H(1,1)
not geom. not bal. (480)
example: VeeringTriangulation("(0,~3,2)(1,5,~2)(3,9,~4)(4,~11,~5)(6,10,~7)(7,11,~8)(8,~10,~9)(~6,~1,~0)", "RBBBBRBRBRBB")

geom. not bal. (22)
example: VeeringTriangulation("(0,~2,1)(2,9,~3)(3,~7,~4)(4,6,~5)(5,11,~6)(7,~11,~8)(8,~10,~9)(10,~1,~0)", "RBRBRBBBBBRR")

geom. bal. (374)
example: VeeringTriangulation("(0,9,~3)(1,10,~2)(2,~1,~0)(3,~7,~4)(4,6,~5)(5,11,~6)(7,~11,~8)(8,~10,~9)", "RBBBRBBBBBRR")


Q(2,2)
not geom. not bal. (598)
example: VeeringTriangulation("(0,9,~8)(1,10,~2)(2,11,~3)(3,~10,~4)(4,~9,~5)(5,8,~6)(6,~11,~7)(7,~1,~0)", "RBBBBRBBBBRR")

geom. not bal. (48)
example: VeeringTriangulation("(0,~9,8)(1,10,~2)(2,~6,~3)(3,5,~4)(4,9,~5)(6,11,~7)(7,~10,~8)(~11,~1,~0)", "RBBBBRRBBBRB")

geom. bal. (200)
example: VeeringTriangulation("(0,7,~3)(1,5,~2)(2,~1,~0)(3,~11,~4)(4,~6,~5)(6,~8,~7)(8,10,~9)(9,11,~10)", "RBBBRRBRRBBR")


Q(2,1,1)
not geom. not bal. (13432)
example: VeeringTriangulation("(0,~11,10)(1,~13,~2)(2,12,~3)(3,~14,~4)(4,13,~5)(5,~9,~6)(6,8,~7)(7,11,~8)(9,~12,~10)(14,~1,~0)", "RBBBBBRBBBBRRRR")

geom. not bal. (764)
example: VeeringTriangulation("(0,~11,10)(1,~13,~2)(2,12,~3)(3,~7,~4)(4,6,~5)(5,11,~6)(7,~14,~8)(8,13,~9)(9,~12,~10)(14,~1,~0)", "RBBBBBRRBBBBRRB")

geom. bal. (2700)
example: VeeringTriangulation("(0,9,~8)(1,13,~2)(2,~14,~3)(3,~13,~4)(4,12,~5)(5,~11,~6)(6,10,~7)(7,~1,~0)(8,~10,~9)(11,14,~12)", "RBBBBRBBBRRBBRR")


Q(1,1,1,1)
not geom. not bal. (92736)
example: VeeringTriangulation("(0,~10,9)(1,17,~2)(2,~16,~3)(3,15,~4)(4,~17,~5)(5,16,~6)(6,~13,~7)(7,~12,~8)(8,11,~9)(10,12,~11)(13,~15,~14)(14,~1,~0)", "RBBBBBBRBBRRBBBRRR")

geom. not bal. (2926)
example: VeeringTriangulation("(0,9,14)(1,~16,~2)(2,15,~3)(3,~11,~4)(4,~8,~5)(5,7,~6)(6,10,~7)(8,~10,~9)(11,~17,~12)(12,16,~13)(13,~15,~14)(17,~1,~0)", "RBBBBBBRRBBRBBBRRB")

geom. bal. (11410)
example: VeeringTriangulation("(0,10,16)(1,11,~2)(2,~17,~3)(3,~14,~4)(4,13,~5)(5,~9,~6)(6,8,~7)(7,15,~8)(9,~11,~10)(12,17,~13)(14,~16,~15)(~12,~1,~0)", "RBBBBRBBRBRRBBRBBR")


genus 3 (holomorphic)
---------------------

H(4)^hyp
not geom. not bal. (6200)
example: VeeringTriangulation("(0,~4,3)(1,~7,~2)(2,6,~3)(4,~13,~5)(5,12,~6)(7,~12,~8)(8,11,~9)(9,~14,~10)(10,13,~11)(14,~1,~0)", "RBBBBBRRBRBBBRB")

geom. not bal. (220)
example: VeeringTriangulation("(0,~11,~3)(1,12,~2)(2,~1,~0)(3,10,~4)(4,~8,~5)(5,7,~6)(6,14,~7)(8,~14,~9)(9,13,~10)(11,~13,~12)", "RBBBRBRBBRRBRBB")

geom. bal. (2696)
example: VeeringTriangulation("(0,~14,~3)(1,12,~2)(2,~1,~0)(3,~7,~4)(4,6,~5)(5,13,~6)(7,~13,~8)(8,~12,~9)(9,11,~10)(10,14,~11)", "RBBBRBBBBBRBRRB")


H(4)^odd
not geom. not bal. (12036)
example: VeeringTriangulation("(0,~3,2)(1,~13,~2)(3,~14,~4)(4,~11,~5)(5,10,~6)(6,~12,~7)(7,11,~8)(8,13,~9)(9,14,~10)(12,~1,~0)", "RBBBBBRBBBBRBRR")

not geom. bal. (40)
example: VeeringTriangulation("(0,~13,12)(1,~11,~2)(2,10,~3)(3,~12,~4)(4,~9,~5)(5,8,~6)(6,~14,~7)(7,13,~8)(9,11,~10)(14,~1,~0)", "RBBRBBBBRRBRBBR")

geom. not bal. (12652)
example: VeeringTriangulation("(0,~11,~5)(1,12,~2)(2,~9,~3)(3,8,~4)(4,~1,~0)(5,10,~6)(6,14,~7)(7,~13,~8)(9,13,~10)(11,~14,~12)", "RBBBBBRBRRBBRBB")

geom. bal. (22824)
example: VeeringTriangulation("(0,~14,~5)(1,12,~2)(2,~9,~3)(3,8,~4)(4,~1,~0)(5,~12,~6)(6,11,~7)(7,~13,~8)(9,13,~10)(10,14,~11)", "RBBBBBRBRRBBRBR")




MORE TO COME
------------

Q(8)

H(3,1)

H(2,2)^hyp

H(2,2)^odd

Q(7, 1)

Q(6, 2)^hyp

Q(6, 2)^nonhyp

Q(5,3)

Q(4,4)

"""


from __future__ import absolute_import, print_function

from collections import defaultdict

import veerer
import surface_dynamics

def display_stats(stratum):
    A = veerer.Automaton.from_stratum(stratum)
    counts = defaultdict(int)
    examples = {}
    for t in A:
        geom = t.is_geometric()
        bal = t.is_balanced()
        stat = (geom,bal)

        counts[stat] += 1

        if stat not in examples:
            examples[stat] = t

    for stat in sorted(examples):
        geom, bal = stat
        print("%s %s (%d)" % ('geom.' if geom else 'not geom.',
                         'bal.' if bal else 'not bal.',
                          counts[stat]))
        print("example: %s" % examples[stat])
        print()
