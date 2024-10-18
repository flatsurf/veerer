.. -*- coding: utf-8 -*-
.. linkall

Ferenczi-Zamboni induction
==========================

Ferenczi-Zamboni induction is a particularly nice way to do veering
flips in hyperelliptic strata of Abelian differentials (or equivalently,
minimal strata of quadratic differentials on the sphere). They
are defined in [FeZa10]_, [CaFeZa11]_ and put into a more geometrical
context in [DeUl15]_.

Torus case or `\mathcal{Q}(-1^4)`
---------------------------------

There are only two possible moves::

    sage: from veerer import VeeringTriangulation, VeeringFlipSequence # random output due to deprecation warnings from realalg
    sage: V = VeeringTriangulation("(0,1,2)", "PBR")
    sage: R = VeeringFlipSequence(V, "0R", [2,1,0])
    sage: L = VeeringFlipSequence(V, "0B", [1,0,2])

The smallest dilatation is the golden rotation::

    sage: assert R.is_closed() and L.is_closed()
    sage: fp = R * L
    sage: fp
    VeeringFlipSequence(VeeringTriangulation("(0,1,2)", "PBR"), "0R 2B", "(0,2,1)")
    sage: a, S = fp.self_similar_surface()
    sage: SS = S.copy(mutable=True)
    sage: SS.flip(0)
    sage: SS.flip(2)
    sage: SS.relabel("(0,2,1)")
    sage: SS.xy_scaling(-a, -1/a)
    sage: SS == S
    True

`\mathcal{H}(2)` or `\mathcal{Q}(1,-1^5)`
-----------------------------------------

The case of 3 triangles correspond in terms of strata to H(2) (which is
the canonical orientation cover of Q(1,-1^5) on the sphere). The
Ferenczi-Zamboni automaton is made of three triangulations. We denote them
below by ``Vc``, ``Vl`` and ``Vr`` where c, l and r respectively stand for
center, left and right.

::

    sage: from veerer import VeeringTriangulation, BLUE, RED, PURPLE

    sage: Vc = VeeringTriangulation("(0,~5,4)(3,5,6)(1,2,~6)", "PPBPRBR")
    sage: Vr = VeeringTriangulation("(0,6,5)(1,2,~6)(3,4,~5)", "BPBBRPR")
    sage: Vl = VeeringTriangulation("(0,~5,4)(1,6,5)(2,3,~6)", "PRBRRBP")

    sage: V = Vc.copy(mutable=True)

    sage: V.flip(1, BLUE)
    sage: V.relabel("(1,2)")
    sage: assert V == Vc

    sage: V.flip(0, RED)
    sage: V.relabel("(0,4)")
    sage: assert V == Vc

    sage: V.flip(0, BLUE)
    sage: V.flip(3, BLUE)
    sage: V.relabel("(0,3)")
    sage: assert V == Vr

    sage: V.flip(1, BLUE)
    sage: V.relabel("(1,2)")
    sage: assert V == Vr

    sage: V.flip(1, RED)
    sage: V.flip(5, RED)
    sage: V.relabel("(0,2,3)(1,4)(5,6)")
    sage: assert V == Vr

    sage: V.flip(5, BLUE)
    sage: assert V == Vc

    sage: V.flip(1, RED)
    sage: V.flip(3, RED)
    sage: V.relabel("(1,3)")
    sage: assert V == Vl

    sage: V.flip(0, RED)
    sage: V.relabel("(0,4)")
    sage: assert V == Vl

    sage: V.flip(0, BLUE)
    sage: V.flip(6, BLUE)
    sage: V.relabel("(0,2)(1,4,3)(5,6,~5,~6)")
    sage: assert V == Vl

    sage: V.flip(6, RED)
    sage: V.relabel("(6,~6)")
    sage: assert V == Vc

Now, instead of modifying a given triangulation we instead define flip sequences

::

    sage: from veerer import VeeringFlipSequence

    sage: CR5 = VeeringFlipSequence(Vc, "1B", "(1,2)")
    sage: CL5 = VeeringFlipSequence(Vc, "0R", "(0,4)")
    sage: R3 = VeeringFlipSequence(Vc, "0B 3B", "(0,3)")
    sage: R5 = VeeringFlipSequence(Vr, "1B", "(1,2)")
    sage: R1 = VeeringFlipSequence(Vr, "1R 5R", "(0,2,3)(1,4)(5,6)")
    sage: R2 = VeeringFlipSequence(Vr, "5B")
    sage: L3 = VeeringFlipSequence(Vc, "1R 3R", "(1,3)")
    sage: L5 = VeeringFlipSequence(Vl, "0R", "(0,4)")
    sage: L1 = VeeringFlipSequence(Vl, "0B 6B", "(0,2)(1,4,3)(5,6,~5,~6)")
    sage: L2 = VeeringFlipSequence(Vl, "6R", "(6,~6)")

    sage: assert CL5.start() == CL5.end() == Vc
    sage: assert CR5.start() == CR5.end() == Vc
    sage: assert R3.start() == Vc and R3.end() == Vr
    sage: assert R5.start() == R5.end() == Vr
    sage: assert R1.start() == R1.end() == Vr
    sage: assert R2.start() == Vr and R2.end() == Vc
    sage: assert L3.start() == Vc and L3.end() == Vl
    sage: assert L5.start() == L5.end() == Vl
    sage: assert L1.start() == L1.end() == Vl
    sage: assert L2.start() == Vl and L2.end() == Vc

They can be composed and one can check whether they define pseudo-Anosov homeomorphism::

    sage: (R3 * R2 * CR5).is_pseudo_anosov()
    False
    sage: (R3 * R5 * R2 * L3 * L5 * L2).is_pseudo_anosov()
    True

Some pseudo-Anosov with small dilatation in H(2)

::

    sage: f = R1 * R5
    sage: assert f.is_pseudo_anosov()
    sage: f.self_similar_surface()
    (a,
     FlatVeeringTriangulation(Triangulation("(0,6,5)(1,2,~6)(3,4,~5)"), [(1, -1), (a, a^3 - a^2 - a - 1), (a^3 - 2*a - 2, a^2), (-a^3 + a^2 + a + 1, -a), (2*a^3 - a^2 - 2*a - 2, a^3 - 2), (-a^3 + a + 1, -a^3 + a + 2), (a^3 - a - 2, a^3 - a - 1), (-a^3 + a + 2, -a^3 + a + 1), (-a^3 + a + 1, -a^3 + a + 2)]))

    sage: f = R1 * R1 * R5
    sage: assert f.is_pseudo_anosov()
    sage: f.self_similar_surface()
    (a,
     FlatVeeringTriangulation(Triangulation("(0,6,5)(1,2,~6)(3,4,~5)"), [(1, -1), (a^2, 2*a^3 - 3*a^2 - 2*a - 4), (a^3 - 2*a^2 - 2, a), (a, a^3 - 2*a^2 - 2), (a^3 - a^2 - a - 1, a^3 - a^2 - a - 3), (-a^3 + a^2 + 1, -2*a^3 + 3*a^2 + a + 5), (a^3 - a^2 - 2, 2*a^3 - 3*a^2 - a - 4), (-a^3 + a^2 + 2, -2*a^3 + 3*a^2 + a + 4), (-a^3 + a^2 + 1, -2*a^3 + 3*a^2 + a + 5)]))

    sage: f = R3 * R1 * R2 * CL5
    sage: assert f.is_pseudo_anosov()
    sage: f.self_similar_surface()
    (a,
     FlatVeeringTriangulation(Triangulation("(0,~5,4)(1,2,~6)(3,5,6)"), [(1, 1), (1, 1), (-1/2*a + 3/2, 1/2*a - 1/2), (1/2*a - 1/2, -1/2*a + 3/2), (a - 4, -a), (-a + 3, a - 1), (1/2*a - 5/2, -1/2*a - 1/2), (1/2*a - 5/2, -1/2*a - 1/2), (-a + 3, a - 1)]))

    sage: f = R3 * R1 * R2 * CL5 * CR5
    sage: assert f.is_pseudo_anosov()
    sage: f.self_similar_surface()
    (a,
     FlatVeeringTriangulation(Triangulation("(0,~5,4)(1,2,~6)(3,5,6)"), [(1, 1), (7/33*a^3 - 23/33*a^2 - 19/33*a - 25/33, -10/33*a^3 + 32/33*a^2 + 37/33*a + 16/33), (-20/33*a^3 + 61/33*a^2 + 92/33*a + 62/33, 5/33*a^3 - 16/33*a^2 - 2/33*a - 8/33), (-1/33*a^3 + 8/33*a^2 - 2/33*a - 20/33, -8/33*a^3 + 19/33*a^2 + 56/33*a + 26/33), (4/11*a^3 - 10/11*a^2 - 25/11*a - 30/11, -1/11*a^3 + 1/11*a^2 + 7/11*a - 5/11), (-4/11*a^3 + 10/11*a^2 + 25/11*a + 19/11, 1/11*a^3 - 1/11*a^2 - 7/11*a - 6/11), (13/33*a^3 - 38/33*a^2 - 73/33*a - 37/33, 5/33*a^3 - 16/33*a^2 - 35/33*a - 8/33), (13/33*a^3 - 38/33*a^2 - 73/33*a - 37/33, 5/33*a^3 - 16/33*a^2 - 35/33*a - 8/33), (-4/11*a^3 + 10/11*a^2 + 25/11*a + 19/11, 1/11*a^3 - 1/11*a^2 - 7/11*a - 6/11)]))
