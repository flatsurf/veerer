.. -*- coding: utf-8 -*-
.. linkall

Ferenczi-Zamboni induction
==========================

Ferenczi-Zamboni induction is a particularly nice way to do veering
flips in hyperelliptic strata of Abelian differentials (or equivalently,
minimal strata of quadratic differentials on the sphere).

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

    sage: V = Vc.copy()

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
    (1.722083805739043?, Flat Triangulation made of 3 triangles)

    sage: f = R1 * R1 * R5
    sage: assert f.is_pseudo_anosov()
    sage: f.self_similar_surface()
    (2.296630262886538?, Flat Triangulation made of 3 triangles)

    sage: f = R3 * R1 * R2 * CL5
    sage: assert f.is_pseudo_anosov()
    sage: f.self_similar_surface()
    (3.732050807568878?, Flat Triangulation made of 3 triangles)

    sage: f = R3 * R1 * R2 * CL5 * CR5
    sage: assert f.is_pseudo_anosov()
    sage: f.self_similar_surface()
    (4.309326438416744?, Flat Triangulation made of 3 triangles)
