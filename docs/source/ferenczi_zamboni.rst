.. -*- coding: utf-8 -*-
.. linkall

Ferenczi-Zamboni induction
==========================

Ferenczi-Zamboni induction is a particularly nice way to do veering
flips in hyperelliptic strata of Abelian differentials (or equivalently,
minimal strata of quadratic differentials on the sphere).

The case of H(2) corresponds to Q(1,-1^5). The automaton is made of three
canonical triangulations. We denote them below by ``V0``, ``Vright`` and
``Vleft``.

::

    sage: from veerer import VeeringTriangulation, BLUE, RED

    sage: V0 = VeeringTriangulation("(0,3,4)(1,~3,5)(2,6,~4)", "PPPBRRB")
    sage: Vright = VeeringTriangulation("(0,4,3)(1,5,~3)(2,6,~4)", "BBPPRRB")
    sage: Vleft = VeeringTriangulation("(0,4,3)(1,~3,5)(2,~4,6)", "RPRBPRB")

    sage: assert V0.angles() == Vright.angles() == Vleft.angles() == [3, 1, 1, 1, 1, 1]

    sage: assert V0.flippable_edges() == [0, 1, 2]
    sage: assert Vright.flippable_edges() == [2, 3]
    sage: assert Vleft.flippable_edges() == [1, 4]

    sage: V = V0.copy()

    sage: V.flip(2, BLUE, reduced=True)
    sage: V.relabel("(2,6)")
    sage: assert V == V0

    sage: V.flip(1, RED, reduced=True)
    sage: V.relabel("(1,5)")
    sage: assert V == V0

    sage: V.flip(0, BLUE, reduced=True)
    sage: V.flip(1, BLUE, reduced=True)
    sage: assert V == Vright

    sage: V.flip(2, RED, reduced=True)
    sage: V.flip(3, RED, reduced=True)
    sage: V.relabel("(0,6,1)(2,5)(3,4)")
    sage: assert V == Vright

    sage: V.flip(2, BLUE, reduced=True)
    sage: V.relabel("(2,6)")
    sage: assert V == Vright

    sage: V.flip(3, BLUE, reduced=True)
    sage: V.relabel("(0,1)")
    sage: assert V == V0

    sage: V.flip(0, RED, reduced=True)
    sage: V.flip(2, RED, reduced=True)
    sage: assert V == Vleft

    sage: V.flip(1, RED, reduced=True)
    sage: V.relabel("(1,5)")
    sage: assert V == Vleft

    sage: V.flip(4, RED, reduced=True)
    sage: V.relabel("(0,2)(4,~4)")
    sage: assert V == V0

::

    sage: from veerer import VeeringFlipSequence

    sage: F0 = VeeringFlipSequence(V0, "2B", "(2,6)", reduced=True)
    sage: F1 = VeeringFlipSequence(V0, "1R", "(1,5)", reduced=True)
    sage: assert F0.is_closed() and F1.is_closed()
    sage: F2 = VeeringFlipSequence(V0, "0B 1B", reduced=True)
    sage: F3 = VeeringFlipSequence(Vright, "3B", "(0,1)", reduced=True)
    sage: assert not F2.is_closed() and not F3.is_closed()
    sage: assert F2._end == Vright and F3._end == V0
    sage: F4 = VeeringFlipSequence(Vright, "2R 3R", "(0,6,1)(2,5)(3,4)", reduced=True)
    sage: F5 = VeeringFlipSequence(Vright, "2B", "(2,6)", reduced=True)
    sage: assert F4.is_closed() and F5.is_closed()
    sage: F6 = VeeringFlipSequence(V0, "0R 2R", reduced=True)
    sage: F7 = VeeringFlipSequence(Vleft, "4R", "(0,2)(4,~4)", reduced=True)
    sage: assert not F6.is_closed() and not F7.is_closed()
    sage: assert F6._end == Vleft and F7._end == V0
    sage: F8 = None
    sage: F9 = VeeringFlipSequence(Vleft, "1R", "(1,5)", reduced=True)
    sage: assert F9.is_closed()

    sage: F0 * F2 * F4 * F3 * F1
    VeeringFlipSequence(VeeringTriangulation("(0,3,4)(1,~3,5)(2,6,~4)", "PPPBRRB"), "2B 0B 1B 6R 3R 4B 1R", (0,6,1,5,2)(3,4)(~4,~3), reduced=True)
