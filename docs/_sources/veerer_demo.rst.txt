.. -*- coding: utf-8 -*-
.. linkall

Veering triangulation demo
==========================


:Authors:
    - Vincent Delecroix
    - Saul Schleimer
:License: CC BY-NC-SA 3.0

-  27 Aug. 2018 Toronto
-  20 Sept. 2018 Temple
-  25 Sept. 2018 Villetaneuse and CUNY

``veerer`` is a Python library for exploration of veering
triangulations. It is written by
`Mark Bell <https://markcbell.github.io>`_,
`Vincent Delecroix <http://www.labri.fr/perso/vdelecro/>`_ and
`Saul Schleimer <http://homepages.warwick.ac.uk/~masgar/>`_. It is
part of a project that also involve
`Vaibhav Gadre <http://www.maths.gla.ac.uk/~vgadre/>`_ and
`Rodolfo Gutiérrez-Romo <http://rodol.fo>`_.

``veerer`` works in conjunction with ``pplpy`` (for rational polytope
computations and optimiation) and ``surface_dynamics`` (for analyzing
stratum components). The plotting part is only available in SageMath.

::

   sage: from veerer import *
   sage: from surface_dynamics import *   # optional - surface_dynamics

To input a triangulation in the program one needs to specify a list of
triangles ``(i, j, k)`` and the pairing of edges is by convention given
by ``i <-> ~i``. The veering coloring is then specified by a list of
colors.

::

    sage: # standard torus made of two triangles
    sage: faces0 = '(0,1,2)(~0,~1,~2)'
    sage: colours0 = 'RRB'
    sage: T0 = VeeringTriangulation(faces0, colours0)

::

    sage: T0.genus()
    1

::

    sage: T0.angles()
    [2]
    sage: T0.stratum()  # optional - surface_dynamics
    H_1(0)

::

    sage: # triangulation is core if it carries some flat structure
    sage: T0.is_core()
    True

::

    sage: FS0 = T0.flat_structure_min()
    sage: FS0.plot().show(figsize=5)

::

    sage: FFS0 = T0.flat_structure_geometric_middle()
    sage: FFS0.plot().show(figsize=5)

::

    sage: # an Abelian genus 2 example
    sage: faces1 = '(0, ~3, 4)(1, 2, ~7)(3, ~1, ~2)(5, ~8, ~4)(6, ~5, 8)(7, ~6, ~0)'
    sage: colours1 = 'RBRRBRBRB'
    sage: T1 = VeeringTriangulation(faces1, colours1)

::

    sage: T1.angles()
    [6]
    sage: T1.stratum()  # optional - surface_dynamics
    H_2(2)

::

    sage: FS1 = T1.flat_structure_min()
    sage: FS1.plot().show(figsize=5)

::

    sage: FFS1 = T1.flat_structure_geometric_middle()
    sage: FFS1.plot().show(figsize=5)

::

    sage: # a quadratic genus 1 example
    sage: faces2 = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
    sage: colours2 = 'BRBRBBBRBR'
    sage: T2 = VeeringTriangulation(faces2, colours2)

::

    sage: T2.angles()
    [3, 3, 1, 1]
    sage: T2.stratum()  # optional - surface_dynamics
    Q_1(1^2, -1^2)

::

    sage: FS2 = T2.flat_structure_min()
    sage: FS2.plot().show(figsize=5)  # not tested (warning from matplotlib)


Viewing train-tracks!
---------------------

Recall that a veering triangulation is just a pair of transversal
train-tracks.

::

    sage: TT_horiz = FS1.plot(horizontal_train_track=True, edge_labels=False)
    sage: TT_vert = FS1.plot(vertical_train_track=True, edge_labels=False)
    sage: graphics_array([TT_horiz, TT_vert], 1, 2).show(figsize=6)

::

    sage: # constructing veering triangulations from a component stratum
    sage: Q = QuadraticStratum(3,3,3,3)     # optional - surface_dynamics
    sage: Qreg = Q.regular_component()      # optional - surface_dynamics
    sage: Qirr = Q.irregular_component()    # optional - surface_dynamics

::

    sage: VeeringTriangulation.from_stratum(Qreg)    # optional - surface_dynamics
    VeeringTriangulation("(0,19,~18)(1,20,~19)(2,21,~20)(3,22,~21)(4,23,~22)(5,25,~24)(6,27,~26)(7,28,~27)(8,29,~28)(9,~16,17)(10,~5,~29)(11,~6,~10)(12,~1,~11)(13,~9,~12)(14,~7,~13)(15,~2,~14)(16,~0,~15)(18,~8,~17)(24,~23,~3)(26,~25,~4)", "RRRRRRRRRRBBBBBBBBBBBBBBBBBBBB")

::

    sage: VeeringTriangulation.from_stratum(Qirr)   # optional - surface_dynamics
    VeeringTriangulation("(0,21,~20)(1,22,~21)(2,23,~22)(3,24,~23)(4,25,~24)(5,26,~25)(6,28,~27)(7,29,~28)(8,~16,17)(9,~14,15)(10,~2,~29)(11,~1,~10)(12,~6,~11)(13,~9,~12)(14,~8,~13)(16,~4,~15)(18,~5,~17)(19,~3,~18)(20,~7,~19)(27,~26,~0)", "RRRRRRRRRRBBBBBBBBBBBBBBBBBBBB")

::

    sage: # constructing a veering triangulation from a pseudo-Anosov homeomorphism
    sage: import flipper                    # optional - flipper
    sage: S_2_1 = flipper.load('S_2_1')     # optional - flipper
    sage: h = S_2_1.mapping_class('abcD')   # optional - flipper
    sage: print(h.nielsen_thurston_type())  # optional - flipper
    Pseudo-Anosov

::

    sage: VeeringTriangulation.from_pseudo_anosov(h)  # optional - flipper
    VeeringTriangulation("(0,~3,~1)(1,2,14)(3,~5,~13)(4,~12,~8)(5,6,~11)(7,8,13)(9,~6,~7)(10,~0,11)(12,~14,~10)(~9,~4,~2)", "RBRBRRBRBBBBRBR")

Core vs not core
----------------

::

    sage: # start from our surface in H(2) and let us flip some edges
    sage: S = T1.copy()
    sage: print(S.is_core())
    True
    sage: print(S.flippable_edges())
    [0, 2, 3, 7, 8]

::

    sage: S.flip(3, BLUE)
    sage: print(S.is_core())
    True
    sage: print(S.flippable_edges())
    [3, 7, 8]

::

    sage: S.flip(8, BLUE)
    sage: print(S.is_core())
    True
    sage: print(S.flippable_edges())
    [3, 4, 7, 8]

::

    sage: S.flip(4, RED)
    sage: print(S.is_core())
    True
    sage: print(S.flippable_edges())
    [4, 7]

::

    sage: FS = S.flat_structure_min()
    sage: FS.plot()
    Graphics object consisting of 37 graphics primitives

::

    sage: # in the geometric setting, the flipped edge is forced to be BLUE
    sage: S.flip(7, RED)
    sage: S.is_core()
    False

::

    sage: print(S.train_track_polytope(HORIZONTAL))
    A 4-dimensional polyhedron in QQ^9 defined as the convex hull of 1 point, 5 rays
    sage: print(S.train_track_polytope(VERTICAL))
    A 3-dimensional polyhedron in QQ^9 defined as the convex hull of 1 point, 3 rays

::

    sage: # check that we indeed started with a core veering triangulation
    sage: print(T1.train_track_polytope(HORIZONTAL))
    A 4-dimensional polyhedron in QQ^9 defined as the convex hull of 1 point, 4 rays
    sage: print(T1.train_track_polytope(VERTICAL))
    A 4-dimensional polyhedron in QQ^9 defined as the convex hull of 1 point, 5 rays


Geometric polytope
------------------


A triangulation is *geometric* if it is the L^infinity-Delaunay triangulation of
some flat structure

::

    sage: # triangulation of some flat structure
    sage: T0.is_geometric()
    True

The geometric polytope that parametrizes the geometric vectors is a sub-polytope
of the product of the two train-track polytopes.

::

    sage: print(T1.is_geometric())
    True
    sage: print(T1.geometric_polytope())
    A 8-dimensional polyhedron in QQ^18 defined as the convex hull of 1 point, 61 rays

Core automaton
--------------

The core automaton of a given triangulations `T_0` is the directed graph whose
vertices are core veering triangulations that can be reached from `T_0` by a
sequence of flips and there is a directed edge `T_i \to T_j` if `T_j` is obtained
from `T_i` by a flip.

::

    sage: # T0 was the torus example
    sage: from veerer import CoreAutomaton
    sage: A0 = CoreAutomaton(T0)
    sage: A0
    Core veering automaton with 2 vertices

::

    sage: print(A0.num_states(), A0.num_transitions())
    2 4
    sage: print(A0.num_geometric_triangulations())
    2
    sage: print(A0.num_cylindrical_triangulations())
    2

::

    sage: # T1 was the genus 2 example in H(2)
    sage: A1 = CoreAutomaton(T1)

::

    sage: print(A1.num_states(), A1.num_transitions())
    86 300
    sage: print(A1.num_geometric_triangulations())
    54
    sage: print(A1.num_cylindrical_triangulations())
    24

::

    sage: # T2 was the genus 1 example in Q(1^2, -1^2)
    sage: A2 = CoreAutomaton(T2)
    sage: print(A2.num_states(), A2.num_transitions())
    1074 3620
    sage: print(A2.num_geometric_triangulations())
    270
    sage: print(A2.num_cylindrical_triangulations())
    196

Some data (orientable case)
---------------------------

+---------------------+-----+---------+-----------+-------------+
| component           | dim | core    | geometric | cylindrical |
+=====================+=====+=========+===========+=============+
| H(0)                | 2   | 2       | 2         | 2           |
+---------------------+-----+---------+-----------+-------------+
| H(2)                | 4   | 86      | 54        | 24          |
+---------------------+-----+---------+-----------+-------------+
| H(1,1)              | 5   | 876     | 396       | 136         |
+---------------------+-----+---------+-----------+-------------+
| H(4)^hyp            | 6   | 9116    | 2916      | 636         |
+---------------------+-----+---------+-----------+-------------+
| H(4)^odd            | 6   | 47552   | 35476     | 1970        |
+---------------------+-----+---------+-----------+-------------+
| H(2,2)^hyp          | 7   | 111732  | 24192     | 3934        |
+---------------------+-----+---------+-----------+-------------+
| H(2,2)^odd          | 7   | 874750  | 711568    | 12740       |
+---------------------+-----+---------+-----------+-------------+
| H(3,1)              | 7   | 2011366 | 1317136   | 33164       |
+---------------------+-----+---------+-----------+-------------+

To give an idea about the complexity and timings when generating the
above data, here are the steps involved. The timings are for the stratum
component H(4)^hyp that is the fourth row in the above array: -
generating the core graph ~20 secs for H(4)^hyp (the graph has 9116
vertices and 44664 edges) - filtering the geometric triangulations
(single test involves a polytope computation) ~20 secs for H(4)^hyp -
filtering cylindrical (single test is cheap) ~2 sec for H(4)^hyp

::

    sage: H = AbelianStratum(4).hyperelliptic_component()  # optional - surface_dynamics
    sage: V = VeeringTriangulation.from_stratum(H)         # optional - surface_dynamics
    sage: AV = CoreAutomaton(V)                          # long time - ~21 secs # optional - surface_dynamics
    sage: print(AV.num_states())                         # long time - ~150 µs # optional - surface_dynamics
    9116
    sage: print(AV.num_geometric_triangulations())       # long time - ~21 secs # optional - surface_dynamics
    2916
    sage: print(AV.num_cylindrical_triangulations())     # long time - ~1.5 secs # optional - surface_dynamics
    636

License
-------

This document is published under the Creative Commons
`CC BY-SA 3.0 <https://creativecommons.org/licenses/by-sa/3.0/>`_.
