r"""
Flip automata.

We consider various automata related to edge flip in graphs on surfaces, namely
- triangulations and more generally flip automaton of combinatorial maps with fixed face degrees
- train-track splitting
- core veering triangulations of Abelian and quadratic differentials
- (L-infinity) Delaunay triangulations of Abelian and quadratic differentials
- Delaunay-Strebel for meromorphic

In order to avoid confusion with graph terminology (that might refer to the
underlying triangulation), the flip graphs are considered as automata.

EXAMPLES::

    sage: from veerer import *

    sage: vt = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
    sage: A = CoreAutomaton()
    sage: A.add_seed(vt)
    1
    sage: A.run()
    0
    sage: len(A)
    2

    sage: vt = VeeringTriangulation('(0,1,2)', 'BBR')
    sage: A = CoreAutomaton()
    sage: A.add_seed(vt)
    1
    sage: A.run()
    0
    sage: len(A)
    2

The examples below is Q(1^2, -1^2) with two folded edges::

    sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
    sage: cols = 'BRBRBBBRBR'
    sage: vt = VeeringTriangulation(fp, cols)

    sage: C = CoreAutomaton()
    sage: C.add_seed(vt)
    1
    sage: C.run()
    0
    sage: C
    Core veering automaton with 1074 vertices

    sage: R = ReducedCoreAutomaton()
    sage: R.add_seed(vt)
    1
    sage: R.run()
    0
    sage: R
    Reduced core veering automaton with 356 vertices

Exploring strata::

    sage: from surface_dynamics import Stratum                       # optional - surface_dynamics
    sage: strata = [Stratum([2], 1), Stratum([2,-1,-1], 2),          # optional - surface_dynamics
    ....:           Stratum([2,2], 2), Stratum([1,1], 1)]
    sage: for stratum in strata:                                     # optional - surface_dynamics
    ....:     print(stratum)
    ....:     vt = VeeringTriangulation.from_stratum(stratum)
    ....:     C = CoreAutomaton()
    ....:     _ = C.add_seed(vt)
    ....:     _ = C.run()
    ....:     R = ReducedCoreAutomaton()
    ....:     _ = R.add_seed(vt)
    ....:     _ = R.run()
    ....:     print(C)
    ....:     print(R)
    H_2(2)
    Core veering automaton with 86 vertices
    Reduced core veering automaton with 28 vertices
    Q_1(2, -1^2)
    Core veering automaton with 160 vertices
    Reduced core veering automaton with 68 vertices
    Q_2(2^2)
    Core veering automaton with 846 vertices
    Reduced core veering automaton with 305 vertices
    H_2(1^2)
    Core veering automaton with 876 vertices
    Reduced core veering automaton with 234 vertices
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2018 Mark Bell
#                     2018-2024 Vincent Delecroix
#                     2018 Saul Schleimer
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************

import collections
import os
import sys
from time import time

from .triangulation import Triangulation
from .veering_triangulation import VeeringTriangulation
from .strebel_graph import StrebelGraph
from .linear_family import VeeringTriangulationLinearFamily
from .constants import RED, BLUE, PURPLE, VERTICAL, HORIZONTAL, PROPERTIES_COLOURS, colour_to_char, colour_to_string
from .permutation import perm_invert

# TODO: when set to True a lot of intermediate checks are performed
CHECK = False


class Automaton:
    r"""
    Abstract class for automaton.

    INPUT:

    - ``backward`` -- whether to explore backward (only useful if the automaton
      is not strongly connected)

    - ``verbosity`` -- the level of verbosity (the higher the more messages
      you get while the automaton is being computed)

    - ``algorithm`` -- either ``'DFS'`` (depths first search) or ``'BFS'``
      (breadth first search)

    - ``extra_kwds`` -- additional argument that are forwarded to the method `_setup`
      that could be overridden in subclasses
    """
    _name = ''

    def __init__(self, backward=False, method='BFS', verbosity=0, **extra_kwds):
        # verbosity level
        self._verbosity = int(verbosity)

        # whether to explore backward flips
        if backward is not True and backward is not False:
            raise TypeError('the argument backward must be a boolean')
        self._backward = backward

        # the automaton is encoded in two dictionaries where the keys are the states
        # and the values are respectively the outgoing or incoming edges.
        # In both cases, the neighbors are encoded by pairs (neighbor, label)
        self._forward_neighbors = {}
        self._backward_neighbors = {}

        # list of seeds to be considered to explore the automaton from
        self._seeds = []

        # the current list of states whose neighbors are still to be explored
        self._branch = []

        # queue of states on which we still have run backward flips from
        # (each time an element is popped its backward flip neighbors will
        #  populate the _seeds list)
        self._backward_flip_queue = []

        if method not in ['BFS', 'DFS']:
            raise ValueError('method should either be \'BFS\' or \'DFS\'')
        self._method = method

        self._setup(**extra_kwds)

    def _check(self):
        r"""
        Some consistency checks.
        """
        if not (self._branch or self._backward_flip_queue or self._seeds):
            d1 = set(self._forward_neighbors).difference(self._backward_neighbors)
            assert not d1, d1

            d2 = set(self._backward_neighbors).difference(self._forward_neighbors)
            assert not d2, d2

            num_incoming_edges = sum(len(v) for v in self._backward_neighbors.values())
            num_outgoing_edges = sum(len(v) for v in self._forward_neighbors.values())
            assert num_outgoing_edges == num_incoming_edges, (num_outgoing_edges, num_incoming_edges)

            for state in self:
                in_neighbors1 = set(x for x, label in self._backward_neighbors[state])
                in_neighbors2 = set(x for x, label in self._in_neighbors(state))
                assert in_neighbors1 == in_neighbors2, (state, in_neighbors1, in_neighbors2)

                out_neighbors1 = set(x for x, label in self._forward_neighbors[state])
                out_neighbors2 = set(x for x, label in self._out_neighbors(state))
                assert out_neighbors1 == out_neighbors2, (state, out_neighbors1, out_neighbors2)

    def __str__(self):
        r"""
        Return a string representing this automaton
        """
        status = None
        if self._branch or self._backward_flip_queue or self._seeds:
            status = 'Partial '
            name = self._name
        else:
            status = ''
            name = self._name[0].upper() + self._name[1:]
        return ("%s%s automaton with %s %s" % (status, name, len(self._forward_neighbors), 'vertex' if len(self._forward_neighbors) <= 1 else 'vertices'))

    def __repr__(self):
        return str(self)

    def __len__(self):
        r"""
        Return the number of states in this automaton.

        EXAMPLES::

            sage: from veerer import *

            sage: vt = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton()
            sage: A.number_of_states()
            0
            sage: A.add_seed(vt)
            1
            sage: A.number_of_states()
            0
            sage: A.run()
            0
            sage: A.number_of_states()
            86
        """
        return len(self._forward_neighbors)

    number_of_states = __len__
    num_states = __len__

    def number_of_transitions(self):
        r"""
        Return the number of transitions of this automaton.

        EXAMPLES::

            sage: from veerer import *

            sage: vt = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton()
            sage: A.number_of_transitions()
            0
            sage: A.add_seed(vt)
            1
            sage: A.number_of_transitions()
            0
            sage: A.run()
            0
            sage: A.number_of_transitions()
            300
        """
        return sum(len(v) for v in self._forward_neighbors.values())

    num_transitions = number_of_transitions

    def __iter__(self):
        r"""
        Run through the states of this automaton.

        EXAMPLES::

            sage: from veerer import *

            sage: vt = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton()
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: sum(T.is_geometric() for T in A)
            54
            sage: sum(T.is_cylindrical() for T in A)
            24

            sage: vt = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton()
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: assert all(t.angles() == [6] for t in A)
        """
        return iter(self._forward_neighbors)

    states = __iter__

    def sources(self):
        r"""
        Iterate through sources (states with no backward neighbor) in this automaton.

        EXAMPLES::

            sage: from veerer import *
            sage: fp = "(0,2,1)(~0,3,~1)"
            sage: bdry = "(~2:2,~3:2)"
            sage: cols = "BRRR"
            sage: vt = VeeringTriangulation(fp, bdry, cols)

        In the Delaunay automaton, the sources are the horizontal-Strebel veering triangulations::

            sage: A = DelaunayAutomaton(backward=True)
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: list(A.sources())
            [VeeringTriangulation("(0,~3,2)(1,3,~2)", boundary="(~1:2,~0:2)", colouring="RRBR")]
            sage: set(A.sources()) == set(vt for vt in A if vt.is_strebel(HORIZONTAL))
            True

        In the Delaunay-Strebel automaton, the sources are the horizontal Strebel graphs::

            sage: A = DelaunayStrebelAutomaton(backward=True)
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: list(A.sources())
            [('horizontal-strebel', StrebelGraph("(0,~1:1,~0,1:1)"))]
        """
        return (x for x, backward_neighbors in self._backward_neighbors.items() if not backward_neighbors)

    def sinks(self):
        r"""
        Iterate through sinks (states with no forward neighbor) in this automaton.

        EXAMPLES::

            sage: from veerer import *
            sage: fp = "(0,2,1)(~0,3,~1)"
            sage: bdry = "(~2:2,~3:2)"
            sage: cols = "BRRR"
            sage: vt = VeeringTriangulation(fp, bdry, cols)

        In the Delaunay automaton, the sinks are the vertical-Strebel veering triangulations::

            sage: A = DelaunayAutomaton(backward=True)
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: list(A.sinks())
            [VeeringTriangulation("(0,~3,2)(1,3,~2)", boundary="(~1:2,~0:2)", colouring="RRRB")]
            sage: set(A.sinks()) == set(vt for vt in A if vt.is_strebel(VERTICAL))
            True

        In the Delaunay-Strebel, the sinks are the vertical Strebel graphs::

            sage: A = DelaunayStrebelAutomaton(backward=True)
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: list(A.sinks())
            [('vertical-strebel', StrebelGraph("(0,~1:1,~0,1:1)"))]
        """
        return (x for x, forward_neighbors in self._forward_neighbors.items() if not forward_neighbors)

    def transitions(self):
        r"""
        Run through the transitions of this automaton.
        """
        for s, out_neighbors in self._forward_neighbors.items():
            for t, label in out_neighbors:
                yield (s, t, label)

    def __contains__(self, state):
        r"""
        Return whether ``state`` is contained in the automaton.
        """
        return state.copy(mutable=False) in self._forward_neighbors

    # TODO: provide vertex_map and edge_map functions
    def to_graph(self, directed=True, multiedges=True, loops=True):
        r"""
        Return the underlying graph as a Sage graph or digraph.

        INPUT:

        - ``directed`` - boolean (default ``True``) - whether to make it directed

        - ``multiedges`` - boolean (default ``True``)

        - ``loops`` - boolean (default ``True``)

        EXAMPLES::

            sage: from veerer import *

            sage: fp = "(0,~7,6)(1,~8,~2)(2,~6,~3)(3,5,~4)(4,8,~5)(7,~1,~0)"
            sage: cols = 'RBRBRBBBB'
            sage: vt = VeeringTriangulation(fp, cols)
            sage: A = CoreAutomaton()
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: A
            Core veering automaton with 86 vertices

            sage: A.to_graph()
            Looped multi-digraph on 86 vertices

            sage: A.to_graph(directed=False, multiedges=False, loops=True)
            Looped graph on 86 vertices
        """
        if directed:
            from sage.graphs.digraph import DiGraph
            G = DiGraph(loops=loops, multiedges=multiedges)
        else:
            from sage.graphs.graph import Graph
            G = Graph(loops=loops, multiedges=multiedges)

        for g, neighb in self._forward_neighbors.items():
            for gg, label in neighb:
                G.add_edge(g, gg, label)

        return G

    # TODO: move or deprecate (not a generic method)
    def rotation_automorphism(self):
        r"""
        Return the automorphism of the vertices that corresponds to rotation.

        Note that this automorphism reverses the edge direction.

        EXAMPLES::

            sage: from veerer import *

            sage: vt = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton()
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: rot = A.rotation_automorphism()

            sage: G = A.to_graph()
            sage: all(G.has_edge(rot[b], rot[a]) for a,b in G.edges(sort=False, labels=False))
            True
        """
        aut = {}
        for vt in self._forward_neighbors:
            vt2 = vt.copy(mutable=True)
            vt2.rotate()
            vt2.set_canonical_labels()
            vt2.set_immutable()
            aut[vt] = vt2
        return aut

    # TODO: move or deprecate (not a generic method)
    def conjugation_automorphism(self):
        """
        Return the automorphism of the vertices that corresponds to complex conjugation.

        EXAMPLES::

            sage: from veerer import *

            sage: vt = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton()
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: conj = A.conjugation_automorphism()

            sage: G = A.to_graph()
            sage: all(G.has_edge(conj[a], conj[b]) for a,b in G.edges(sort=False, labels=False))
            True

        Conjugation and rotation commutes::

            sage: rot = A.rotation_automorphism()
            sage: all(rot[conj[a]] == conj[rot[a]] for a in G)
            True
        """
        aut = {}
        for vt in self._forward_neighbors:
            vt2 = vt.copy(mutable=True)
            vt2.conjugate()
            vt2.set_canonical_labels()
            vt2.set_immutable()
            aut[vt] = vt2
        return aut

    def export_dot(self, filename=None, triangulations=False):
        r"""
        Write dot data into the file ``filename``.

        To compile with graphviz, use


            $ PROG -Tpdf -o OUTPUT_FILE INPUT_FILE

        where

        - PROG is one of dot, neato, twopi, circo, fdp, sfdp, patchwork, osage
        - OUTPUT_FILE the name of the output file
        - INPUT_FILE the name of the dot file

        INPUT:

        - ``filename`` - an optional filename

        - ``triangulations`` - boolean - if set to ``True`` then an svg image is
          generated for each triangulation with proper link. If the dot file is
          compiled in svg format then clicking on a vertex will open the picture
          of the triangulation.

        EXAMPLES::

            sage: from veerer import *

            sage: filename = tmp_filename() + '.dot'
            sage: vt = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton()
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: A.export_dot(filename)
        """
        if filename is not None:
            if not filename.endswith('.dot'):
                raise ValueError('the filename should ends with .dot')
            f = open(filename, 'w')
            path = filename[:-4]
            rel_path = os.path.split(path)[-1]
        elif triangulations:
            raise ValueError
        else:
            f = sys.stdout

        if triangulations:
            os.mkdir(path)

        f.write('/**********************************************************************/\n')
        f.write('/* Automaton of core Veering triangulations                           */\n')
        f.write('/*                                                                    */\n')
        f.write('/* To compile the file to pdf run                                     */\n')
        f.write('/*                                                                    */\n')
        f.write('/*    $ sfdp -Tpdf -o file.pdf file.dot                               */\n')
        f.write('/*                                                                    */\n')
        seed_line = '/* seed: %s' % min(self._forward_neighbors)
        seed_line += ' ' * (70 - len(seed_line)) + '*/\n'
        f.write(seed_line)
        f.write('/*                                                                    */\n')
        f.write('/**********************************************************************/\n')
        f.write('\n')

        f.write('digraph MyGraph {\n')
        f.write(' node [shape=circle style=filled margin=0.1 width=0 height=0]\n')
        for T in self._forward_neighbors:
            g = T.to_string()
            if triangulations:
                t_filename = os.path.join(path, g + '.svg')
                t_rel_filename = os.path.join(rel_path, g + '.svg')
                F = T.flat_structure_min()
                F.set_pos(cylinders=T.cylinders(BLUE) + T.cylinders(RED))
                F.plot().save(t_filename)

            typ = T.properties_code()
            colour = PROPERTIES_COLOURS.get(typ, '#000000')

            aut_size = len(T.automorphisms())
            if triangulations:
                f.write("""    %s [label="%d" style=filled color="%s" tooltip="%s" href="%s"];\n""" % (g, aut_size, colour, g, t_rel_filename))
            else:
                f.write('    %s [label="%d" color="%s"];\n' % (g, aut_size, colour))

            for TT, flip_data in self._forward_neighbors[T]:
                gg = TT.to_string()
                # TODO: restore coloring options
                # f.write('    %s -> %s [color="%s;%f:%s"];\n' % (g, gg, old_col, 0.5, new_col))
                f.write('    %s -> %s;\n' % (g, gg))
        f.write('}\n')

        if filename is not None:
            f.close()

    # TODO: move or deprecate (not a generic method)
    def statistics(self):
        r"""
        Return detailed statistics about the properties of the veering
        triangulations.

        .. SEEALSO::

            :meth:`print_statistics`

        EXAMPLES::

            sage: from veerer import *

            sage: vt = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton()
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: st = A.statistics()
            sage: st
            {0: 24, 1: 4, 2: 4, 16: 28, 17: 5, 18: 5, 24: 10, 29: 3, 30: 3}
        """
        from collections import defaultdict
        d = defaultdict(int)
        for vt in self:
            d[vt.properties_code()] += 1
        return dict(d)

    # TODO: move or deprecate (not a generic method)
    def print_statistics(self, f=None):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: vt = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton()
            sage: A.add_seed(vt)
            1
            sage: A.run()
            0
            sage: A.print_statistics()
                    red square-tiled 3
                   blue square-tiled 3
            quadrangulable geometric 10
                       red geometric 5
                      blue geometric 5
                           geometric 28
                                 red 4
                                blue 4
                                none 24
        """
        if f is None:
            from sys import stdout
            f = stdout
        from veerer.constants import properties_to_string, key_property
        st = self.statistics()
        for k in sorted(st, key=key_property):
            v = st[k]
            f.write("%24s %d\n" % (properties_to_string(k), v))

    @classmethod
    def from_triangulation(cls, T, *args, **kwds):
        A = cls()
        A.add_seed(T)
        A.run(**kwds)
        return A

    # TODO: it is a bit absurd to have this in the generic class since this does not make
    # sense for triangulations
    # TODO: this is broken for DelaunayAutomaton with QuadraticStratum(8)
    # TODO: one should be more careful whether one wants the poles to be half edges or vertices
    @classmethod
    def from_stratum(self, stratum, **kwds):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import Stratum         # optional - surface_dynamics
            sage: CoreAutomaton.from_stratum(Stratum([2], 1))  # optional - surface_dynamics
            Core veering automaton with 86 vertices
            sage: DelaunayAutomaton.from_stratum(Stratum([2], 1))  # optional - surface_dynamics
            Delaunay automaton with 54 vertices

            sage: Q = Stratum([8], 2)                             # optional - surface_dynamics
            sage: A = CoreAutomaton.from_stratum(Q, max_size=100) # optional - surface_dynamics
            sage: A                                               # optional - surface_dynamics
            Partial core veering automaton with 101 vertices
        """
        return self.from_triangulation(VeeringTriangulation.from_stratum(stratum), **kwds)

    def out_neighbors(self, state):
        state = state.copy(mutable=False)
        neighbors = self._forward_neighbors.get(state, None)
        if neighbors is None:
            raise ValueError('state not in the automaton')
        return neighbors

    ######################
    # search implementation #
    ######################

    def add_seed(self, state, setup=True):
        r"""
        Add the seed ``state`` to the search.

        Return ``0`` if the state is already present in the graph and ``1``
        otherwise.
        """
        if setup:
            state = self._seed_setup(state)

        if state in self._forward_neighbors:
            if self._verbosity >= 2:
                print('[add_seed] state=%s already in the graph' % (state,))
            return 0

        if self._verbosity >= 2:
            print('[add_seed] adding state=%s to the list of seeds' % (state,))

        self._seeds.append(state)
        return 1

    def next_seed(self):
        r"""
        Return the next seed or ``None`` if there is none.

        This function modifies the `_seeds` and `_backward_flip_queue` attributes.
        """
        while self._seeds or self._backward_flip_queue:
            while self._seeds:
                seed = self._seeds.pop()
                if seed not in self._forward_neighbors:
                    # not explored forward yet
                    return seed

            if not self._backward_flip_queue:
                return

            state = self._backward_flip_queue.pop()
            if self._verbosity >= 2:
                print('[next_seed] found state to feed backward_flip_queue %s' % (state,))
            for back_neighbor, label in self._in_neighbors(state):
                if self._verbosity >= 2:
                    print('[next_seed] add back_neighbor %s' % (back_neighbor,))
                self.add_seed(back_neighbor, setup=False)

    def run(self, max_size=None):
        r"""
        Discover new states and transitions in the automaton.

        INPUT:

        - ``max_size`` -- an optional bound on the number of new states to
          compute

        EXAMPLES::

            sage: from veerer import *

        The torus::

            sage: T = VeeringTriangulation("(0,2,~1)(1,~0,~2)", "RBB")
            sage: A = CoreAutomaton()
            sage: A.add_seed(T)
            1
            sage: A.run()
            0
            sage: A
            Core veering automaton with 2 vertices

            sage: A = ReducedCoreAutomaton()
            sage: A.add_seed(T)
            1
            sage: A.run()
            0
            sage: A
            Reduced core veering automaton with 1 vertex

        A more complicated surface in Q_1(1^2, -1^2)::

            sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
            sage: cols = 'BRBRBBBRBR'
            sage: T = VeeringTriangulation(fp, cols)

            sage: C = CoreAutomaton()
            sage: C.add_seed(T)
            1
            sage: C.run(10)
            1
            sage: C
            Partial core veering automaton with 10 vertices
            sage: C.run()
            0
            sage: C
            Core veering automaton with 1074 vertices

            sage: C = ReducedCoreAutomaton()
            sage: C.add_seed(T)
            1
            sage: C.run(10)
            1
            sage: C
            Partial reduced core veering automaton with 10 vertices
            sage: C.run()
            0
            sage: C
            Reduced core veering automaton with 356 vertices

        TESTS::

            sage: from veerer import VeeringTriangulation, CoreAutomaton
            sage: T = VeeringTriangulation("(0,2,~1)(1,~0,~2)", "RBB")
            sage: A = CoreAutomaton()
            sage: A.add_seed(T)
            1
            sage: A
            Partial core veering automaton with 0 vertex
            sage: A.run(1)
            1
            sage: A
            Partial core veering automaton with 1 vertex
            sage: A.run(1)
            1
            sage: A
            Partial core veering automaton with 2 vertices
            sage: A.run(1)
            0
            sage: A
            Core veering automaton with 2 vertices
        """
        forward_neighbors = self._forward_neighbors
        backward_neighbors = self._backward_neighbors
        backward_flip_queue = self._backward_flip_queue
        branch = collections.deque(self._branch)

        count = 0
        old_size = len(self._forward_neighbors)
        while max_size is None or count < max_size:
            if self._verbosity >= 2:
                print('[automaton] new loop')
                sys.stdout.flush()

            while not branch:
                # forward search is over... find a new seed
                state = self.next_seed()
                if state is None:
                    # no seed available anymore
                    if self._verbosity >= 2:
                        print('[automaton] done')
                        sys.stdout.flush()
                    self._branch = []
                    return 0

                if self._verbosity >= 2:
                    print('[automaton] seed %s' % (state,))
                    sys.stdout.flush()

                forward_neighbors[state] = []
                backward_neighbors[state] = []
                if self._backward:
                    self._backward_flip_queue.append(state)
                branch.clear()
                branch.append(state)
                count += 1

                if max_size is not None and count >= max_size:
                    self._branch = list(branch)
                    return 1

            # next step of forward search
            if self._method == 'DFS':
                state = branch.pop()
            elif self._method == 'BFS':
                state = branch.popleft()
            else:
                raise RuntimeError

            if self._verbosity >= 2:
                print('[automaton] at %s' % (state,))
                sys.stdout.flush()

            for out_neighbor, label in self._out_neighbors(state):
                if out_neighbor not in forward_neighbors:
                    assert out_neighbor not in backward_neighbors
                    forward_neighbors[out_neighbor] = []
                    backward_neighbors[out_neighbor] = []
                    branch.append(out_neighbor)

                    if self._backward:
                        self._backward_flip_queue.append(out_neighbor)

                    count += 1

                forward_neighbors[state].append((out_neighbor, label))
                backward_neighbors[out_neighbor].append((state, label))

        self._branch = list(branch)
        return 1

    ############################################################
    # Custom methods that have to be implemented in subclasses #
    ############################################################

    def _setup(self, **extra_kwds):
        r"""
        Handling of extra setup arguments (called once during ``__init__``).
        """
        pass

    def _seed_setup(self, state):
        return state

    def _out_neighbors(self, state):
        r"""
        Return an iterable of out-neighbors of ``state``.

        Each element consists of a pair ``(new_state, edge_label)``.
        """
        raise NotImplementedError

    def _in_neighbors(self, state):
        r"""
        Return an iterable of in-neighbors of ``state``.

        Each element consists of a pair ``(new_state, edge_label)``.
        """
        raise NotImplementedError


# TODO: add the examples of the triangulations of the n-gon. To match with
# the literature, one would need to consider canonical labels where we do
# not modify edge names on the boundary.
class FlipGraph(Automaton):
    r"""
    The flip graph of triangulations.

    EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation([(0,1,2),(-1,-2,-3)])
            sage: A = FlipGraph()
            sage: A.add_seed(T)
            1
            sage: A.run()
            0
            sage: A
            Triangulation automaton with 1 vertex

            sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
            sage: T = Triangulation(fp)
            sage: A = FlipGraph()
            sage: A.add_seed(T)
            1
            sage: A.run()
            0
            sage: A
            Triangulation automaton with 236 vertices
    """
    _name = 'triangulation'

    def _seed_setup(self, state):
        state = state.copy(mutable=True)
        state.set_canonical_labels()
        state.set_immutable()
        return state

    def _out_neighbors(self, state):
        for e in state.flippable_edges():
            new_state = state.copy(mutable=True)
            new_state.flip(e)
            new_state.set_canonical_labels()
            new_state.set_immutable()
            yield (new_state, e)

    _in_neighbors = _out_neighbors


class CoreAutomaton(Automaton):
    r"""
    Automaton of core veering triangulations.
    """
    _name = 'core veering'

    def _seed_setup(self, state):
        state = state.copy(mutable=True)
        state.set_canonical_labels()
        state.set_immutable()
        return state

    def _out_neighbors(self, state):
        state = state.copy(mutable=True)
        for e in state.forward_flippable_edges():
            for col in (BLUE, RED):
                old_col = state.colour(e)
                state.flip(e, col, check=CHECK)
                if state.edge_has_curve(e):
                    new_state = state.copy(mutable=True)
                    new_state.set_canonical_labels()
                    new_state.set_immutable()
                    yield (new_state, (e, col))
                state.flip_back(e, old_col, check=CHECK)

    def _in_neighbors(self, state):
        state = state.copy(mutable=True)
        for e in state.backward_flippable_edges():
            for col in (BLUE, RED):
                old_col = state.colour(e)
                state.flip_back(e, col, check=CHECK)
                if state.edge_has_curve(e):
                    new_state = state.copy(mutable=True)
                    new_state.set_canonical_labels()
                    new_state.set_immutable()
                    yield (new_state, (e, col))
                state.flip(e, col, check=CHECK)


class ReducedCoreAutomaton(Automaton):
    r"""
    Automaton of reduced core veering triangulations.

    This automaton can equivalently be thought as the (birecurrent) train track
    automaton in a fixed stratum.

    EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,2,~1)(1,~0,~2)", "RBB")
            sage: A = ReducedCoreAutomaton()
            sage: A.add_seed(T)
            1
            sage: A.run()
            0
            sage: A
            Reduced core veering automaton with 1 vertex

        A more complicated surface in Q_1(1^2, -1^2)::

            sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
            sage: cols = 'BRBRBBBRBR'
            sage: T = VeeringTriangulation(fp, cols)

            sage: C_BFS = ReducedCoreAutomaton(method='BFS')
            sage: C_BFS.add_seed(T)
            1
            sage: C_BFS.run()
            0
            sage: C_BFS
            Reduced core veering automaton with 356 vertices

            sage: C_DFS = ReducedCoreAutomaton(method='DFS')
            sage: C_DFS.add_seed(T)
            1
            sage: C_DFS.run()
            0
            sage: C_DFS
            Reduced core veering automaton with 356 vertices
    """
    _name = 'reduced core veering'

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulation):
            raise TypeError
        state = state.copy(mutable=True)
        state.forgot_forward_flippable_colour()
        state.set_canonical_labels()
        state.set_immutable()
        return state

    # TODO: move this method to VeeringTriangulation
    def _flip(self, state, e, col, check=CHECK):
        r"""
        Implement the edge flip for reduced veering triangulations.

        The flip is the same as for veering triangulation but one additional
        step is needed to turn some of the edges into purple edges.
        """
        if CHECK:
            assert state._mutable
            assert state.is_forward_flippable(e)

        state.flip(e, col, reduced=False, check=CHECK)

        if not state.edge_has_curve(e):
            return (False, ())

        recolorings = []
        a, b, c, d = state.square_about_edge(e)

        # assertions to be removed
        assert state._colouring[a] == RED, (a, colour_to_string(state._colouring[a]))
        assert state._colouring[b] == BLUE, (b, colour_to_string(state._colouring[b]))
        assert state._colouring[c] == RED, (c, colour_to_string(state._colouring[c]))
        assert state._colouring[d] == BLUE, (d, colour_to_string(state._colouring[d]))

        if col == BLUE:
            if state.is_forward_flippable(b):
                recolorings.append((b, BLUE))
                state._colouring[b] = PURPLE
                state._colouring[state._ep[b]] = PURPLE
            if b != d and state.is_forward_flippable(d):
                recolorings.append((d, BLUE))
                state._colouring[d] = PURPLE
                state._colouring[state._ep[d]] = PURPLE

            # assertions to be removed
            if CHECK:
                assert not state.is_forward_flippable(e)
                assert not state.is_forward_flippable(a)
                assert not state.is_forward_flippable(c)

        elif col == RED:
            if state.is_forward_flippable(a):
                recolorings.append((a, RED))
                state._colouring[a] = PURPLE
                state._colouring[state._ep[a]] = PURPLE
            if a != c and state.is_forward_flippable(c):
                recolorings.append((c, RED))
                state._colouring[c] = PURPLE
                state._colouring[state._ep[c]] = PURPLE

            if CHECK:
                assert not state.is_forward_flippable(e)
                assert not state.is_forward_flippable(b)
                assert not state.is_forward_flippable(d)

        return (True, recolorings)

    def _flip_back(self, state, e, recolorings, check=CHECK):
        for ee, ccol in recolorings:
            state._colouring[ee] = ccol
            state._colouring[state._ep[ee]] = ccol
        state.flip_back(e, PURPLE, check=CHECK)

    def _out_neighbors(self, state):
        state = state.copy(mutable=True)
        for e in state.purple_edges():
            for col in (BLUE, RED):
                status, recolorings = self._flip(state, e, col, check=CHECK)
                if state.edge_has_curve(e):
                    new_state = state.copy(mutable=True)
                    new_state.set_canonical_labels()
                    new_state.set_immutable()
                    yield (new_state, (e, col))
                self._flip_back(state, e, recolorings, check=CHECK)

    # TODO: implement in_neighbors


# TODO: this should be renamed DelaunayAutomaton
class DelaunayAutomaton(Automaton):
    r"""
    Automaton of Delaunay (veering) triangulations.

    A veering triangulation is Delaunay (for the L-infinity metric), if the set
    of associated L^oo-vector data is full dimensional in the ambient stratum.

    EXAMPLES::

        sage: from veerer import *

        sage: fp = "(0,~7,6)(1,~8,~2)(2,~6,~3)(3,5,~4)(4,8,~5)(7,~1,~0)"
        sage: cols = "RBRBRBBBB"
        sage: vt = VeeringTriangulation(fp, cols)
        sage: A = DelaunayAutomaton()
        sage: A.add_seed(vt)
        1
        sage: A.run()
        0
        sage: A
        Delaunay automaton with 54 vertices

    One can check that the cardinality is indeed correct::

        sage: C = CoreAutomaton()
        sage: C.add_seed(vt)
        1
        sage: C.run()
        0
        sage: sum(x.is_delaunay() for x in C)
        54

    A meromorphic example in H(2,-2) that consists of three veering
    triangulations. Since the automaton has sinks and sources, one
    needs to run it with the option ``backward=True`` to get the
    full connected component that is made of 3 states::

        sage: fp = "(0,2,1)(~0,3,~1)"
        sage: bdry = "(~2:2,~3:2)"
        sage: cols0 = "BRRR"
        sage: cols1 = "BBRR"
        sage: cols2 = "RBRR"
        sage: vt0 = VeeringTriangulation(fp, bdry, cols0)
        sage: vt1 = VeeringTriangulation(fp, bdry, cols1)
        sage: vt2 = VeeringTriangulation(fp, bdry, cols2)

        sage: A0 = DelaunayAutomaton()
        sage: A0.add_seed(vt0)
        1
        sage: A0.run()
        0
        sage: A0
        Delaunay automaton with 3 vertices

        sage: A1 = DelaunayAutomaton()
        sage: A1.add_seed(vt1)
        1
        sage: A1.run()
        0
        sage: A1
        Delaunay automaton with 2 vertices

        sage: A2 = DelaunayAutomaton()
        sage: A2.add_seed(vt2)
        1
        sage: A2.run()
        0
        sage: A2
        Delaunay automaton with 1 vertex

        sage: B0 = DelaunayAutomaton(backward=True)
        sage: B0.add_seed(vt0)
        1
        sage: B0.run()
        0
        sage: B0
        Delaunay automaton with 3 vertices

        sage: B1 = DelaunayAutomaton(backward=True)
        sage: B1.add_seed(vt1)
        1
        sage: B1.run()
        0
        sage: B1
        Delaunay automaton with 3 vertices

        sage: B2 = DelaunayAutomaton(backward=True)
        sage: B2.add_seed(vt2)
        1
        sage: B2.run()
        0
        sage: B2
        Delaunay automaton with 3 vertices

    An example with linear constraint : the L-shape surface in the stratum H(2)
    made of three squares::

        sage: vt, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
        sage: f = VeeringTriangulationLinearFamily(vt, [s, t])
        sage: A = DelaunayAutomaton()
        sage: A.add_seed(f)
        1
        sage: A.run()
        0
        sage: A
        Delaunay automaton with 6 vertices

    Some more L-shape surfaces::

        sage: for n in range(3, 7):
        ....:     vt, s, t = VeeringTriangulations.L_shaped_surface(1, n-2, 1, 1)
        ....:     f = VeeringTriangulationLinearFamily(vt, [s, t])
        ....:     A = DelaunayAutomaton()
        ....:     _ = A.add_seed(f)
        ....:     _ = A.run()
        ....:     print('n={}: {}'.format(n, A))
        n=3: Delaunay automaton with 6 vertices
        n=4: Delaunay automaton with 86 vertices
        n=5: Delaunay automaton with 276 vertices
        n=6: Delaunay automaton with 800 vertices
    """
    _name = 'Delaunay'

    def _setup(self, backend=None):
        self._backend = backend

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulation):
            raise TypeError('invalid seed of type {}'.format(type(state).__name__))
        if self._backend is None:
            from .polyhedron.cone import default_backend
            self._backend = default_backend(state.base_ring())
        if not state.is_delaunay(backend=self._backend):
            raise ValueError('seed is not Delaunay')
        state = state.copy(mutable=True)
        state.set_canonical_labels()
        state.set_immutable()
        return state

    def _out_neighbors(self, state, check=CHECK):
        r"""
        Run through the list of out neighbors.
        """
        for edges, col in state.delaunay_flips(backend=self._backend):
            assert all(state.colour(e) == state.colour(edges[0]) for e in edges)
            out_neighbor = state.copy(mutable=True)
            for e in edges:
                out_neighbor.flip(e, col, check=CHECK)
            if CHECK:
                out_neighbor._check(RuntimeError)
                if not out_neighbor.is_delaunay(backend=self._backend):
                    raise RuntimeError
            out_neighbor.set_canonical_labels()
            out_neighbor.set_immutable()
            yield (out_neighbor, (edges, col))

    def _in_neighbors(self, state, check=CHECK):
        """
        Run through the list of in neighbors.
        """
        for edges, col in state.backward_delaunay_flips(backend=self._backend):
            assert all(state.colour(e) == state.colour(edges[0]) for e in edges)
            in_neighbor = state.copy(mutable=True)
            for e in edges:
                in_neighbor.flip_back(e, col, check=CHECK)
                if CHECK:
                    in_neighbor._check(RuntimeError)
                    if not in_neighbor.is_delaunay(backend=self._backend):
                        raise RuntimeError
            in_neighbor.set_canonical_labels()
            in_neighbor.set_immutable()
            yield (in_neighbor, (edges, col))

    def cylinder_diagrams(self, col=RED):
        r"""
        A cylinder diagram is a cylindrical veering triangulation up to twist action.

        EXAMPLES::

            sage: from veerer import RED, BLUE, DelaunayAutomaton
            sage: from surface_dynamics import Stratum                       # optional - surface_dynamics
            sage: A = DelaunayAutomaton.from_stratum(Stratum([2], 1))        # optional - surface_dynamics
            sage: len(A.cylinder_diagrams(RED))                              # optional - surface_dynamics
            2
            sage: A = DelaunayAutomaton.from_stratum(Stratum([1, 1], 1))     # optional - surface_dynamics
            sage: len(A.cylinder_diagrams(RED))                              # optional - surface_dynamics
            4

        For Veech surfaces, the cylinder diagrams are in bijection with cusps::

            sage: from veerer import VeeringTriangulations, VeeringTriangulationLinearFamily
            sage: t, x, y = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: f = VeeringTriangulationLinearFamily(t, [x, y])
            sage: A = f.delaunay_automaton()
            sage: len(A.cylinder_diagrams())
            2

        We check that in H(2) prototypes coincide with cylinder diagrams made of two cylinders::

            sage: from veerer.linear_family import VeeringTriangulationLinearFamilies
            sage: for D in [5, 8, 9, 12, 13, 16]:  # long time
            ....:     if D % 8 == 1:
            ....:         spins = [0, 1]
            ....:     else:
            ....:         spins = [None]
            ....:     for spin in spins:
            ....:         prototypes = list(VeeringTriangulationLinearFamilies.H2_prototype_parameters(D, spin))
            ....:         if not prototypes:
            ....:             continue
            ....:         a, b, c, e = prototypes[0]
            ....:         F = VeeringTriangulationLinearFamilies.prototype_H2(a, b, c, e)
            ....:         A = F.delaunay_automaton()
            ....:         cylsRED = A.cylinder_diagrams(RED)
            ....:         n1RED = sum(len(next(iter(vts)).cylinders(RED)) == 1 for vts in cylsRED)
            ....:         n2RED = sum(len(next(iter(vts)).cylinders(RED)) == 2 for vts in cylsRED)
            ....:         cylsBLUE = A.cylinder_diagrams(BLUE)
            ....:         n1BLUE = sum(len(next(iter(vts)).cylinders(BLUE)) == 1 for vts in cylsBLUE)
            ....:         n2BLUE = sum(len(next(iter(vts)).cylinders(BLUE)) == 2 for vts in cylsBLUE)
            ....:         assert n1RED == n1BLUE
            ....:         assert n2RED == n2BLUE
            ....:         assert n2RED == len(prototypes)
            ....:         print(D, spin, n1RED, n2RED)
            5 None 0 1
            8 None 0 2
            9 0 1 1
            12 None 0 3
            13 None 0 3
            16 None 1 2
        """
        cylindricals = set()
        orbits = []
        for x in self._forward_neighbors:
            if x in cylindricals or not x.is_cylindrical(col):
                continue
            orbit = set()
            todo = [x]
            while todo:
                x = todo.pop()
                orbit.add(x)
                for y, (edges, new_col) in self._forward_neighbors[x]:
                    if new_col == col and y not in orbit:
                        assert y.is_cylindrical(col)
                        orbit.add(y)
                        todo.append(y)
            orbits.append(orbit)
            if not all(x not in cylindricals for x in orbit):
                print("PROBLEM")
            cylindricals.update(orbit)
        return orbits


# TODO: for now the implementation is only partial since we are missing the
# transitions from Strebel graphs to Delauany triangulations
class DelaunayStrebelAutomaton(Automaton):
    r"""
    Delaunay-Strebel automaton.

    The states of the Delaunay-Strebel automaton are Delaunay triangulations,
    vertical Strebel graphs and horizontal Strebel graphs. They are stored
    as pairs ``(kind, state)`` where

    - ``kind`` is either one of the strings ``'delaunay'``,
      ``'vertical-strebel'`` or ``'horizontal-strebel'``
    - ``state`` is either a :class:`VeeringTriangulation` or
      :class:`StrebelGraph`

    EXAMPLES:

    The example of the meromorphic stratum H(2,-2) that consists of 9 Delaunay
    veering triangulations, 1 horizontal Strebel graph and 1 vertical Strebel
    graph::

        sage: from veerer import *

        sage: fp = "(0,2,1)(~0,3,~1)"
        sage: bdry = "(~2:2,~3:2)"
        sage: cols0 = "BBRR"
        sage: cols1 = "RRBB"
        sage: vt0 = VeeringTriangulation(fp, bdry, cols0)
        sage: vt1 = VeeringTriangulation(fp, bdry, cols1)
        sage: DS = DelaunayStrebelAutomaton(backward=True)
        sage: DS.add_seed(vt0)
        1
        sage: DS.add_seed(vt1)
        1
        sage: DS.run()
        0
        sage: DS
        Delaunay-Strebel automaton with 11 vertices
        sage: sum(kind == 'vertical-strebel' for kind, state in DS)
        1
        sage: sum(kind == 'horizontal-strebel' for kind, state in DS)
        1
        sage: sum(kind == 'delaunay' for kind, state in DS)
        9

    Some one and two dimensional examples in genus 0::

        sage: examples = [StrebelGraph("(0,~1)(1,~0)"), StrebelGraph("(0:2,~1)(1,~0)"),
        ....:             StrebelGraph("(0:2,~1)(1,~0:2)"), StrebelGraph("(0,~1)(1)(~0)"),
        ....:             StrebelGraph("(0:2,~1)(1)(~0)"), StrebelGraph("(0,~0,~1:1,1)"),
        ....:             StrebelGraph("(0,~0,~1)(1)"), StrebelGraph("(0,2,~0,~1)(1)(~2)"),
        ....:             StrebelGraph("(0,2,~1)(1)(~2,~0)"), StrebelGraph("(0:2,2,~1)(1,~0)(~2)")]
        sage: for G in examples:
        ....:     print(G)
        ....:     DS = DelaunayStrebelAutomaton(backward=True)
        ....:     _ = DS.add_seed(G)
        ....:     _ = DS.run()
        ....:     DS._check()
        ....:     print(DS)
        ....:     n_hs = sum(kind == 'horizontal-strebel' for kind, state in DS)
        ....:     n_vs = sum(kind == 'vertical-strebel' for kind, state in DS)
        ....:     n_d = sum(kind == 'delaunay' for kind, state in DS)
        ....:     print(n_hs, n_vs, n_d)
        StrebelGraph("(0,~1)(1,~0)")
        Delaunay-Strebel automaton with 10 vertices
        1 1 8
        StrebelGraph("(0:2,~1)(1,~0)")
        Delaunay-Strebel automaton with 26 vertices
        3 3 20
        StrebelGraph("(0:2,~1)(1,~0:2)")
        Delaunay-Strebel automaton with 16 vertices
        2 2 12
        StrebelGraph("(0,~1)(1)(~0)")
        Delaunay-Strebel automaton with 6 vertices
        1 1 4
        StrebelGraph("(0:2,~1)(1)(~0)")
        Delaunay-Strebel automaton with 20 vertices
        3 3 14
        StrebelGraph("(0,~0,~1:1,1)")
        Delaunay-Strebel automaton with 13 vertices
        1 1 11
        StrebelGraph("(0,~0,~1)(1)")
        Delaunay-Strebel automaton with 10 vertices
        1 1 8
        StrebelGraph("(0,2,~0,~1)(1)(~2)")
        Delaunay-Strebel automaton with 34 vertices
        2 2 30
        StrebelGraph("(0,2,~1)(1)(~2,~0)")
        Delaunay-Strebel automaton with 52 vertices
        1 1 50
        StrebelGraph("(0:2,2,~1)(1,~0)(~2)")
        Delaunay-Strebel automaton with 126 vertices
        6 6 114
    """
    _name = 'Delaunay-Strebel'

    def _check(self):
        super()._check()

        from .features import surface_dynamics_feature
        if surface_dynamics_feature.is_present():
            s = set(state.stratum() for kind, state in self)
            if len(s) != 1:
                raise ValueError('got different strata: {}'.format(sorted(s)))

        horizontal_strebel = set()
        vertical_strebel = set()
        for kind, state in self:
            if kind == 'horizontal-strebel':
                horizontal_strebel.add(state)
            elif kind == 'vertical-strebel':
                vertical_strebel.add(state)

        assert horizontal_strebel == vertical_strebel, horizontal_strebel.symmetric_difference(vertical_strebel)

    def _setup(self, backend=None):
        self._backend = backend

    def _seed_setup(self, state):
        if isinstance(state, tuple) and len(state) == 2:
            kind, state = state
        elif isinstance(state, VeeringTriangulation):
            kind = 'delaunay'
        elif isinstance(state, StrebelGraph):
            kind = 'vertical-strebel'
        else:
            raise TypeError('invalid state')

        if self._backend is None:
            from .polyhedron.cone import default_backend
            self._backend = default_backend(state.base_ring())

        if kind == 'delaunay':
            if not state.is_geometric(backend=self._backend):
                raise ValueError('invalid seed: non-Delaunay veering triangulation {}'.format(state))

        state = state.copy(mutable=True)
        state.set_canonical_labels()
        state.set_immutable()
        return (kind, state)

    def _out_neighbors(self, state, check=CHECK):
        r"""
        Return the list of out-neighbors.
        """
        kind, state = state
        if kind == 'delaunay':
            flips = state.delaunay_flips(backend=self._backend)
            if not flips:
                if CHECK:
                    assert state.is_strebel(VERTICAL)
                out_neighbor = state.strebel_graph(VERTICAL, mutable=True)
                out_neighbor.set_canonical_labels()
                out_neighbor.set_immutable()
                yield (('vertical-strebel', out_neighbor), 'strebel')
                return

            for edges, col in flips:
                assert all(state.colour(e) == state.colour(edges[0]) for e in edges)
                out_neighbor = state.copy(mutable=True)
                for e in edges:
                    out_neighbor.flip(e, col, check=CHECK)
                if CHECK:
                    out_neighbor._check(RuntimeError)
                    if not out_neighbor.is_delaunay(backend=self._backend):
                        raise RuntimeError
                out_neighbor.set_canonical_labels()
                out_neighbor.set_immutable()
                yield ((kind, out_neighbor), (edges, col))

        elif kind == 'vertical-strebel':
            return

        elif kind == 'horizontal-strebel':
            for colouring in state.colourings():
                for vt in state.delaunay_triangulations(colouring, HORIZONTAL, mutable=True):
                    vt.set_canonical_labels()
                    vt.set_immutable()
                    yield (('delaunay', vt), colouring)

        else:
            raise RuntimeError

    def _in_neighbors(self, state):
        """
        Return the list of Delaunay backward flippable edges from ``state``.
        """
        kind, state = state
        if kind == 'delaunay':
            flips = state.backward_delaunay_flips(backend=self._backend)
            if not flips:
                if CHECK:
                    assert state.is_strebel(HORIZONTAL)
                in_neighbor = state.strebel_graph(HORIZONTAL, mutable=True)
                in_neighbor.set_canonical_labels()
                in_neighbor.set_immutable()
                yield (('horizontal-strebel', in_neighbor), 'strebel')
                return

            for edges, col in flips:
                assert all(state.colour(e) == state.colour(edges[0]) for e in edges)
                in_neighbor = state.copy(mutable=True)
                for e in edges:
                    in_neighbor.flip_back(e, col, check=CHECK)
                    if CHECK:
                        in_neighbor._check(RuntimeError)
                        if not in_neighbor.is_delaunay(backend=self._backend):
                            raise RuntimeError
                in_neighbor.set_canonical_labels()
                in_neighbor.set_immutable()
                yield ((kind, in_neighbor), (edges, col))

        elif kind == 'vertical-strebel':
            for colouring in state.colourings():
                for vt in state.delaunay_triangulations(colouring, VERTICAL, mutable=True):
                    vt.set_canonical_labels()
                    vt.set_immutable()
                    yield (('delaunay', vt), colouring)

        elif kind == 'horizontal-strebel':
            return

        else:
            raise RuntimeError
