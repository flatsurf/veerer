r"""
Train-track and triangulations automata.
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2018 Mark Bell
#                     2018-2023 Vincent Delecroix
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

import os
import sys
from time import time

from .triangulation import Triangulation
from .veering_triangulation import VeeringTriangulation
from .linear_family import VeeringTriangulationLinearFamily
from .constants import RED, BLUE, PURPLE, PROPERTIES_COLOURS, colour_to_char, colour_to_string
from .permutation import perm_invert

# TODO: when set to True a lot of intermediate checks are performed
CHECK = False


class Automaton(object):
    r"""
    Automaton of veering triangulations.

    EXAMPLES::

        sage: from veerer import *  # random output due to deprecation warnings from realalg

        sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
        sage: A = CoreAutomaton(T)
        sage: len(A)
        2

        sage: T = VeeringTriangulation('(0,1,2)', 'BBR')
        sage: A = CoreAutomaton(T)
        sage: len(A)
        2

    The examples below is Q(1^2, -1^2) with two folded edges::

        sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
        sage: cols = 'BRBRBBBRBR'
        sage: T = VeeringTriangulation(fp, cols)
        sage: CoreAutomaton(T)
        Core veering automaton with 1074 vertices
        sage: ReducedCoreAutomaton(T)
        Reduced core veering automaton with 356 vertices

    Exploring strata::

        sage: from surface_dynamics import *                             # optional - surface_dynamics
        sage: strata = [AbelianStratum(2), QuadraticStratum(2,-1,-1),    # optional - surface_dynamics
        ....:           QuadraticStratum(2,2), AbelianStratum(1,1)]
        sage: for stratum in strata:                                     # optional - surface_dynamics
        ....:     print(stratum)
        ....:     vt = VeeringTriangulation.from_stratum(stratum)
        ....:     print(CoreAutomaton(vt))
        ....:     print(ReducedCoreAutomaton(vt))
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
    _name = ''

    def __init__(self, seed=None, verbosity=0, **extra_kwds):
        # verbosity level
        self._verbosity = 0

        # we encode the graph in a dictionary where the keys are the states
        # and the outgoing edges from a given state are encoded in the
        # corresponding value by a list of pairs (end_state, flip_data)
        self._graph = {}

        # current state (mutable triangulations)
        self._state = None

        # current branches for running the DFS
        self._state_branch = []
        self._flip_branch = []

        # the flips and relabelling performed along the branch
        self._flips = []
        self._relabellings = []

        self._setup(**extra_kwds)

        if seed is not None:
            self.set_seed(seed)
            self.run()

    def __str__(self):
        status = None
        if not self._flip_branch:
            status = 'uninitialized '
        elif len(self._flip_branch) != 1 or self._flip_branch[-1]:
            status = 'partial '
        else:
            status = ''
        return ("%s%s automaton with %s %s" % (status, self._name, len(self._graph), 'vertex' if len(self._graph) <= 1 else 'vertices')).capitalize()

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self._graph)

    def one_triangulation(self):
        return next(iter(self))

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
            sage: sum(T.is_geometric() for T in A)
            54
            sage: sum(T.is_cylindrical() for T in A)
            24
        """
        return iter(self._graph)

    def __contains__(self, item):
        if not isinstance(item, Triangulation):
            return False
        if item._mutable:
            item = item.copy(mutable=False)
        return item in self._graph

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
            sage: T = VeeringTriangulation(fp, cols)
            sage: A = CoreAutomaton(T)
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

        for g, neighb in self._graph.items():
            for gg, flip_data in neighb:
                G.add_edge(g, gg, flip_data)

        return G

    def rotation_automorphism(self):
        r"""
        Return the automorphism of the vertices that corresponds to rotation.

        Note that this automorphism reverses the edge direction.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
            sage: rot = A.rotation_automorphism()

            sage: G = A.to_graph()
            sage: all(G.has_edge(rot[b], rot[a]) for a,b in G.edges(sort=False, labels=False))
            True
        """
        aut = {}
        for vt in self._graph:
            vt2 = vt.copy(mutable=True)
            vt2.rotate()
            vt2.set_canonical_labels()
            vt2.set_immutable()
            aut[vt] = vt2
        return aut

    def conjugation_automorphism(self):
        """
        Return the automorphism of the vertices that corresponds to complex conjugation.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
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
        for vt in self._graph:
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
            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
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
        seed_line = '/* seed: %s' % min(self._graph)
        seed_line += ' ' * (70 - len(seed_line)) + '*/\n'
        f.write(seed_line)
        f.write('/*                                                                    */\n')
        f.write('/**********************************************************************/\n')
        f.write('\n')

        f.write('digraph MyGraph {\n')
        f.write(' node [shape=circle style=filled margin=0.1 width=0 height=0]\n')
        for T in self._graph:
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

            for TT, flip_data in self._graph[T]:
                gg = TT.to_string()
                # TODO: restore coloring options
                # f.write('    %s -> %s [color="%s;%f:%s"];\n' % (g, gg, old_col, 0.5, new_col))
                f.write('    %s -> %s;\n' % (g, gg))
        f.write('}\n')

        if filename is not None:
            f.close()

    def triangulations(self):
        r"""
        Return an iterator over the veering triangulations in this automaton.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
            sage: for t in A.triangulations():
            ....:     assert t.angles() == [6]
        """
        return iter(self)

    def num_triangulations(self):
        r"""
        Return the number of triangulations (= states of the automaton).

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
            sage: A.num_triangulations()
            86
        """
        return len(self._graph)

    num_states = num_triangulations

    def num_transitions(self):
        r"""
        Return the number of transitions of this automaton.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
            sage: A.num_transitions()
            300
        """
        return sum(len(flips) for flips in self._graph.values())

    def statistics(self):
        r"""
        Return detailed statistics about the properties of the veering
        triangulations.

        .. SEEALSO::

            :meth:`print_statistics`

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
            sage: st = A.statistics()
            sage: st
            {0: 24, 1: 4, 2: 4, 16: 28, 17: 5, 18: 5, 24: 10, 29: 3, 30: 3}
        """
        from collections import defaultdict
        d = defaultdict(int)
        for vt in self:
            d[vt.properties_code()] += 1
        return dict(d)

    def print_statistics(self, f=None):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
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
    def from_triangulation(self, T, reduced=False, max_size=None, verbosity=0):
        if reduced:
            A = ReducedCoreAutomaton(verbosity=verbosity)
        else:
            A = CoreAutomaton(verbosity=verbosity)
        A.set_seed(T)
        A.run(max_size)
        return A

    @classmethod
    def from_stratum(self, stratum, reduced=False, **kwds):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import *                 # optional - surface_dynamics
            sage: CoreAutomaton.from_stratum(AbelianStratum(2))  # optional - surface_dynamics
            Core veering automaton with 86 vertices

            sage: Q = QuadraticStratum(8)                         # optional - surface_dynamics
            sage: A = CoreAutomaton.from_stratum(Q, max_size=100) # optional - surface_dynamics
            sage: A                                               # optional - surface_dynamics
            Partial core veering automaton with 101 vertices
        """
        return self.from_triangulation(VeeringTriangulation.from_stratum(stratum), reduced=reduced, **kwds)

    ######################
    # DFS implementation #
    ######################

    def _setup(self):
        pass

    def _seed_setup(self, state):
        raise NotImplementedError

    def _forward_flips(self, state):
        r"""
        Return the list of forward flip_data from ``state``.

        To be implemented in subclasses.
        """
        raise NotImplementedError

    def _flip(self, flip_data):
        r"""
        Perform the flips given in ``flip_data``.

        To be implemented in subclasses
        """
        raise NotImplementedError

    def _flip_back(self, flip_back_data):
        r"""
        Perform back the flips given in ``flip_back_data``

        To be implemented in subclasses
        """
        raise NotImplementedError

    def set_seed(self, state):
        if self._state is not None:
            raise ValueError('seed already set')

        state = self._seed_setup(state)

        # TODO: check to be removed
        if CHECK:
            state._check()

        if self._verbosity:
            print('[automaton] stratum: %s' % state.stratum())
            print('[automaton] stratum dimension: %d' % state.stratum_dimension())
        if self._verbosity == 1:
            print('[automaton] size(graph)   size(path)    time')

        self._state = state
        immutable_state = state.copy(mutable=False)

        self._graph[immutable_state] = []

        self._state_branch.append(immutable_state)
        self._flip_branch.append(self._forward_flips(self._state))

    def run(self, max_size=None):
        r"""
        Discover the automaton via depth first search.

        INPUT:

        - ``max_size`` -- an optional bound on the number of new states to compute

        EXAMPLES::

            sage: from veerer import *

        The torus::

            sage: T = VeeringTriangulation("(0,2,~1)(1,~0,~2)", "RBB")
            sage: A = CoreAutomaton()
            sage: A.set_seed(T)
            sage: A.run()
            sage: A
            Core veering automaton with 2 vertices

            sage: A = ReducedCoreAutomaton()
            sage: A.set_seed(T)
            sage: A.run()
            sage: A
            Reduced core veering automaton with 1 vertex

        A more complicated surface in Q_1(1^2, -1^2)::

            sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
            sage: cols = 'BRBRBBBRBR'
            sage: T = VeeringTriangulation(fp, cols)

            sage: C = CoreAutomaton()
            sage: C.set_seed(T)
            sage: C.run(10)
            sage: C
            Partial core veering automaton with 11 vertices
            sage: C.run()
            sage: C
            Core veering automaton with 1074 vertices

            sage: C = ReducedCoreAutomaton()
            sage: C.set_seed(T)
            sage: C.run(10)
            sage: C
            Partial reduced core veering automaton with 11 vertices
            sage: C.run()
            sage: C
            Reduced core veering automaton with 356 vertices

        TESTS::

            sage: from veerer import VeeringTriangulation, CoreAutomaton
            sage: T = VeeringTriangulation("(0,2,~1)(1,~0,~2)", "RBB")
            sage: A = CoreAutomaton()
            sage: A
            Uninitialized core veering automaton with 0 vertex
            sage: A.set_seed(T)
            sage: A
            Partial core veering automaton with 1 vertex
            sage: A.run(1)
            sage: A
            Partial core veering automaton with 2 vertices
            sage: A.run(1)
            sage: A
            Core veering automaton with 2 vertices

            sage: A0 = CoreAutomaton()
            sage: A0.set_seed(T)
            sage: A0.run()
            sage: A0._flip_branch
            [[]]
            sage: A1 = CoreAutomaton()
            sage: A1.set_seed(T)
            sage: A1.run(1)
            sage: A1.run(1)
            sage: A1._flip_branch
            [[]]
        """
        if not self._flip_branch:
            raise ValueError('uninitialized automaton; call set_seed() first')

        if len(self._flip_branch) == 1 and not self._flip_branch[0]:
            return

        T = self._state
        state_branch = self._state_branch
        flip_branch = self._flip_branch
        flips = self._flips
        relabellings = self._relabellings
        graph = self._graph
        recol = None

        count = 0
        old_size = len(graph)
        while max_size is None or count < max_size:
            assert T == state_branch[-1]
            immutable_state = state_branch[-1]
            new_flip = flip_branch[-1].pop()

            # TODO: check to be removed
            if CHECK:
                T._check()

            if self._verbosity >= 2:
                print('[automaton] NEW LOOP')
                print('[automaton] at %s' % T)
                print('[automaton] going along %s' % new_edge)
            elif self._verbosity:
                if len(graph) > old_size + 500:
                    print('\r[automaton] %8d      %8d      %.3f      ' % (len(graph), len(branch), time() - t0), end='')
                    sys.stdout.flush()
                    old_size = len(graph)

            # perform a flip
            is_target_valid, flip_back_data = self._flip(new_flip)

            if self._verbosity >= 2:
                print('[automaton] ... landed at %s' % T)
                sys.stdout.flush()

            if is_target_valid: # = T is a valid vertex
                r, _ = T.best_relabelling()
                rinv = perm_invert(r)
                T.relabel(r, check=False)
                new_state = T.copy(mutable=False)
                new_state.set_immutable()
                graph[state_branch[-1]].append((new_state, new_flip))

                if self._verbosity >= 2:
                    print('[automaton] valid new state')
                    sys.stdout.flush()

                if new_state not in graph:
                    # new target vertex
                    if self._verbosity >= 2:
                        print('[automaton] new vertex')
                        sys.stdout.flush()

                    # new core
                    flips.append(flip_back_data)
                    relabellings.append(rinv)
                    graph[new_state] = []
                    state_branch.append(new_state)

                    # TODO: one can optimize the computation of forward flippable edges knowing
                    # how the current state was built
                    flip_branch.append(self._forward_flips(T))
                    count += 1
                    continue

                else:
                    # existing target vertex
                    if self._verbosity >= 2:
                        print('[automaton] known vertex')
                        sys.stdout.flush()
                    T.relabel(rinv, check=False)

            else:  # = T is not a valid vertex
                if self._verbosity >= 2:
                    print('[automaton] invalid vertex')
                    sys.stdout.flush()

            # not core or already visited
            self._flip_back(flip_back_data)

            # TODO: check to be removed
            if CHECK:
                T._check()
            assert T == state_branch[-1], (T, state_branch[-1])

            while flips and not flip_branch[-1]:
                rinv = relabellings.pop()
                T.relabel(rinv, check=False)
                flip_back_data = flips.pop()
                self._flip_back(flip_back_data)
                state_branch.pop()
                flip_branch.pop()
                assert T == state_branch[-1]
            if not flip_branch[-1]:
                # done
                assert len(flip_branch) == 1
                break

        if self._verbosity == 1:
            print('\r[automaton] %8d      %8d      %.3f      ' % (len(graph), len(flip_branch), time() - t0))
            sys.stdout.flush()

        assert self._state is T


class FlipGraph(Automaton):
    r"""
    The flip graph

    EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation([(0,1,2),(-1,-2,-3)])
            sage: FlipGraph(T)
            Triangulation automaton with 1 vertex

            sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
            sage: T = Triangulation(fp)
            sage: FlipGraph(T)
            Triangulation automaton with 236 vertices
    """
    _name = 'triangulation'

    def _seed_setup(self, state):
        if not isinstance(state, Triangulation):
            raise ValueError('invalid seed')
        state = state.copy(mutable=True)
        state.set_canonical_labels()
        return state

    def _forward_flips(self, state):
        r"""
        Return the list of forward flippable edges from ``state``
        """
        return state.flippable_edges()

    def _flip(self, flip_data):
        e = flip_data
        self._state.flip(e, check=CHECK)
        return True, e

    def _flip_back(self, flip_back_data):
        e = flip_back_data
        self._state.flip_back(e, check=CHECK)


class CoreAutomaton(Automaton):
    r"""
    Automaton of core veering triangulations.
    """
    _name = 'core veering'

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulation) or not state.is_core():
            raise ValueError('invalid seed')
        state = state.copy(mutable=True)
        state.set_canonical_labels()
        return state

    def _forward_flips(self, state):
        r"""
        Return the list of forward flippable edges from ``state``
        """
        ffe = state.forward_flippable_edges()
        return [(x, col) for x in ffe for col in (BLUE, RED)]

    def _flip(self, flip_data):
        e, col = flip_data
        old_col = self._state.colour(e)
        self._state.flip(e, col, check=CHECK)
        flip_back_data = (e, old_col)
        return self._state.edge_has_curve(e), flip_back_data

    def _flip_back(self, flip_back_data):
        e, old_col = flip_back_data
        self._state.flip_back(e, old_col, check=CHECK)


class ReducedCoreAutomaton(Automaton):
    _name = 'reduced core veering'

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulation) or not state.is_core():
            raise ValueError('invalid seed')

        state = state.copy(mutable=True)
        state.forgot_forward_flippable_colour()
        state.set_canonical_labels()
        return state

    def _forward_flips(self, state):
        r"""
        Return the list of forward flippable edges from ``state``
        """
        ffe = state.purple_edges()
        # TODO: check to be removed
        if CHECK:
            assert ffe == state.forward_flippable_edges()
        return [(x, col) for x in ffe for col in (BLUE, RED)]

    def _flip(self, flip_data):
        T = self._state
        e, col = flip_data

        # TODO: check to be removed
        if CHECK:
            assert T.is_forward_flippable(e)
        old_col = self._state.colour(e)
        self._state.flip(e, col, reduced=False, check=CHECK)

        if not self._state.edge_has_curve(e):
            return False, (e, old_col, ())

        recolorings = []
        # only two edges are succeptible to become forward flippable
        a, b, c, d = T.square_about_edge(e)

        # assertions to be removed
        assert T._colouring[a] == RED, (a, colour_to_string(T._colouring[a]))
        assert T._colouring[b] == BLUE, (b, colour_to_string(T._colouring[b]))
        assert T._colouring[c] == RED, (c, colour_to_string(T._colouring[c]))
        assert T._colouring[d] == BLUE, (d, colour_to_string(T._colouring[d]))
        if col == BLUE:
            if T.is_forward_flippable(b):
                recolorings.append((b, BLUE))
                T._colouring[b] = PURPLE
                T._colouring[T._ep[b]] = PURPLE
            if b != d and T.is_forward_flippable(d):
                recolorings.append((d, BLUE))
                T._colouring[d] = PURPLE
                T._colouring[T._ep[d]] = PURPLE

            # assertions to be removed
            if CHECK:
                assert not T.is_forward_flippable(e)
                assert not T.is_forward_flippable(a)
                assert not T.is_forward_flippable(c)

        elif col == RED:
            if T.is_forward_flippable(a):
                recolorings.append((a, RED))
                T._colouring[a] = PURPLE
                T._colouring[T._ep[a]] = PURPLE
            if a != c and T.is_forward_flippable(c):
                recolorings.append((c, RED))
                T._colouring[c] = PURPLE
                T._colouring[T._ep[c]] = PURPLE

            # assertions to be removed
            if CHECK:
                assert not T.is_forward_flippable(e)
                assert not T.is_forward_flippable(b)
                assert not T.is_forward_flippable(d)

        return True, (e, old_col, recolorings)

    def _flip_back(self, flip_back_data):
        e, old_col, recolorings = flip_back_data
        for ee, ccol in recolorings:
            self._state._colouring[ee] = ccol
            self._state._colouring[self._state._ep[ee]] = ccol
        self._state.flip_back(e, old_col, check=CHECK)


class GeometricAutomaton(Automaton):
    r"""
    Automaton of geometric veering triangulations.

    A veering triangulation is called geometric, if the set of
    associated L^oo-vector data is full dimensional in the
    ambient stratum.

    EXAMPLES::

        sage: from veerer import *
        sage: fp = "(0,~7,6)(1,~8,~2)(2,~6,~3)(3,5,~4)(4,8,~5)(7,~1,~0)"
        sage: cols = "RBRBRBBBB"
        sage: vt = VeeringTriangulation(fp, cols)
        sage: GeometricAutomaton(vt)
        Geometric veering automaton with 54 vertices

    One can check that the cardinality is indeed correct::

        sage: sum(x.is_geometric() for x in CoreAutomaton(vt))
        54
    """
    _name = 'geometric veering'

    def _setup(self, backend=None):
        self._backend = backend

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulation) or not state.is_geometric(backend=self._backend):
            raise ValueError('invalid seed')
        if self._backend is None:
            from .polyhedron.cone import default_backend
            self._backend = default_backend(state.base_ring())
        state = state.copy(mutable=True)
        state.set_canonical_labels()
        return state

    def _forward_flips(self, state):
        r"""
        Return the list of forward flippable edges from ``state``
        """
        return state.geometric_flips(backend=self._backend)

    def _flip(self, flip_data):
        edges, col = flip_data
        assert all(self._state.colour(e) == self._state.colour(edges[0]) for e in edges)
        flip_back_data = (edges, self._state.colour(edges[0]))
        for e in edges:
            self._state.flip(e, col, check=CHECK)
        if CHECK and not self._state.is_geometric(backend=self._backend):
            raise RuntimeError('that was indeed possible!')
        return True, flip_back_data

    def _flip_back(self, flip_back_data):
        edges, old_col = flip_back_data
        for e in edges:
            self._state.flip_back(e, old_col, check=CHECK)


class GeometricAutomatonSubspace(Automaton):
    r"""
    Automaton of geometric veering triangulations with a linear subspace constraint.

    A veering triangulation is called geometric, if the set of
    associated L^oo-vector data is full dimensional in the
    ambient stratum.

    EXAMPLES::

        sage: from veerer import *

        sage: vt, s, t = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
        sage: f = VeeringTriangulationLinearFamily(vt, [s, t])
        sage: GeometricAutomatonSubspace(f)
        Geometric veering linear constraint automaton with 6 vertices

    One can check that the cardinality is indeed correct::

        sage: sum(x.is_geometric() for x in CoreAutomaton(vt))
        54

    Some L-shape surfaces::

        sage: for n in range(3, 7):
        ....:     vt, s, t = VeeringTriangulations.L_shaped_surface(1, n-2, 1, 1)
        ....:     f = VeeringTriangulationLinearFamily(vt, [s, t])
        ....:     print('n={}: {}'.format(n, GeometricAutomatonSubspace(f)))
        n=3: Geometric veering linear constraint automaton with 6 vertices
        n=4: Geometric veering linear constraint automaton with 86 vertices
        n=5: Geometric veering linear constraint automaton with 276 vertices
        n=6: Geometric veering linear constraint automaton with 800 vertices
    """
    _name = 'geometric veering linear constraint'

    def _setup(self, backend=None):
        self._backend = backend

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulationLinearFamily) or not state.is_geometric(backend=self._backend):
            raise ValueError('invalid seed')
        if self._backend is None:
            from .polyhedron.cone import default_backend
            self._backend = default_backend(state.base_ring())
        state = state.copy(mutable=True)
        state.set_canonical_labels()
        return state

    def _forward_flips(self, state):
        r"""
        Return the list of forward flippable edges from ``state``
        """
        return state.geometric_flips(backend=self._backend)

    def _flip(self, flip_data):
        edges, col = flip_data
        assert all(self._state.colour(e) == self._state.colour(edges[0]) for e in edges)
        flip_back_data = (edges, self._state.colour(edges[0]))
        for e in edges:
            self._state.flip(e, col, check=CHECK)
        if CHECK:
            self._state._check(RuntimeError)
            if not self._state.is_geometric(backend=self._backend):
                raise RuntimeError
        return True, flip_back_data

    def _flip_back(self, flip_back_data):
        edges, old_col = flip_back_data
        for e in edges:
            self._state.flip_back(e, old_col, check=CHECK)
