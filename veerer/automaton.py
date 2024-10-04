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

        sage: from surface_dynamics import Stratum                       # optional - surface_dynamics
        sage: strata = [Stratum([2], 1), Stratum([2,-1,-1], 2),           # optional - surface_dynamics
        ....:           Stratum([2,2], 2), Stratum([1,1], 1)]
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

    def __init__(self, seed=None, verbosity=0, backward=False, **extra_kwds):
        # verbosity level
        self._verbosity = int(verbosity)

        # whether to explore backward flips
        self._backward = backward

        # we encode the graph in two dictionaries where the keys are the states
        # and the values are respectively the outgoing or incoming edges
        # In both cases, the neighbors are encoded by pairs (neighbor, flip_data)
        self._forward_neighbors = {}
        self._backward_neighbors = {}

        # list of seeds to be considered
        self._seeds = []

        # current state (mutable triangulations)
        self._state = None

        # current branches for running the DFS
        self._state_branch = []
        self._flip_branch = []

        # queue of states on which we still have run backward flips from
        # (each time an element is poped its backward flip neighbors will
        #  populate the _seeds list)
        self._backward_flip_queue = []

        # the flips and relabelling performed along the branch
        self._flips = []
        self._relabellings = []

        self._setup(**extra_kwds)

        if seed is not None:
            self.add_seed(seed)
            self.run()

    def __str__(self):
        status = None
        if self._backward_flip_queue or len(self._flip_branch) != 1 or self._flip_branch[-1]:
            status = 'partial '
        else:
            status = ''
        return ("%s%s automaton with %s %s" % (status, self._name, len(self._forward_neighbors), 'vertex' if len(self._forward_neighbors) <= 1 else 'vertices')).capitalize()

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self._forward_neighbors)

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
        return iter(self._forward_neighbors)

    def __contains__(self, item):
        if not isinstance(item, Triangulation):
            return False
        if item._mutable:
            item = item.copy(mutable=False)
        return item in self._forward_neighbors

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

        for g, neighb in self._forward_neighbors.items():
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
        for vt in self._forward_neighbors:
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
        return len(self._forward_neighbors)

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
        return sum(len(flips) for flips in self._forward_neighbors.values())

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
    def from_triangulation(cls, T, *args, **kwds):
        A = cls()
        A.add_seed(T)
        A.run(**kwds)
        return A

    # TODO: it is a bit absurd to have this in the generic class since this does not make
    # sense for triangulations
    # TODO: this is broken for GeometricAutomaton with QuadraticStratum(8)
    # TODO: one should be more careful whether one wants the poles to be half edges or vertices
    @classmethod
    def from_stratum(self, stratum, **kwds):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import Stratum         # optional - surface_dynamics
            sage: CoreAutomaton.from_stratum(Stratum([2], 1))  # optional - surface_dynamics
            Core veering automaton with 86 vertices
            sage: GeometricAutomaton.from_stratum(Stratum([2], 1))  # optional - surface_dynamics
            Geometric veering automaton with 54 vertices

            sage: Q = Stratum([8], 2)                             # optional - surface_dynamics
            sage: A = CoreAutomaton.from_stratum(Q, max_size=100) # optional - surface_dynamics
            sage: A                                               # optional - surface_dynamics
            Partial core veering automaton with 101 vertices
        """
        return self.from_triangulation(VeeringTriangulation.from_stratum(stratum), **kwds)

    ######################
    # DFS implementation #
    ######################

    def _setup(self):
        pass

    def _seed_setup(self, state):
        r"""
        Transform ``state`` into a canonical triangulation for the automaton.

        To be implemented in subclasses.
        """
        raise NotImplementedError

    def _forward_flips(self, state):
        r"""
        Return the list of forward flip_data from ``state``.

        To be implemented in subclasses.
        """
        raise NotImplementedError

    def _backward_flips(self, state):
        r"""
        Return the list of backward flips data from ``state``.

        By default it return an empty tuple (ie backward flips are not
        explored). Can be overriden in subclasses.
        """
        raise NotImplementedError

    def _flip(self, state, flip_data, status=False):
        r"""
        Perform the flips to ``state`` given by ``flip_data``.

        If ``status`` is ``False`` the method return ``flip_back_data``
        that can be used to flip back ``state`` to its previous form
        using :meth:`_flip_back`. If it is ``True`` the method return
        a pair ``(is_valid, flip_back_data)`` where ``flip_back_data``
        is as before and ``is_valid`` test whether the obtained triangulation
        is valid.

        To be implemented in subclasses
        """
        raise NotImplementedError

    def _flip_back(self, flip_back_data, status=False):
        r"""
        Perform back the flips on ``state`` given by ``flip_back_data``.

        To be implemented in subclasses
        """
        raise NotImplementedError

    def add_seed(self, state, setup=True):
        r"""
        Add the seed ``state`` to the search.

        Return ``0`` if the state is already present in the graph and ``1``
        otherwise.
        """
        if setup:
            state = self._seed_setup(state)

            # TODO: check to be removed
            if CHECK:
                state._check()
            state = state.copy(mutable=False)

        assert not state._mutable

        if state in self._forward_neighbors:
            if self._verbosity >= 2:
                print('[add_seed] state=%s already in the graph' % state)
            return 0

        if self._verbosity >= 2:
            print('[add_seed] adding state=%s to the list of seeds' % state)

        self._seeds.append(state)
        return 1

    def reset_current_state(self, state):
        """
        The method assumes that ``state`` is immutable and is not present in the graph.
        """
        if not state._mutable:
            immutable_state = state
            state = immutable_state.copy(mutable=True)
        else:
            immutable_state = state.copy(mutable=False)

        if immutable_state in self._forward_neighbors:
            raise ValueError('state already in the graph')

        self._forward_neighbors[immutable_state] = []
        self._backward_neighbors[immutable_state] = []
        self._state = state
        self._state_branch.clear()
        self._state_branch.append(immutable_state)
        self._flip_branch.clear()
        self._flip_branch.append(self._forward_flips(state))
        if self._backward:
            self._backward_flip_queue.append(immutable_state)

    def next_seed(self):
        r"""
        Return the next seed or ``None`` if there is none.

        This function modifies the `_seeds` and `_backward_flip_queue` attributes.
        """
        while self._seeds or self._backward_flip_queue:
            while self._seeds:
                seed = self._seeds.pop()
                assert not seed._mutable
                if seed not in self._forward_neighbors:
                    # not explored forward yet
                    return seed
                elif self._verbosity >= 2:
                    print('[next_seed] seed already in the graph')

            if not self._backward_flip_queue:
                return

            state = self._backward_flip_queue.pop().copy(mutable=True)
            for flip_data in self._backward_flips(state):
                is_target_valid, flip_back_data = self._flip_back(state, flip_data, True)
                if is_target_valid:
                    if self._verbosity >= 2:
                        print('[next_seed] explore backward flip')
                    r, _ = state.best_relabelling()
                    new_state = state.copy(mutable=True)
                    new_state.relabel(r, check=False)
                    new_state.set_immutable()
                    self.add_seed(new_state, setup=False)
                self._flip(state, flip_back_data, False)

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
            sage: A.add_seed(T)
            1
            sage: A.run()
            sage: A
            Core veering automaton with 2 vertices

            sage: A = ReducedCoreAutomaton()
            sage: A.add_seed(T)
            1
            sage: A.run()
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
            sage: C
            Partial core veering automaton with 10 vertices
            sage: C.run()
            sage: C
            Core veering automaton with 1074 vertices

            sage: C = ReducedCoreAutomaton()
            sage: C.add_seed(T)
            1
            sage: C.run(10)
            sage: C
            Partial reduced core veering automaton with 10 vertices
            sage: C.run()
            sage: C
            Reduced core veering automaton with 356 vertices

        TESTS::

            sage: from veerer import VeeringTriangulation, CoreAutomaton
            sage: T = VeeringTriangulation("(0,2,~1)(1,~0,~2)", "RBB")
            sage: A = CoreAutomaton()
            sage: A
            Partial core veering automaton with 0 vertex
            sage: A.add_seed(T)
            1
            sage: A.run(1)
            sage: A
            Partial core veering automaton with 1 vertex
            sage: A.run(1)
            sage: A
            Partial core veering automaton with 2 vertices
            sage: A.run(1)
            sage: A
            Core veering automaton with 2 vertices

            sage: A0 = CoreAutomaton()
            sage: A0.add_seed(T)
            1
            sage: A0.run()
            sage: A0._flip_branch
            [[]]
            sage: A1 = CoreAutomaton()
            sage: A1.add_seed(T)
            1
            sage: A1.run(1)
            sage: A1.run(1)
            sage: A1.run(1)
            sage: A1._flip_branch
            [[]]
        """
        backward_flip_queue = self._backward_flip_queue
        T = self._state
        state_branch = self._state_branch
        flip_branch = self._flip_branch
        flips = self._flips
        relabellings = self._relabellings
        fgraph = self._forward_neighbors
        bgraph = self._backward_neighbors
        recol = None

        count = 0
        old_size = len(fgraph)
        while max_size is None or count < max_size:
            if self._verbosity >= 2:
                print('[automaton] NEW LOOP')

            while not flip_branch or (len(flip_branch) == 1 and not flip_branch[0]):
                # forward search is over... find a new start
                seed = self.next_seed()
                if seed is None:
                    # no seed available anymore
                    if self._verbosity >= 2:
                        print('[automaton] done')
                    return
                if self._verbosity >= 2:
                    print('[automaton] seed %s' % seed)
                sys.stdout.flush()
                self.reset_current_state(seed)
                T = self._state
                count += 1
                if max_size is not None and count >= max_size:
                    return

            # next step of forward search
            assert T == state_branch[-1]
            immutable_state = state_branch[-1]
            assert not immutable_state._mutable
            new_flip = flip_branch[-1].pop()

            # TODO: check to be removed
            if CHECK:
                T._check()
            if self._verbosity >= 2:
                print('[automaton] at %s' % T)
                print('[automaton] exploring flip %s' % str(new_flip))
                sys.stdout.flush()
            elif self._verbosity:
                if len(fgraph) > old_size + 500:
                    print('\r[automaton] %8d      %8d      %.3f      ' % (len(fgraph), len(branch), time() - t0), end='')
                    sys.stdout.flush()
                    old_size = len(fgraph)

            # perform a flip on self._state
            is_target_valid, flip_back_data = self._flip(T, new_flip, True)

            if self._verbosity >= 2:
                print('[automaton] ... landed at %s' % T)
                sys.stdout.flush()

            if is_target_valid: # = T is a valid vertex
                r, _ = T.best_relabelling()
                rinv = perm_invert(r)
                T.relabel(r, check=False)
                new_state = T.copy(mutable=False)
                new_state.set_immutable()
                fgraph[state_branch[-1]].append((new_state, new_flip))

                if self._verbosity >= 2:
                    print('[automaton] valid new state')
                    sys.stdout.flush()

                if new_state not in fgraph:
                    # new target vertex
                    if self._verbosity >= 2:
                        print('[automaton] new vertex')
                        sys.stdout.flush()

                    fgraph[new_state] = []
                    bgraph[new_state] = [(state_branch[-1], new_flip)]

                    if self._backward:
                        if self._verbosity >= 2:
                            print('[automaton] adding it to backward_flip_queue')
                        self._backward_flip_queue.append(new_state)

                    # TODO: one can optimize the computation of forward flippable edges knowing
                    # how the current state was built
                    forward_flips = self._forward_flips(T)
                    if forward_flips:
                        flips.append(flip_back_data)
                        relabellings.append(rinv)
                        state_branch.append(new_state)
                        flip_branch.append(forward_flips)
                        count += 1
                        continue
                else:
                    bgraph[new_state].append((state_branch[-1], new_flip))

                # dead end or existing target vertex: backtrack the DFS search
                if self._verbosity >= 2:
                    print('[automaton] dead end or known state')
                    sys.stdout.flush()
                T.relabel(rinv, check=False)


            else:  # = T is not a valid vertex
                if self._verbosity >= 2:
                    print('[automaton] invalid vertex')
                    sys.stdout.flush()

            # not core or already visited
            self._flip_back(self._state, flip_back_data, False)

            # TODO: check to be removed
            if CHECK:
                T._check()
            assert T == state_branch[-1], (T, state_branch[-1])

            while flips and not flip_branch[-1]:
                rinv = relabellings.pop()
                T.relabel(rinv, check=False)
                flip_back_data = flips.pop()
                self._flip_back(self._state, flip_back_data, False)
                state_branch.pop()
                flip_branch.pop()
                assert T == state_branch[-1]

        if self._verbosity == 1:
            print('\r[automaton] %8d      %8d      %.3f      ' % (len(fgraph), len(flip_branch), time() - t0))
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

    def _flip(self, state, flip_data, status=False):
        e = flip_data
        state.flip(e, check=CHECK)
        return (True, e) if status else e

    def _flip_back(self, state, flip_back_data, status=False):
        e = flip_back_data
        state.flip_back(e, check=CHECK)
        return (True, e) if status else e


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

    def _flip(self, state, flip_data, status=False):
        e, col = flip_data
        old_col = state.colour(e)
        state.flip(e, col, check=CHECK)
        flip_back_data = (e, old_col)
        return (state.edge_has_curve(e), flip_back_data) if status else flip_back_data

    def _flip_back(self, state, flip_back_data, status=False):
        e, old_col = flip_back_data
        col = state.colour(e)
        state.flip_back(e, old_col, check=CHECK)
        flip_data = (e, col)
        return (state.edge_has_curve(e), flip_data) if status else flip_data

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

    def _flip(self, state, flip_data, status=False):
        e, col = flip_data

        # TODO: check to be removed
        if CHECK:
            assert state.is_forward_flippable(e)
        old_col = state.colour(e)
        state.flip(e, col, reduced=False, check=CHECK)

        if not self._state.edge_has_curve(e):
            return False, (e, old_col, ())

        recolorings = []
        # only two edges are succeptible to become forward flippable
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

            # assertions to be removed
            if CHECK:
                assert not state.is_forward_flippable(e)
                assert not state.is_forward_flippable(b)
                assert not state.is_forward_flippable(d)

        return True, (e, old_col, recolorings)

    def _flip_back(self, state, flip_back_data, status=False):
        e, old_col, recolorings = flip_back_data
        for ee, ccol in recolorings:
            state._colouring[ee] = ccol
            state._colouring[self._state._ep[ee]] = ccol
        state.flip_back(e, old_col, check=CHECK)


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

        sage: GeometricAutomaton(vt0)
        Geometric veering automaton with 3 vertices
        sage: GeometricAutomaton(vt1)
        Geometric veering automaton with 2 vertices
        sage: GeometricAutomaton(vt2)
        Geometric veering automaton with 1 vertex

        sage: GeometricAutomaton(vt0, backward=True)
        Geometric veering automaton with 3 vertices
        sage: GeometricAutomaton(vt1, backward=True)
        Geometric veering automaton with 3 vertices
        sage: GeometricAutomaton(vt2, backward=True)
        Geometric veering automaton with 3 vertices
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
        Return the list of Delaunay forward flippable edges from ``state``.
        """
        return state.delaunay_flips(backend=self._backend)

    def _backward_flips(self, state):
        """
        Return the list of Delaunay backward flippable edges from ``state``.
        """
        return state.backward_delaunay_flips(backend=self._backend)

    def _flip(self, state, flip_data, status=False):
        edges, col = flip_data
        assert all(state.colour(e) == state.colour(edges[0]) for e in edges)
        flip_back_data = (edges, state.colour(edges[0]))
        for e in edges:
            state.flip(e, col, check=CHECK)
        if CHECK:
            state._check(RuntimeError)
            if not state.is_geometric(backend=self._backend):
                raise RuntimeError
        return (True, flip_back_data) if status else flip_back_data

    def _flip_back(self, state, flip_back_data, status=False):
        edges, old_col = flip_back_data
        assert all(state.colour(e) == state.colour(edges[0]) for e in edges)
        flip_data = (edges, state.colour(edges[0]))
        for e in edges:
            state.flip_back(e, old_col, check=CHECK)
        if CHECK:
            state._check(RuntimeError)
            if not state.is_geometric(backend=self._backend):
                raise RuntimeError
        return (True, flip_data) if status else flip_data

    def sources(self):
        r"""
        Iterate through sources (states with no backward neighbor) in this automaton.
        """

    def sinks(self):
        r"""
        Iterate through sinks (states with no forward neighbor) in this automaton.
        """
        return (x for x, forward_neighbors in self._forward_neighbors.items() if not forward_neighbors)

    def cylinder_diagrams(self, col=RED):
        r"""
        A cylinder diagram is a cylindrical veering triangulation up to twist action.

        EXAMPLES::

            sage: from veerer import RED, BLUE, GeometricAutomaton
            sage: from surface_dynamics import Stratum                       # optional - surface_dynamics
            sage: A = GeometricAutomaton.from_stratum(Stratum([2], 1))       # optional - surface_dynamics
            sage: len(A.cylinder_diagrams(RED))                              # optional - surface_dynamics
            2
            sage: A = GeometricAutomaton.from_stratum(Stratum([1, 1], 1))    # optional - surface_dynamics
            sage: len(A.cylinder_diagrams(RED))                              # optional - surface_dynamics
            4

        For Veech surfaces, the cylinder diagrams are in bijection with cusps::

            sage: from veerer import VeeringTriangulations, VeeringTriangulationLinearFamily
            sage: t, x, y = VeeringTriangulations.L_shaped_surface(1, 1, 1, 1)
            sage: f = VeeringTriangulationLinearFamily(t, [x, y])
            sage: A = f.geometric_automaton()
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
            ....:         A = F.geometric_automaton()
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


class GeometricAutomatonSubspace(GeometricAutomaton):
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

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulationLinearFamily) or not state.is_geometric(backend=self._backend):
            raise ValueError('invalid seed')
        if self._backend is None:
            from .polyhedron.cone import default_backend
            self._backend = default_backend(state.base_ring())
        state = state.copy(mutable=True)
        state.set_canonical_labels()
        return state
