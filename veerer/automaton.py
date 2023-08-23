r"""
Train-track and triangulations automata.
"""

import os
import sys
from time import time

from .veering_triangulation import VeeringTriangulation
from .constants import RED, BLUE, PURPLE, PROPERTIES_COLOURS, colour_to_char, colour_to_string

from . import env
from .env import sage, require_package


class Automaton(object):
    r"""
    Automaton of veering triangulations.

    EXAMPLES::

        sage: from veerer import *

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
        ....:     vt = VeeringTriangulation.from_stratum(stratum)
        ....:     print(CoreAutomaton(vt))
        ....:     print(ReducedCoreAutomaton(vt))
        H_2(2)
        Core veering automaton with 86 vertices
        Core veering automaton with 28 vertices
        Q_1(2, -1^2)
        Core veering automaton with 160 vertices
        Core veering automaton with 68 vertices
        Q_2(2^2)
        Core veering automaton with 846 vertices
        Core veering automaton with 305 vertices
        H_2(1^2)
        Core veering automaton with 876 vertices
        Core veering automaton with 234 vertices
    """
    _name = ''

    def __init__(self, seed=None, verbosity=0):
        # verbosity level
        self._verbosity = 0

        # we encode the graph in a dictionary where the keys are the states
        # (serialized with self._serialize) and the outgoing edges from a
        # given state are encoded in the corresponding value by a list of
        # pairs (end_state, flip_data)
        self._graph = {}

        # current state
        self._state = None
        self._serialized_state = None

        # current branch for running the DFS
        self._branch = []
        self._serialized_branch = []
        self._flips = []

        if seed is not None:
            self.set_seed(seed)
            self.run()

    def __str__(self):
        status = None
        if not self._branch:
            status = 'uninitialized '
        elif len(self._branch) != 1 or self._branch[-1]:
            status = 'partial '
        else:
            status = ''
        return ("%s%s veering automaton with %s %s" % (status, self._name, len(self._graph), 'vertex' if len(self._graph) <= 1 else 'vertices')).capitalize()

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self._graph)

    def one_triangulation(self):
        return next(iter(self))

    def __iter__(self):
        for s in sorted(self._graph):
            yield VeeringTriangulation.from_string(s)

    def __contains__(self, item):
        return item in self._graph

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
        require_package('sage', 'to_graph')
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
        for a in self._graph:
            T = VeeringTriangulation.from_string(a)
            T.rotate()
            aut[a] = T.iso_sig()
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
        for a in self._graph:
            T = VeeringTriangulation.from_string(a)
            T.conjugate()
            aut[a] = T.iso_sig()
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
        for g in self._graph:
            T = VeeringTriangulation.from_string(g)
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

            for gg, flip_data in self._graph[g]:
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
            sage: A.triangulations()
            <generator object ...>
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

    def geometric_triangulations(self, method=None):
        r"""
        Return an iterator over the pairs (veering triangulation,
        geometric polytope) when the triangulation is geometric.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
            sage: vt, P = next(A.geometric_triangulations())
            sage: vt.geometric_polytope() == P
            True
            sage: P.affine_dimension() == 8
            True

            sage: sum(1 for _ in A.geometric_triangulations())
            54
        """
        dim = self.one_triangulation().stratum_dimension()
        for vt in self:
            p = vt.geometric_polytope()
            if p.affine_dimension() == 2 * dim:
                yield vt, p

    def num_geometric_triangulations(self):
        r"""
        Return the number of geometric triangulations (among the states).

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
            sage: A.num_geometric_triangulations()
            54
        """
        return sum(vt.is_geometric() for vt in self)

    def cylindrical_triangulations(self):
        r"""
        Return an iterator over cylindrical configurations.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)

            sage: vt = next(A.cylindrical_triangulations())
            sage: sum(len(u) for u,_,_,_ in vt.cylinders(BLUE)) == 6 or \
            ....: sum(len(u) for u,_,_,_ in vt.cylinders(RED)) == 6
            True

            sage: sum(1 for _ in A.cylindrical_triangulations())
            24
        """
        for vt in self:
            if vt.is_cylindrical():
                yield vt

    def num_cylindrical_triangulations(self):
        r"""
        Return the number of cylindrical triangulations (among the states).

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton(T)
            sage: A.num_cylindrical_triangulations()
            24
        """
        return sum(vt.is_cylindrical() for vt in self)

    def set_seed(self, state):
        if self._state is not None:
            raise ValueError('seed already set')

        state = self._seed_setup(state)

        # TODO: check to be removed
        if env.CHECK:
            state._check()

        if self._verbosity:
            print('[automaton] stratum: %s' % state.stratum())
            print('[automaton] stratum dimension: %d' % state.stratum_dimension())
        if self._verbosity == 1:
            print('[automaton] size(graph)   size(path)    time')

        self._state = state
        self._serialized_state = iso_sig = state.iso_sig()

        self._serialized_branch = [iso_sig]

        self._graph[iso_sig] = []
        self._branch.append(self._forward_flips(self._state))

    @classmethod
    def from_triangulation(self, T, reduced=False, max_size=None, verbosity=0):
        if reduced:
            A = ReducedCoreAutomaton(verbosity=verbosity)
        else:
            A = CoreAutomaton(verbosity=verbosity)
        A.set_seed(T)
        A.run(max_size, verbosity)
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
            Core veering automaton with 100 vertices

        """
        return self.from_triangulation(VeeringTriangulation.from_stratum(stratum), reduced=reduced, **kwds)

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

    def _serialize(self, state):
        raise NotImplementedError

    def _unserialize(self, string):
        raise NotImplementedError

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
        """
        if not self._branch:
            return

        T = self._state
        iso_sig = self._serialized_state
        iso_sigs = self._serialized_branch
        branch = self._branch
        flips = self._flips
        graph = self._graph
        recol = None

        count = 0
        old_size = len(graph)
        while max_size is None or count < max_size:
            new_flip = branch[-1].pop()

            # TODO: check to be removed
            if env.CHECK:
                T._check()

            if self._verbosity >= 2:
                print('[automaton] NEW LOOP')
                print('[automaton] at %s with iso_sig %s' % (T.to_string(), iso_sig))
                print('[automaton] going along %s' % new_edge)
            elif self._verbosity:
                if len(graph) > old_size + 500:
                    print('\r[automaton] %8d      %8d      %.3f      ' % (len(graph), len(branch), time() - t0), end='')
                    sys.stdout.flush()
                    old_size = len(graph)

            # perform a flip
            is_target_valid, flip_back_data = self._flip(new_flip)

            if self._verbosity >= 2:
                print('[automaton] ... landed at %s' % T.to_string())
                sys.stdout.flush()

            if is_target_valid: # = T is a valid vertex
                new_iso_sig = self._serialize(T)
                graph[iso_sig].append((new_iso_sig, new_flip))

                if self._verbosity >= 2:
                    print('[automaton] it is core with iso_sig %s' % new_iso_sig)
                    sys.stdout.flush()

                if new_iso_sig not in graph:
                    # new target vertex
                    if self._verbosity >= 2:
                        print('[automaton] new vertex')
                        sys.stdout.flush()

                    # new core
                    flips.append(flip_back_data)
                    graph[new_iso_sig] = []

                    iso_sig = new_iso_sig
                    iso_sigs.append(iso_sig)

                    # TODO: one can optimize the computation of forward flippable edges knowing
                    # how the current state was built
                    branch.append(self._forward_flips(T))
                    count += 1
                    continue

                else:
                    # existing target vertex
                    if self._verbosity >= 2:
                        print('[automaton] known vertex')
                        sys.stdout.flush()

            else:  # = T is not a valid vertex
                if self._verbosity >= 2:
                    print('[automaton] invalid vertex')
                    sys.stdout.flush()

            # not core or already visited
            self._flip_back(flip_back_data)

            # TODO: check to be removed
            if env.CHECK:
                T._check()
            assert self._serialize(T) == iso_sigs[-1], (T.iso_sig(), iso_sigs[-1])

            while flips and not branch[-1]:
                flip_back_data = flips.pop()
                self._flip_back(flip_back_data)
                branch.pop()
                iso_sigs.pop()
                assert self._serialize(T) == iso_sigs[-1]
            if not branch[-1]:
                break
            iso_sig = iso_sigs[-1]

        if self._verbosity == 1:
            print('\r[automaton] %8d      %8d      %.3f      ' % (len(graph), len(branch), time() - t0))
            sys.stdout.flush()

        assert self._state is T
        self._serialized_state = iso_sig


class CoreAutomaton(Automaton):
    r"""
    Automaton of core veering triangulations.
    """
    _name = 'core'

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulation) or not state.is_core():
            raise ValueError('invalid seed')
        return state.copy()

    def _serialize(self, state):
        return state.iso_sig()

    def _unserialize(self):
        return VeeringTriangulation.from_string(string)

    def _forward_flips(self, state):
        r"""
        Return the list of forward flippable edges from ``state``
        """
        ffe = state.forward_flippable_edges()
        return [(x, col) for x in ffe for col in (BLUE, RED)]

    def _flip(self, flip_data):
        e, col = flip_data
        old_col = self._state.colour(e)
        self._state.flip(e, col)
        flip_back_data = (e, old_col)
        return self._state.edge_has_curve(e), flip_back_data

    def _flip_back(self, flip_back_data):
        e, old_col = flip_back_data
        self._state.flip_back(e, old_col)


class ReducedCoreAutomaton(Automaton):
    _name = 'reduced core'

    def _serialize(self, state):
        return state.iso_sig()

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulation) or not state.is_core():
            raise ValueError('invalid seed')

        state = state.copy() 
        state.forgot_forward_flippable_colour()
        return state

    def _forward_flips(self, state):
        r"""
        Return the list of forward flippable edges from ``state``
        """
        ffe = state.purple_edges()
        # TODO: check to be removed
        if env.CHECK:
            assert ffe == state.forward_flippable_edges()
        return [(x, col) for x in ffe for col in (BLUE, RED)]

    def _flip(self, flip_data):
        T = self._state
        e, col = flip_data

        # TODO: check to be removed
        if env.CHECK:
            assert T.is_forward_flippable(e)
        old_col = self._state.colour(e)
        self._state.flip(e, col, reduced=False)

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
            if env.CHECK:
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
            if env.CHECK:
                assert not T.is_forward_flippable(e)
                assert not T.is_forward_flippable(b)
                assert not T.is_forward_flippable(d)

        return True, (e, old_col, recolorings)

    def _flip_back(self, flip_back_data):
        e, old_col, recolorings = flip_back_data
        for ee, ccol in recolorings:
            self._state._colouring[ee] = ccol
            self._state._colouring[self._state._ep[ee]] = ccol
        self._state.flip_back(e, old_col)


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
    _name = 'geometric'

    def _seed_setup(self, state):
        if not isinstance(state, VeeringTriangulation) or not state.is_geometric():
            raise ValueError('invalid seed')
        return state.copy()

    def _serialize(self, state):
        return state.iso_sig()

    def _unserialize(self):
        return VeeringTriangulation.from_string(string)

    def _forward_flips(self, state):
        r"""
        Return the list of forward flippable edges from ``state``
        """
        return state.geometric_flips()

    def _flip(self, flip_data):
        flip_back_data = tuple((e, self._state.colour(e)) for e, _ in flip_data)
        for e, col in flip_data:
            self._state.flip(e, col)
        return True, flip_back_data

    def _flip_back(self, flip_back_data):
        for e, old_col in flip_back_data:
            self._state.flip_back(e, old_col)

class GeometricAutomatonWithLinearConstraint(Automaton):
    r"""
    Automaton of core veering triangulations with a linear constraint.

    This class can be used to certify linear invariant suborbifolds.

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
    _name = 'geometric'

    def _seed_setup(self, state):
        if not isinstance(state, (tuple, list)):
            raise TypeError('seed must be a pair (triangulation, linear constraint)')
        vt, L = state
        if not isinstance(vt, VeeringTriangulation) or not vt.is_geometric():
            raise ValueError('invalid triangulation')

    def _serialize(self, state):
        return state.iso_sig()

    def _unserialize(self):
        return VeeringTriangulation.from_string(string)

    def _forward_flips(self, state):
        r"""
        Return the list of forward flippable edges from ``state``
        """
        return state.geometric_flips()

    def _flip(self, flip_data):
        flip_back_data = tuple((e, self._state.colour(e)) for e, _ in flip_data)
        for e, col in flip_data:
            self._state.flip(e, col)
        return True, flip_back_data

    def _flip_back(self, flip_back_data):
        for e, old_col in flip_back_data:
            self._state.flip_back(e, old_col)


