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


class CoreAutomaton(object):
    r"""
    Automaton of core veering triangulations.

    EXAMPLES::

        sage: from veerer import *

        sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
        sage: A = CoreAutomaton.from_triangulation(T)
        sage: len(A)
        2

        sage: T = VeeringTriangulation('(0,1,2)', 'BBR')
        sage: A = CoreAutomaton.from_triangulation(T)
        sage: len(A)
        2

    The examples below is Q(1^2, -1^2) with two folded edges::

        sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
        sage: cols = 'BRBRBBBRBR'
        sage: T = VeeringTriangulation(fp, cols)
        sage: CoreAutomaton.from_triangulation(T)
        Core veering automaton with 1074 vertices
        sage: CoreAutomaton.from_triangulation(T, reduced=True)
        Core veering automaton with 356 vertices

    Exploring strata::

        sage: from surface_dynamics import *                             # optional - surface_dynamics
        sage: strata = [AbelianStratum(2), QuadraticStratum(2,-1,-1),    # optional - surface_dynamics
        ....:           QuadraticStratum(2,2), AbelianStratum(1,1)]
        sage: for stratum in strata:                                     # optional - surface_dynamics
        ....:     print(stratum)
        ....:     print(CoreAutomaton.from_stratum(stratum))
        ....:     print(CoreAutomaton.from_stratum(stratum, reduced=True))
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
    def __init__(self, graph):
        self._graph = graph
        self._iso_sigs = sorted(graph)
        self._index = {sig: index for index, sig in enumerate(self._iso_sigs)}

    def __str__(self):
        return "Core veering automaton with %s vertices" % len(self._graph)

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self._iso_sigs)

    def one_triangulation(self):
        return next(iter(self))

    def __iter__(self):
        for s in self._iso_sigs:
            yield VeeringTriangulation.from_string(s)

    def __contains__(self, item):
        return item in self._iso_sigs

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
            sage: A = CoreAutomaton.from_triangulation(T)
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
            for gg, e, old_col, col in neighb:
                G.add_edge(g, gg, (e, old_col, col))

        return G

    def rotation_automorphism(self):
        r"""
        Return the automorphism of the vertices that corresponds to rotation.

        Note that this automorphism reverses the edge direction.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton.from_triangulation(T)
            sage: rot = A.rotation_automorphism()

            sage: G = A.to_graph()
            sage: all(G.has_edge(rot[b], rot[a]) for a,b in G.edges(labels=False))
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
            sage: A = CoreAutomaton.from_triangulation(T)
            sage: conj = A.conjugation_automorphism()

            sage: G = A.to_graph()
            sage: all(G.has_edge(conj[a], conj[b]) for a,b in G.edges(labels=False))
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
            sage: A = CoreAutomaton.from_triangulation(T)
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
        for g in self._iso_sigs:
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

            for gg, e, old_col, new_col in self._graph[g]:
                f.write('    %s -> %s [color="%s;%f:%s"];\n' % (g, gg, old_col, 0.5, new_col))
        f.write('}\n')

        if filename is not None:
            f.close()

    def triangulations(self):
        r"""
        Return an iterator over the veering triangulations in this automaton.

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,6,~5)(1,8,~7)(2,7,~6)(3,~1,~8)(4,~2,~3)(5,~0,~4)", "RRRBBBBBB")
            sage: A = CoreAutomaton.from_triangulation(T)
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
            sage: A = CoreAutomaton.from_triangulation(T)
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
            sage: A = CoreAutomaton.from_triangulation(T)
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
            sage: A = CoreAutomaton.from_triangulation(T)
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
            sage: A = CoreAutomaton.from_triangulation(T)
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
            sage: A = CoreAutomaton.from_triangulation(T)
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
            sage: A = CoreAutomaton.from_triangulation(T)
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
            sage: A = CoreAutomaton.from_triangulation(T)

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
            sage: A = CoreAutomaton.from_triangulation(T)
            sage: A.num_cylindrical_triangulations()
            24
        """
        return sum(vt.is_cylindrical() for vt in self)

    @classmethod
    def from_triangulation(self, T, reduced=False, max_size=None, verbosity=0):
        r"""
        Build the core automaton of ``T``.

        INPUT:

        - ``T`` - a veering triangulation (no GREEN allowed)

        - ``reduced`` - boolean - whether we forgot colours of forward flippable edges

        - ``verbosity`` - integer (default ``0``) - the verbosity level. If nonzero print
          information during the computation.
        """
        # TODO: if we have boundaries we also need to go backward!
        assert T.is_core()
        if reduced:
            T.forgot_forward_flippable_colour()

            # TODO: check to be removed
            if env.CHECK:
                T._check()

        if verbosity:
            print('[automaton] stratum: %s' % T.stratum())
            print('[automaton] stratum dimension: %d' % T.stratum_dimension())
        if verbosity == 1:
            print('[automaton] size(graph)   size(path)    time')

        T = T.copy()
        iso_sig = T.iso_sig()
        d = T.stratum_dimension()

        graph = {iso_sig: []}
        iso_sigs = [iso_sig]  # iso_sig sequence
        flips = []       # backward flip sequence and recolouring
        branch = []    # forward flips available

        if reduced:
            ffe = T.purple_edges()
            assert ffe == T.forward_flippable_edges()
        else:
            ffe = T.forward_flippable_edges()

        branch.append([])
        branch[-1].extend((x, BLUE) for x in ffe)
        branch[-1].extend((x, RED) for x in ffe)
        e, col = branch[-1].pop()
        recol = None
        old_size = 0
        t0 = time()
        while max_size is None or len(graph) < max_size:
            # TODO: check to be removed
            if env.CHECK:
                T._check()

            if verbosity >= 2:
                print('[automaton] NEW LOOP')
                print('[automaton] at %s with iso_sig %s' % (T.to_string(), iso_sig))
                print('[automaton] going along (%s, %s)' % (e, col))
            elif verbosity:
                if len(graph) > old_size + 500:
                    print('\r[automaton] %8d      %8d      %.3f      ' % (len(graph), len(branch), time() - t0), end='')
                    sys.stdout.flush()
                    old_size = len(graph)

            # some safety check to be disabled
            # assert len(flips) + 1 == len(iso_sigs) == len(branch)
            # assert iso_sig == iso_sigs[-1]
            # assert T.iso_sig() == iso_sig

            # explore the automaton in DFS
            old_col = T.colour(e)

            # assertion to be removed
            assert not reduced or old_col == PURPLE
            T.flip(e, col, reduced=False)

            if verbosity >= 2:
                print('[automaton] ... landed at %s' % T.to_string())
                sys.stdout.flush()

            if T.edge_has_curve(e):  # = T is core
                if reduced:
                    recol = []
                    # only two edges are succeptible to become forward flippable
                    a, b, c, d = T.square_about_edge(e)

                    # assertions to be removed
                    assert T._colouring[a] == RED, (a, colour_to_string(T._colouring[a]))
                    assert T._colouring[b] == BLUE, (b, colour_to_string(T._colouring[b]))
                    assert T._colouring[c] == RED, (c, colour_to_string(T._colouring[c]))
                    assert T._colouring[d] == BLUE, (d, colour_to_string(T._colouring[d]))
                    if col == BLUE:
                        if T.is_forward_flippable(b):
                            recol.append((b, BLUE))
                            T._colouring[b] = PURPLE
                            T._colouring[T._ep[b]] = PURPLE
                        if b != d and T.is_forward_flippable(d):
                            recol.append((d, BLUE))
                            T._colouring[d] = PURPLE
                            T._colouring[T._ep[d]] = PURPLE

                        # assertions to be removed
                        assert not T.is_forward_flippable(e)
                        assert not T.is_forward_flippable(a)
                        assert not T.is_forward_flippable(c)

                    elif col == RED:
                        if T.is_forward_flippable(a):
                            recol.append((a, RED))
                            T._colouring[a] = PURPLE
                            T._colouring[T._ep[a]] = PURPLE
                        if a != c and T.is_forward_flippable(c):
                            recol.append((c, RED))
                            T._colouring[c] = PURPLE
                            T._colouring[T._ep[c]] = PURPLE

                        # assertions to be removed
                        assert not T.is_forward_flippable(e)
                        assert not T.is_forward_flippable(b)
                        assert not T.is_forward_flippable(d)

                    if verbosity >= 2:
                        print('[automaton] recol = %s' % (recol,))
                        sys.stdout.flush()

                new_iso_sig = T.iso_sig()
                if verbosity >= 2:
                    print('[automaton] it is core with iso_sig %s' % new_iso_sig)
                    sys.stdout.flush()

                if new_iso_sig not in graph:
                    if verbosity >= 2:
                        print('[automaton] new vertex')
                        sys.stdout.flush()

                    # new core
                    flips.append((e, old_col, recol))
                    graph[new_iso_sig] = []
                    graph[iso_sig].append((new_iso_sig, e, old_col, col))

                    # assertion to be removed
                    assert not reduced or old_col == PURPLE

                    iso_sig = new_iso_sig
                    iso_sigs.append(iso_sig)
                    branch.append([])

                    # Computing the new forward flippable edges should not be that complicated
                    # The edge we flip is not available anymore and we have at most two newly
                    # flippable edges in the neighborhood
                    if reduced:
                        ffe = T.purple_edges()
                        assert ffe == T.forward_flippable_edges()
                    else:
                        ffe = T.forward_flippable_edges()
                    branch[-1].extend((x, BLUE) for x in ffe)
                    branch[-1].extend((x, RED) for x in ffe)
                    e, col = branch[-1].pop()
                    continue
                else:
                    # (e,col) leads to an already visited vertex
                    if verbosity >= 2:
                        print('[automaton] known vertex')
                        sys.stdout.flush()
                    # assertion to be removed
                    assert not reduced or old_col == PURPLE
                    graph[iso_sig].append((new_iso_sig, e, old_col, col))

            else:  # = T is not core
                if reduced:
                    recol = []
                if verbosity >= 2:
                    print('[automaton] not core')
                    sys.stdout.flush()

            # not core or already visited
            if reduced:
                for ee, ccol in recol:
                    T._colouring[ee] = ccol
                    T._colouring[T._ep[ee]] = ccol
            T.flip_back(e, old_col)
            # TODO: check to be removed
            if env.CHECK:
                T._check()
            assert T.iso_sig() == iso_sigs[-1], (T.iso_sig(), iso_sigs[-1])

            while flips and not branch[-1]:
                e, old_col, recol = flips.pop()
                if reduced:
                    for ee, ccol in recol:
                        if verbosity >= 2:
                            print('[automaton] recolour %d as %s' % (ee, colour_to_char(ccol)))
                            sys.stdout.flush()
                        T._colouring[ee] = ccol
                        T._colouring[T._ep[ee]] = ccol
                T.flip_back(e, old_col)
                branch.pop()
                iso_sigs.pop()
                assert T.iso_sig() == iso_sigs[-1]
            if not branch[-1]:
                break
            iso_sig = iso_sigs[-1]
            e, col = branch[-1].pop()

        if verbosity == 1:
            print('\r[automaton] %8d      %8d      %.3f      ' % (len(graph), len(branch), time() - t0))
            sys.stdout.flush()
        return CoreAutomaton(graph)

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
