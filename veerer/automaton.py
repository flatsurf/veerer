r"""
Train-track and triangulations automata.
"""
from __future__ import print_function, absolute_import

import os, sys, shutil
from time import time

from .veering_triangulation import VeeringTriangulation
from .constants import *

from .env import sage

class Automaton(object):
    r"""
    Automaton of core veering triangulations.

    EXAMPLES::

        sage: from veerer import *
        sage: from surface_dynamics import *

        sage: T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
        sage: A = Automaton.from_triangulation(T)
        sage: len(A)
        2

        sage: Automaton.from_stratum(AbelianStratum(2))
        Core veering automaton with 86 vertices

        sage: Automaton.from_stratum(QuadraticStratum(2,-1,-1))
        Core veering automaton with 160 vertices

        sage: Automaton.from_stratum(QuadraticStratum(2,2))
        Core veering automaton with 846 vertices

        sage: Automaton.from_stratum(AbelianStratum(1,1))
        Core veering automaton with 876 vertices

    Quadratic strata::

        sage: T = VeeringTriangulation('(0,1,2)', 'BBR')
        sage: A = Automaton.from_triangulation(T)
        sage: len(A)
        2

    The examples below is Q(1^2, -1^2) with two folded edges::

        sage: fp = '(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)'
        sage: cols = 'BRBRBBBRBR'
        sage: T = VeeringTriangulation(fp, cols)
        sage: A = Automaton.from_triangulation(T)
        sage: len(A)
        1074
    """
    def __init__(self, graph):
        self._graph = graph
        self._iso_sigs = sorted(graph)
        self._index = dict((sig, index) for index, sig in enumerate(self._iso_sigs))
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
            sage: from surface_dynamics import *

            sage: fp = "(0,~7,6)(1,~8,~2)(2,~6,~3)(3,5,~4)(4,8,~5)(7,~1,~0)"
            sage: cols = 'RBRBRBBBB'
            sage: T = VeeringTriangulation(fp, cols)
            sage: A = Automaton.from_triangulation(T)
            sage: A
            Core veering automaton with 86 vertices

            sage: A.to_graph()
            Looped multi-digraph on 86 vertices

            sage: A.to_graph(directed=False, multiedges=False, loops=True)
            Looped graph on 86 vertices
        """
        if sage is None:
            raise ValueError('Only available in SageMath')
        elif directed:
            from sage.graphs.digraph import DiGraph
            G = DiGraph(loops=loops, multiedges=multiedges)
        else:
            from sage.graphs.graph import Graph
            G = Graph(loops=loops, multiedges=multiedges)

        for g, neighb in self._graph.items():
            for gg,_,_ in neighb:
                G.add_edge(g, gg)

        return G

    def rotation_automorphism(self):
        r"""
        Return the automorphism of the vertices that corresponds to rotation.

        Note that this automorphism reverses the edge direction.

        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: T = VeeringTriangulation.from_stratum(AbelianStratum(2))
            sage: A = Automaton.from_triangulation(T)
            sage: rot = A.rotation_automorphism()

            sage: G = A.to_graph()
            sage: all(G.has_edge(rot[b], rot[a]) for a,b in G.edges(False))
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
            sage: from surface_dynamics import *

            sage: T = VeeringTriangulation.from_stratum(AbelianStratum(2))
            sage: A = Automaton.from_triangulation(T)
            sage: conj = A.conjugation_automorphism()

            sage: G = A.to_graph()
            sage: all(G.has_edge(conj[a], conj[b]) for a,b in G.edges(False))
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
            sage: from surface_dynamics import *

            sage: filename = tmp_filename() + '.dot'
            sage: T = VeeringTriangulation.from_stratum(AbelianStratum(2))
            sage: A = Automaton.from_triangulation(T)
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

        seed = min(self._graph)
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

            typ = 0
            if T.is_geometric():
                typ |= GEOMETRIC
            if T.is_cylindrical():
                typ |= CYLINDRICAL
            colour = TYPE_COLOURS[typ]

            aut_size = len(T.automorphisms())
            if triangulations:
                f.write("""    %s [label="%d" style=filled color="%s" tooltip="%s" href="%s"];\n""" % (g, aut_size, colour, g, t_rel_filename))
            else:
                 f.write('    %s [label="%d" color="%s"];\n' % (g, aut_size, colour))

            for gg, old_col, new_col in self._graph[g]:
                f.write('    %s -> %s [color="%s;%f:%s"];\n' % (g, gg, old_col, 0.5, new_col))
        f.write('}\n')

        if filename is not None:
            f.close()

    def geometric_triangulations(self, method=None):
        r"""
        Return an iterator over the pairs (veering triangulation,
        geometric polytope) when the triangulation is geometric.

        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: A = Automaton.from_stratum(AbelianStratum(2))
            sage: vt, P = next(A.geometric_triangulations())
            sage: vt.geometric_polytope() == P
            True
            sage: P.affine_dimension() == 8
            True

            sage: sum(1 for _ in A.geometric_triangulations())
            54
        """
        dim = self.one_triangulation().stratum().dimension()
        for vt in self:
            p = vt.geometric_polytope()
            if p.affine_dimension() == 2*dim:
                yield vt, p

    def cylindrical_triangulations(self):
        r"""
        Return an iterator over cylindrical configurations.

        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import *

            sage: T = VeeringTriangulation.from_stratum(AbelianStratum(2))
            sage: A = Automaton.from_triangulation(T)

            sage: vt = next(A.cylindrical_triangulations())
            sage: sum(len(u) for u in vt.cylinders(BLUE)) == 6 or \
            ....: sum(len(u) for u in vt.cylinders(RED)) == 6
            True

            sage: sum(1 for _ in A.cylindrical_triangulations())
            24
        """
        for vt in self:
            if vt.is_cylindrical():
                yield vt

    @classmethod
    def from_triangulation(self, T, verbose=0, mode='core', **kwargs):
        assert T.is_core()

        if verbose:
            print('[automaton] stratum: %s' % T.stratum())
            print('[automaton] stratum dimension: %d' % T.stratum().dimension())
        if verbose == 1:
            print('[automaton] size(graph)   size(path)    time')

        T = T.copy()
        iso_sig = T.iso_sig()
        d = T.stratum_dimension()
        count = 0

        graph = {iso_sig: []}
        iso_sigs = [iso_sig]  # iso_sig sequence
        flips = []     # backward flip sequence
        branch = []    # forward flips available

        ffe = T.forward_flippable_edges()
        ffe_orb = []
        branch.append([])
        branch[-1].extend((x,BLUE) for x in ffe)
        branch[-1].extend((x,RED) for x in ffe)
        e,col = branch[-1].pop()
        old_size = 0
        t0 = time()
        while True:
            if verbose >= 2:
                print('[automaton] NEW LOOP')
                print('[automaton] at %s with iso_sig %s' % (T.to_string(), iso_sig))
                print('[automaton] going along (%s, %s)' % (e, col))
            elif verbose:
                if len(graph) > old_size + 500:
                    print('\r[automaton] %8d      %8d      %.3f      ' % (len(graph),len(branch),time()-t0), end='')
                    sys.stdout.flush()
                    old_size = len(graph)

            # some safety check to be disabled
            #assert len(flips) + 1 == len(iso_sigs) == len(branch)
            #assert iso_sig == iso_sigs[-1]
            #assert T.iso_sig() == iso_sig

            # go down in DFS
            old_col = T.colour(e)
            T.flip(e, col)

            if verbose >= 2:
                print('[automaton] ... landed at %s' % T.to_string())

            if T.edge_has_curve(e): # = T is core
                new_iso_sig = T.iso_sig()
                if verbose >= 2:
                    print('[automaton] it is core with iso_sig %s' % new_iso_sig)


                if new_iso_sig not in graph:
                    if verbose >= 2:
                        print('[automaton] new vertex')
                    # new core
                    flips.append((e,old_col))
                    graph[new_iso_sig] = []
                    graph[iso_sig].append((new_iso_sig,old_col,col))

                    iso_sig = new_iso_sig
                    iso_sigs.append(iso_sig)
                    branch.append([])
                    ffe = T.forward_flippable_edges()
                    branch[-1].extend((x,BLUE) for x in ffe)
                    branch[-1].extend((x,RED) for x in ffe)
                    e,col = branch[-1].pop()
                    continue
                else:
                    # (e,col) leads to an already visited vertex
                    if verbose >= 2:
                        print('[automaton] known vertex')
                    graph[iso_sig].append((new_iso_sig,old_col,col))
            elif verbose >= 2:
                print('[automaton] not core')

            # not core or already visited
            T.flip_back(e, old_col)
            assert T.iso_sig() == iso_sigs[-1]

            while flips and not branch[-1]:
                e,old_col = flips.pop()
                T.flip_back(e, old_col)
                branch.pop()
                iso_sigs.pop()
                assert T.iso_sig() == iso_sigs[-1]
            if not branch[-1]:
                break
            iso_sig = iso_sigs[-1]
            e,col = branch[-1].pop()

        if verbose == 1:
            print('\r[automaton] %8d      %8d      %.3f      ' % (len(graph),len(branch),time()-t0))
        return Automaton(graph)

    @classmethod
    def from_stratum(self, stratum):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: from surface_dynamics import *
            sage: Automaton.from_stratum(AbelianStratum(2))
            Core veering automaton with 86 vertices
        """
        return self.from_triangulation(VeeringTriangulation.from_stratum(stratum))
