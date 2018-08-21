
from __future__ import print_function, absolute_import

import sys
from time import time

from .coloured_triangulation import ColouredTriangulation, ngon
from .constants import *

class Automaton(object):
    r"""
    EXAMPLES::

        sage: from veerer import *
        sage: from surface_dynamics import *

        sage: T = ColouredTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
        sage: A = Automaton.from_triangulation(T)
        sage: len(A)
        2

        sage: T = ColouredTriangulation.from_stratum(AbelianStratum(2))
        sage: A = Automaton.from_triangulation(T)
        sage: len(A)
        86

        sage: T = ColouredTriangulation.from_stratum(QuadraticStratum(2,-1,-1))
        sage: A = Automaton.from_triangulation(T)
        sage: len(A)
        160

        sage: T = ColouredTriangulation.from_stratum(QuadraticStratum(2,2))
        sage: A = Automaton.from_triangulation(T)
        sage: len(A)
        846

        sage: T = ColouredTriangulation.from_stratum(AbelianStratum(1,1))
        sage: A = Automaton.from_triangulation(T)
        sage: len(A)
        876
    """
    def __init__(self, graph):
        self._graph = graph
        self._iso_sigs = sorted(graph)
        self._index = dict((sig, index) for index, sig in enumerate(self._iso_sigs))
    def __str__(self):
        return '\n'.join('%s\t%s' % (g, self._graph[g][0]) for g in self._iso_sigs)
    def __repr__(self):
        return str(self)
    def __len__(self):
        return len(self._iso_sigs)
    def __iter__(self):
        return iter(self._iso_sigs)
    def __contains__(self, item):
        return item in self._iso_sigs
    def export(self, filepath):
        with open(filepath, 'w') as open_filepath:
            open_filepath.write('digraph MyGraph {\n')
            for g in self._iso_sigs:
                neighbours, node_type = self._graph[g]
                G = ColouredTriangulation.from_iso_sig(g)
                aut_size = len(G.self_isometries()) 
#                open_filepath.write('    %s [label="", style=filled, color="%s"];\n' % (g, TYPE_COLOURS[node_type]))
                open_filepath.write('    %s [label="%d", style=filled, color="%s"];\n' % (g, aut_size, TYPE_COLOURS[node_type]))
                for n, c, w in neighbours:
                    open_filepath.write('    %s -> %s [color="%s", penwidth="%s", arrowsize="%s"];\n' % (g, n, c, w, w))
            open_filepath.write('}\n')
    
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
        n = T._triangulation.num_edges()
        d = T.stratum_dimension()
        count = 0

        graph = {iso_sig: []}
        iso_sigs = [iso_sig]  # iso_sig sequence
        flips = []     # backward flip sequence
        branch = []    # forward flips available

        ffe = T.forward_flippable_edges()
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
                    graph[iso_sig].append(new_iso_sig)

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
                    graph[iso_sig].append(new_iso_sig)
            elif verbose >= 2:
                print('[automaton] not core')

            # not core or already visited
            T.back_flip(e, old_col)
            assert T.iso_sig() == iso_sigs[-1]

            while flips and not branch[-1]:
                e,old_col = flips.pop()
                T.back_flip(e, old_col)
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
