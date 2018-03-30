
from __future__ import print_function

import sys
try:
    from Queue import Queue
except ImportError:
    from queue import queue
from numpy import matrix as Matrix
from math import log

from constants import *
from colouredtriangulation import ColouredTriangulation

import flipper
norm = flipper.norm

class Automaton(object):
    def __init__(self, graph, types):
        self.graph = graph
        self.types = types
        self.iso_sigs = sorted(graph)
        self.index = dict((sig, index) for index, sig in enumerate(self.iso_sigs))
    def __str__(self):
        return '\n'.join('\n'.join('%s\t%s\t%s\t%s' % (g, n, ';'.join(','.join(map(str, row)) for row in M.tolist()), t) for n, M, t in self.graph[g]) for g in self.iso_sigs)
    def __repr__(self):
        return str(self)
    def __len__(self):
        return len(self.iso_sigs)
    def __iter__(self):
        return iter(self.iso_sigs)
    def __contains__(self, item):
        return item in self.iso_sigs
    def export(self, filepath):
        with open(filepath, 'w') as open_filepath:
            open_filepath.write('digraph MyGraph {\n')
            for g in self.iso_sigs:
                open_filepath.write('    %s [style=filled, fillcolor=%s];\n' % (g, TYPE_COLOURS[self.types[g]]))
                for n, M, t in self.graph[g]:
                    open_filepath.write('    %s -> %s;\n' % (g, n))  # , ';'.join(','.join(map(str, row)) for row in M.tolist()), t))
            open_filepath.write('}')
    
    @classmethod
    def from_triangulation(self, T, slope=VERTICAL, verbose=False):
        if verbose:
            print('Computing stratum: %s' % T.stratum())
            print('Stratum dimension: %d' % T.stratum_dimension())
            print('Orientable: %s' % T.is_abelian())
            print('Good | Bad | current | ratio (current/good):')
        
        zeta = T.zeta
        T = T.canonical()
        current = Queue()
        current.put(T)
        count = 0
        d = T.stratum_dimension()
        
        def is_full_dimension(t): # This is broken
            # assert(not t.is_core(d) or t.is_geometric(d))  # Test does core ==> geometric?
            return t.is_core(d)
            # return t.is_geometric(d)
        
        if not is_full_dimension(T):
            print('Starting triangulation is NOT full dimension!')
            print('Searching for full dimensional one.')
            seen = set()
            start = None
            while not current.empty():
                count += 1
                T = current.get()
                for i in T.flippable_edges():
                    for c in COLOURS:
                        neighbour = T.flip_edge(i, c)
                        s, _ = neighbour.best_translation()
                        if s not in seen:
                            if is_full_dimension(neighbour):
                                start = neighbour.canonical()
                            else:
                                current.put(neighbour.canonical())
                        seen.add(s)
                        if start is not None: break
                    if start is not None: break
                if start is not None: break
                
                if verbose and count % 10 == 0:
                    print('\r%d %d              ' % (len(seen), current.qsize()), end='')
                    sys.stdout.flush()
            if start is None: raise ValueError('Stratum is empty.')
            if verbose: print('Found: %s' % start)
            T = start
            current = Queue()
            current.put(T)
            del seen
        
        graph = dict()
        good = set([T.iso_sig()])
        bad = set()
        types = dict([(T.iso_sig(), T.type())])
        while not current.empty():
            count += 1
            T = current.get()
            neighbours = []
            for i in T.mostly_sloped_edges(slope):
                for c in COLOURS:
                    neighbour = T.flip_edge(i, c)
                    s, t = neighbour.best_translation()
                    if s not in good and s not in bad:
                        if is_full_dimension(neighbour):
                            good.add(s)
                            types[s] = neighbour.type()
                            current.put(neighbour.canonical())
                        else:
                            bad.add(s)
                    if s in good:
                        corner1 = T.triangulation.corner_lookup[i]
                        corner2 = T.triangulation.corner_lookup[~i]
                        M = Matrix([[1 if (j != i and k == j) or (j == i and k in [corner1.indices[1], corner2.indices[2]]) else 0 for k in range(zeta)] for j in range(zeta)])
                        N = Matrix([[1 if j == norm(t[k]) else 0 for k in range(zeta)] for j in range(zeta)])
                        
                        neighbours.append((s, N * M, 'A'))
            
            for isom in T.self_isometries():
                M = Matrix([[1 if j == isom.index_map[i] else 0 for j in range(zeta)] for i in range(zeta)])
                neighbours.append((T.iso_sig(), M, 'B'))
            
            graph[T.iso_sig()] = sorted(neighbours, key=lambda (s, A, t): (s, A.tolist(), t))
            
            if verbose and count % 1 == 0:
                print('\r%d %d %d %0.3f               ' % (len(good), len(bad), current.qsize(), float(current.qsize()) / len(good)), end='')
                sys.stdout.flush()
        
        if verbose:
            print('\r%d %d %d %0.3f               ' % (len(good), len(bad), current.qsize(), float(current.qsize()) / len(good)))
        
        return Automaton(graph, types)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Build automata')
    parser.add_argument('--sig', type=str, required=True, help='signature to load')
    args = parser.parse_args()
    
    T = ColouredTriangulation.from_iso_sig(args.sig)
    A = Automaton.from_triangulation(T, verbose=True)
    print(len(A))
    A.export('graphs/' + args.sig + '.dot')

