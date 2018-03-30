
from __future__ import print_function

import sys
try:
    from Queue import Queue
except ImportError:
    from queue import Queue
from numpy import matrix as Matrix

from constants import *

import flipper
norm = flipper.norm

class Automaton(object):
    def __init__(self, graph):
        self.graph = graph
        self.iso_sigs = sorted(graph)
        self.index = dict((sig, index) for index, sig in enumerate(self.iso_sigs))
    def __str__(self):
        return '\n'.join('\n'.join('%s\t%s' % (g, n) for n, in self.graph[g]) for g in self.iso_sigs)
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
                for n, in self.graph[g]:
                    open_filepath.write('    %s -> %s;\n' % (g, n))  # , ';'.join(','.join(map(str, row)) for row in M.tolist()), t))
            open_filepath.write('}\n')
    
    @classmethod
    def from_triangulation(self, T, verbose=False, **kwargs):
        if verbose:
            print('Computing stratum: %s' % T.stratum())
            print('Stratum dimension: %d' % T.stratum_dimension())
            print('Orientable: %s' % T.is_abelian())
            print('Is Core: %s' % T.is_core())
            print('Good | current:')
        
        zeta = T.zeta
        T = T.canonical()
        current = Queue()
        current.put(T)
        count = 0
        d = T.stratum_dimension()
        
        graph = dict()
        seen = set([T.iso_sig()])
        while not current.empty():
            count += 1
            T = current.get()
            neighbours = []
            for neighbour in T.neighbours(**kwargs):
                s = neighbour.iso_sig()
                if s not in seen:
                    seen.add(s)
                    current.put(neighbour.canonical())
                
                neighbours.append((s,))
            
            for isom in T.self_isometries():
                neighbours.append((T.iso_sig(),))
            
            graph[T.iso_sig()] = sorted(neighbours) #, key=lambda (s, A, t): (s, A.tolist(), t))
            
            if verbose and count % 1 == 0:
                print('\r%d %d               ' % (len(seen), current.qsize(), end='')
                sys.stdout.flush()
        
        if verbose:
            print('\r%d %d               ' % (len(seen), current.qsize(), end='')
            sys.stdout.flush()
        
        return Automaton(graph)

