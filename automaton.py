
from __future__ import print_function

import sys
try:
    from Queue import Queue
except ImportError:
    from queue import Queue

from colouredtriangulation import ColouredTriangulation, ngon
from constants import *

import flipper
norm = flipper.norm

class Automaton(object):
    def __init__(self, graph):
        self.graph = graph
        self.iso_sigs = sorted(graph)
        self.index = dict((sig, index) for index, sig in enumerate(self.iso_sigs))
    def __str__(self):
        return '\n'.join('%s\t%s' % (g, self.graph[g][0]) for g in self.iso_sigs)
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
                neighbours, node_type = self.graph[g]
                G = ColouredTriangulation.from_iso_sig(g)
                aut_size = len(G.self_isometries()) 
#                open_filepath.write('    %s [label="", style=filled, color="%s"];\n' % (g, TYPE_COLOURS[node_type]))
                open_filepath.write('    %s [label="%d", style=filled, color="%s"];\n' % (g, aut_size, TYPE_COLOURS[node_type]))
                for n, c, w in neighbours:
                    open_filepath.write('    %s -> %s [color="%s", penwidth="%s", arrowsize="%s"];\n' % (g, n, c, w, w))
            open_filepath.write('}\n')
    
    @classmethod
    def from_triangulation(self, T, verbose=False, mode='core', **kwargs):
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
            if mode == 'core':
                T_neighbours = T.neighbours_core(**kwargs)
            else:
                T_neighbours = T.neighbours_geometric(**kwargs)
            for neighbour, change in T_neighbours:
                s = neighbour.iso_sig()
                if s not in seen:
                    seen.add(s)
                    current.put(neighbour.canonical())
                
                neighbours.append((s, change, 1))
            
#            for isom in T.self_isometries():
#                neighbours.append((T.iso_sig(), 'Gray', 0.1))
            
            graph[T.iso_sig()] = (sorted(neighbours), T.type())
            
            if verbose and count % 1 == 0:
                print('\r%d %d               ' % (len(seen), current.qsize()), end='')
                sys.stdout.flush()
        
        if verbose:
            print('\r%d %d               \n' % (len(seen), current.qsize()), end='')
            sys.stdout.flush()
        
        return Automaton(graph)

if __name__ == '__main__':
    # Abelian strata, ordered by dimension
    # Ordered by dimension, abelian first
    # Dimension 2 
    T = ColouredTriangulation.from_pA(flipper.load('S_1_1').mapping_class('aB'))  # H(0) # 2 nodes
    # T = ColouredTriangulation.from_pA(flipper.load('SB_4').mapping_class('s_0S_1'))  # Q(-1^4) # 6 nodes
    # Dimension 3
    # T = ColouredTriangulation.from_pA(flipper.load('S_1_2').mapping_class('aaaBc'))  # H(0^2) # 16 nodes
    # Q(2, -1^2) #
    # Dimension 4
    # T = ngon(8)  # [4] # 86 nodes    ????
    # Dimension 5
    # T = ngon(10)  # [2,2] # 876 nodes  ????
    
    # T = ColouredTriangulation.from_pA(flipper.load('S_1_2').mapping_class('abC'))  # [1, 1, -1, -1] # 1658 nodes
    # T = ColouredTriangulation.from_pA(flipper.load('S_2_1').mapping_class('a.a.a.a.B.c.d'))  # [2, 1, 1] # 16896 nodes
    ##
    ## Tri = flipper.create_triangulation([(~4,5,0),(~5,1,6),(~6,~1,7),(~7,2,8),(~8,~2,9),(~9,3,10),(~10,~3,11),(~11,4,~0)])
    ## colours = [BLUE]*4 + [RED]*8
    ## T = ColouredTriangulation(Tri, colours) # [1, -1, -1, -1, -1, -1] # 530 nodes
    ##
    ## from surface_dynamics import *
    ## Q = QuadraticStratum(5,-1)
    ## c, = Q.components()
    ## T, C = cyl_diag_to_veering.to_veering_triangulation(c.one_cylinder_diagram())
    ## T = ColouredTriangulation(T, C)  # [5, -1] # 4736
    ##
    ## Q = QuadraticStratum(3,-1,-1,-1)
    ## c, = Q.components()
    ## T, C = cyl_diag_to_veering.to_veering_triangulation(c.one_cylinder_diagram())
    ## T = ColouredTriangulation(T, C)  # [3, -1, -1, -1] # 4246
    A = Automaton.from_triangulation(T, verbose=True)
    A.export('graphs/%s.dot' % ','.join(str(x) if x >= 0 else 'm%s' % abs(x) for x in T.stratum()))



