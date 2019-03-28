"""
Look at nice patterns of vertex cycles on a train-track.

We should probably be looking at nice patterns of *carried* curves!
"""
from veerer import *
from surface_dynamics import *

from sage.graphs.generic_graph_pyx import SubgraphSearch

def longest_chain(T):
    '''Return the longest chain of vertex cycles on the veering triangulation T'''
    cycles = T.vertex_cycles()
    n = len(cycles)
    adj = [x.intersection(y) for y in cycles for x in cycles]

    I = matrix(ZZ, n, adj)
    J = matrix(ZZ, n, [bool(x) for x in adj])
    G = Graph(J)

    m = 2
    path = None
    while True:
        for V in SubgraphSearch(G, graphs.PathGraph(m), induced=True):
            P = G.subgraph(V)
            edges = P.edges(False)
            if all(I[i,j] == 1 for (i,j) in edges):
                print("found a path of length {}".format(m))
                m += 1
                path = edges[:]
                break
        else:
            print("no path of length {}".format(m))
            return path
