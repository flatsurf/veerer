
import flipper

from colouredtriangulation import ColouredTriangulation
from structures import UnionFind
from constants import *
from automaton import Automaton

class SquaredTriangulation(ColouredTriangulation):
    def __init__(self, triangulation, colouring, slope=HORIZONTAL, sanity=False):
        super(SquaredTriangulation, self).__init__(triangulation, colouring, sanity)
        self.diagonals = set(self.mostly_sloped_edges(slope))
        self.colour = self.colouring[list(self.diagonals)[0]]
        assert all(self.colouring[diagonal] == self.colour for diagonal in self.diagonals)
        
        classes = UnionFind(self.triangulation.indices)
        for diagonal in self.diagonals:
            classes.union(diagonal, *[edge.index for edge in self.triangulation.square_about_edge(diagonal) if self.colouring[edge.index] == self.colour])
        self.curves = [cls for cls in classes if len(cls) > 1]
    
    def twist(self, curve):
        ''' Perform a (root of a) twist about this curve. '''
        assert curve in self.curves
        
        tmp = self
        for edge in self.diagonals.intersection(curve):
            tmp = tmp.flip_edge(edge, self.colour)
        
        return SquaredTriangulation(tmp.triangulation, tmp.colouring)
    
    def swap(self):
        tmp = self
        for edge in self.diagonals:
            tmp = tmp.flip_edge(edge, RED if self.colour == BLUE else BLUE)
        
        return SquaredTriangulation(tmp.triangulation, tmp.colouring)
    
    def neighbours(self):
        return [self.twist(curve) for curve in self.curves] + [self.swap()]

def test():
    T = SquaredTriangulation(flipper.create_triangulation([(~5, ~3, 4), (~4, 5, 0), (~2, 1, 3), (~1, ~0, 2)]), [BLUE, BLUE, RED, BLUE, RED, BLUE])
    print(T.iso_sig())

if __name__ == '__main__':
    # RBBBRB_acbdklejfghi
    import argparse
    parser = argparse.ArgumentParser(description='Build automata')
    parser.add_argument('--sig', type=str, required=True, help='signature to load')
    args = parser.parse_args()
    
    T = SquaredTriangulation.from_iso_sig(args.sig)
    A = Automaton.from_triangulation(T, verbose=True)
    A.export('graphs/' + args.sig + '.dot')

