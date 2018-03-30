
import math
import flipper
from vector import Vector2
from utilities import UnionFind
from constants import *

norm = flipper.norm
EPSILON = 10**-8
CUTOFF = 1.0

class FlatTriangulation(object):
    def __init__(self, triangulation, vectors):
        self.triangulation = triangulation
        self.zeta = self.triangulation.zeta
        self.vectors = vectors  # A dict mapping labels --> vectors.
        
        assert(all(sum([self(edge) for edge in triangle], Vector2(0,0)).approx(Vector2(0, 0)) for triangle in self.triangulation))  # FPE.
        assert(all(self(edge).approx(-self(~edge)) or self(edge).approx(self(~edge)) for edge in self.triangulation.edges))  # FPE.
    
    @classmethod
    def from_pA(self, h):
        F = h.flat_structure()
        T = F.triangulation
        vectors = dict((edge.label, Vector2(float(F.edge_vectors[edge].x), float(F.edge_vectors[edge].y))) for edge in T.edges)
        F = FlatTriangulation(T, vectors)
        
        while True:  # Flip to make Delaunay.
            for i in F.forwards_flippable_edges():
                if math.log(F.aspect_ratio(i)) < 0:
                    F = F.flip_edge(i)
                    break
            else:
                break
        
        F = F.orient()
        
        return F
    
    def __str__(self):
        return str(self.triangulation) + '\n\t' + '\n\t'.join('%d: %s  %s' % (i, self.vectors[i], self.vectors[~i]) for i in range(self.zeta))
    def __repr__(self):
        return str(self)
    def __iter__(self):
        return iter(self.vectors)
    def __call__(self, edge):
        return self.vectors[edge.label]
    
    def dual_tree(self, possible=None):
        if possible is None: possible = self.triangulation.indices
        
        # Kruskal's algorithm.
        dual_tree = [False] * self.zeta
        classes = UnionFind(self.triangulation)
        for index in possible:
            a, b = self.triangulation.triangle_lookup[index], self.triangulation.triangle_lookup[~index]
            if classes(a) != classes(b):
                classes.union(a, b)
                dual_tree[index] = True
        
        return dual_tree
    
    def orient(self):
        dual_tree = self.dual_tree()
        
        vectors = dict(self.vectors)
        todo = [self.triangulation.triangles[0]]  # Stack.
        fixed = set(todo)
        while todo:
            triangle = todo.pop()
            for side, index in zip(triangle.labels, triangle.indices):
                if dual_tree[index]:
                    neighbour = self.triangulation.triangle_lookup[~side]
                    if neighbour not in fixed:
                        if not vectors[~side].approx(-vectors[side]):
                            for label in neighbour.labels:
                                vectors[label] = -vectors[label]
                        todo.append(neighbour)
                        fixed.add(neighbour)
        
        return FlatTriangulation(self.triangulation, vectors)
    
    def tighten(self):
        vectors = dict()
        
        reversed_gluing = [i for i in self.triangulation.indices if self.vectors[i].approx(-self.vectors[~i])]
        same_gluing = [i for i in self.triangulation.indices if self.vectors[i].approx(self.vectors[~i])]
        
        dual_tree = self.dual_tree(possible=reversed_gluing)
        
        for index in self.triangulation.indices:
            if not dual_tree[index]:
                vectors[index] = self.vectors[index]
                vectors[~index] = -vectors[index] if self.vectors[index].approx(-self.vectors[~index]) else vectors[index]
        
        if not (all(i in vectors for i in same_gluing)):
            print(self)
            print(vectors)
            print(reversed_gluing)
            print(same_gluing)
            print(dual_tree)
        
        assert(all(i in vectors for i in same_gluing))
        
        if same_gluing:  # Correct boundary, happens iff non-Abelian case.
            vector = sum(vectors.values(), Vector2(0, 0))
            index = same_gluing[0]
            vectors[index] -= 0.5 * vector
            vectors[~index] -= 0.5 * vector
        
        while True:  # Could do a topological sort to avoid quadratic behaviour.
            for label in self.triangulation.labels:
                if label not in vectors:
                    corner = self.triangulation.corner_lookup[label]
                    a, b = corner.labels[1:]
                    if a in vectors and b in vectors:
                        vectors[label] = -vectors[a] - vectors[b]
                        vectors[~label] = -vectors[label] if self.vectors[label].approx(-self.vectors[~label]) else vectors[label]
                        break
            else:
                break
        
        #for index in range(self.zeta):
            #print('\t\tTweaked %d by %0.10f   %s' % (norm(index), (vectors[index] - self.vectors[index]).norm(), vectors[index] - self.vectors[index]))
        
        return FlatTriangulation(self.triangulation, vectors)
    
    def colours_about_edge(self, i):
        return [self.vectors[edge.label].colour() for edge in self.triangulation.square_about_edge(i)]
    
    def alternating_square(self, i):
        colours = self.colours_about_edge(i)
        return all(colours[j] != colours[(j+1) % 4] for j in range(4))
    
    def is_flippable(self, i):
        return self.triangulation.is_flippable(i) and self.alternating_square(i)
    
    def flippable_edges(self):
        return [i for i in self.triangulation.flippable_edges() if self.is_flippable(i)]
    
    def forwards_flippable_edges(self):
        return [i for i in range(self.zeta) if self.colours_about_edge(i) == [BLUE, RED, BLUE, RED]]
    
    def flip_edge(self, i):
        assert(self.is_flippable(i))
        
        T = self.triangulation.flip_edge(i)
        
        a, b, c, d = [edge.label for edge in self.triangulation.square_about_edge(i)]
        vectors = dict(self.vectors)
        if not vectors[i].approx(-vectors[~i]):  # Make these two triangles compatible.
            vectors[~i] = -vectors[~i]
            vectors[c] = -vectors[c]
            vectors[d] = -vectors[d]
        vectors[i] = -vectors[d] - vectors[a]
        vectors[~i] = -vectors[b] - vectors[c]
        
        return FlatTriangulation(T, vectors)
    
    def aspect_ratio(self, i):
        a, b, c, d = self.triangulation.square_about_edge(i)
        
        h = abs(self(b).y) + abs(self(c).y)
        w = abs(self.vectors[i].x)
        
        return h / w
    
    def flow(self, time):
        vectors = dict((k, v.flow(time)) for k, v in self.vectors.items())
        return FlatTriangulation(self.triangulation, vectors)
    
    def max_flowable_time(self):
        times = [(math.log(self.aspect_ratio(i)), i) for i in self.forwards_flippable_edges()]
        assert(all(time > 0 for time in times))
        max_time = min(time for time, i in times)
        return max_time, [i for time, i in times if abs(time - max_time) < EPSILON**2]  # FPE.
    
    def Siegel_Veech_estimates(self, num_steps=None, max_time=None, debug=False):
        clock = 0.0
        count = 0
        last_changed = [0.0 for _ in range(self.zeta)]
        max_completed = 0.0
        sum_completed = 0.0
        lifetimes = []
        F = self.orient()
        while True:
            count += 1
            if debug:
                print('=================== %d at %0.4fs ===================' % (count, clock))
                print('Changes: %s' % last_changed)
                print('Sum completed: %f' % sum_completed)
                print('Max completed: %f' % max_completed)
                print('====================================================')
            t, to_flip = F.max_flowable_time()
            if debug: print('Need to flow for %0.4fs' % t)
            clock += t
            F = F.flow(t)
            if debug:
                print('After flow:')
                print(F)
            F = F.tighten()
            if debug:
                print('After tighen:')
                print(F)
            if debug: print('Need to flip %s' % to_flip)
            for index in to_flip:
                vector = F.vectors[index]
                birth, death = vector.short_times(clock)
                real_birth, real_death = max(last_changed[index], birth), min(clock, death)
                min_length2 = vector.min_length2() # if real_birth <= vector.min_time(clock) <= real_death else min(??)
                
                if min_length2 < CUTOFF**2:
                                        # This will do for now, but we may want to estimate the shear instead.
                    excursion = 2 * math.sqrt((CUTOFF**4/min_length2**2) - 1)
                    max_completed = max(max_completed, excursion)
                    sum_completed = sum_completed + excursion
                last_changed[index] = clock
                F = F.flip_edge(index)
            if debug:
                print('After flip:')
                print(F)
            F = F.orient()
            F = F.tighten()
            if debug:
                print('After retighten:')
                print(F)
            
            
            sum_still = 0.0
            max_still = 0.0
            for index in range(self.zeta):
                if index not in to_flip:
                    vector = F.vectors[index]
                    birth, death = vector.short_times(clock)
                    real_birth, real_death = max(last_changed[index], birth), min(clock, death)
                    if real_birth <= vector.min_time(clock) <= real_death:
                        excursion = 2*math.sqrt(1/vector.min_length2()**2 - 1)
                        max_still = max(max_still, excursion)
                        sum_still = sum_still + excursion
            
            yield clock, max(max_completed, max_still), (sum_completed + sum_still - max(max_completed, max_still)) / (clock * math.log(clock)) / CUTOFF  # Could also yield error term.
            
            if num_steps is not None and count > num_steps: break
            if max_time is not None and clock > max_time: break
        
        return


if __name__ == '__main__':
    h = flipper.load('S_1_1').mapping_class('aB')
    # h = flipper.load('S_1_2').mapping_class('aBBcxxcBB')
    F = FlatTriangulation.from_pA(h)
    for index, (clock, mx, estimate) in enumerate(F.Siegel_Veech_estimates(num_steps=5000000, max_time=5000000, debug=False)):
        if index % 100 == 0:
            print(index, '%0.4fs' % clock, mx, estimate)

