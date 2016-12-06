
from __future__ import print_function

import sys
import string
from Queue import Queue


try:
	from pyparma import Polyhedron
	from pyparma.utils import intize
except ImportError:
	print('PyParma unavailable.')
	pass

import flipper
norm = flipper.norm

RED, BLUE = 'Red', 'Blue'
HORIZONTAL, VERTICAL = 'Horizontal', 'Vertical'

def best_rotation(X):
	return min(X[i:] + X[:i] for i in range(len(X)))

class ColouredTriangulation(object):
	def __init__(self, triangulation, colouring, sanity=False):
		self.triangulation = triangulation
		self.colouring = colouring  # A dict : edge_indices --> {Red, Blue}
		self.zeta = self.triangulation.zeta
		if sanity:
			assert(all(self.colouring[i] in [RED, BLUE] for i in self.triangulation.indices))
			assert(all(0 < [self.colouring[i] for i in t.indices].count(RED) < 3 for t in self.triangulation))
	
	@classmethod
	def from_pA(self, h):
		F = h.flat_structure()
		T = F.triangulation
		colouring = dict((edge.index, RED if (F.edge_vectors[edge].x > 0) == (F.edge_vectors[edge].y > 0) else BLUE) for edge in F.edge_vectors)
		return ColouredTriangulation(T, colouring)
	
	@classmethod
	def from_QD(self, QD):
		return NotImplemented
	
	@classmethod
	def from_iso_sig(self, sig):
		# sig == (Edge colourings | triangles)
		# Hence len(sig) = zeta + F(S) = 5 / 3 zeta
		zeta = 3 * len(sig) // 5
		T = flipper.create_triangulation(sig[zeta:])
		colouring = dict(zip(range(zeta), sig[:zeta]))
		return ColouredTriangulation(T, colouring)
	
	def __str__(self):
		return str(self.triangulation) + ', ' + str(self.colouring)
	def __repr__(self):
		return str(self)
	
	def full_dimension(self):
		return 1 + self.zeta // 3  # Check this with Saul.
	
	def colours_about_edge(self, i):
		return [self.colouring[e.index] for e in self.triangulation.square_about_edge(i)]
	
	def alternating_square(self, i):
		colours = self.colours_about_edge(i)
		return all(colours[j] != colours[(j+1) % 4] for j in range(4))
	
	def is_flippable(self, i):
		return self.triangulation.is_flippable(i) and self.alternating_square(i)
	
	def flippable_edges(self):
		return [i for i in self.triangulation.flippable_edges() if self.is_flippable(i)]
	
	def flip_edge(self, i, colour):
		assert(self.is_flippable(i))
		
		T = self.triangulation.flip_edge(i)
		colouring = dict(self.colouring)
		colouring[norm(i)] = colour
		return ColouredTriangulation(T, colouring)
	
	def matrix(self, orientation=VERTICAL):
		M = []
		for t in self.triangulation:
			corner = [corner for corner in t.corners if self.colouring[corner.indices[1]] == self.colouring[corner.indices[2]]][0]
			I = corner.indices
			corner_colour = self.colouring[I[1]]
			if (orientation == VERTICAL) != (corner_colour == RED):
				# Vertical and Blue.
				# Horizontal and Red.
				# I[1] == I[0] + I[2].
				row = [1 if i == I[1] else -1 if i == I[0] or i == I[2] else 0 for i in range(self.zeta)]
			else:
				# Horizontal and Red.
				# Vertical and Blue.
				# I[2] == I[0] + I[1].
				row = [1 if i == I[2] else -1 if i == I[0] or i == I[1] else 0 for i in range(self.zeta)]
			M.append(row)
		return M
	
	def dimension(self, orientation=VERTICAL):
		# Uses PyParma.
		M = self.matrix(orientation)
		A = [[0] + [0] * i + [1] + [0] * (self.zeta - i - 1) for i in range(self.zeta)]
		for row in M:
			A.append([0] + [i for i in row])
			A.append([-0] + [-i for i in row])
		
		return Polyhedron(hrep=intize(A)).poly.affine_dimension()
	
	def is_full_dimensional(self, orientation=None):
		if orientation is None:
			return all(self.dimension(orientation) == 1 + self.zeta // 3 for orientation in [HORIZONTAL, VERTICAL])
		else:
			return self.dimension(orientation) == 1 + self.zeta // 3
	
	def vertex_data_dict(self):
		return dict((CC[0].vertex, (len([c for c in CC if self.colouring[c.indices[1]] == BLUE and self.colouring[c.indices[2]] == RED]) - 2, len(CC))) for CC in self.triangulation.corner_classes)
	
	def stratum(self):
		return sorted([d for d, v in self.vertex_data_dict().values()], reverse=True)
	
	def good_starts(self):
		VD = self.vertex_data_dict()
		
		best = min(VD.values())
		def test(edge):
			if VD[edge.source_vertex] != best: return False
			if self.colouring[edge.index] != RED: return False
			if self.colouring[self.triangulation.corner_lookup[edge.label].indices[2]] != BLUE: return False
			return True
		
		return [edge.label for edge in self.triangulation.edges if test(edge)]
	
	def iso_sig(self, start_edges=None):
		best = None
		if start_edges is None: start_edges = self.triangulation.labels
		for start_edge in start_edges:
			translate = {start_edge: 0, ~start_edge: ~0}
			num_seen_edges = 1
			
			to_process = Queue()
			to_process.put(start_edge)
			to_process.put(~start_edge)
			while not to_process.empty():
				current_edge = to_process.get()
				
				for child in self.triangulation.corner_lookup[current_edge].labels[1:]:
					if child not in translate:
						translate[child] = num_seen_edges
						translate[~child] = ~num_seen_edges
						num_seen_edges += 1
						to_process.put(~child)
			
			X = tuple(sorted(tuple(best_rotation([translate[edge.label] for edge in t])) for t in self.triangulation))
			
			inv_translate = dict((v, k) for (k, v) in translate.items())
			Y = tuple([self.colouring[norm(inv_translate[i])] for i in range(self.triangulation.zeta)])
			
			if best is None or Y + X < best:
				best = Y + X
		
		return best

def test():
	T = ColouredTriangulation(flipper.create_triangulation([(0,1,2), (~0,~1,~2)]), {0: RED, 1: BLUE, 2: RED})
	
	print(T)
	print(T.stratum())
	print(T.iso_sig())
	T2 = ColouredTriangulation.from_iso_sig(T.iso_sig())
	print(T2.iso_sig())
	M = T.matrix()
	for row in M:
		print(row)
	print(T.dimension())
	
	print(T.flippable_edges())
	T2 = T.flip_edge(0, BLUE)
	print(T2)
	print(T2.stratum())
	print(T2.iso_sig())

def ngon(n):
	assert(n > 4 and n % 2 == 0)
	m = (n - 2) // 2
	T = flipper.create_triangulation( \
		[(i, 2*m+i, ~(i+1)) for i in range(m)] + \
		[(~0, ~(2*m), ~(m+1))] + \
		[(m+i, ~(2*m+i), ~(m+i+1)) for i in range(1, m-1)] + \
		[(2*m-1, ~(3*m-1), m)]
		)
	colouring = dict([(i, RED) for i in range(2*m)] + [(i, BLUE) for i in range(2*m, 3*m)])
	return ColouredTriangulation(T, colouring)

def build(T, skip_dimension=False):
	
	print('Computing stratum: %s' % T.stratum())
	print('Full dimension: %d' % T.full_dimension())
	print('Good | Bad | current | ratio (current/good):')
	
	current = Queue()
	current.put(T)
	count = 0
	if not skip_dimension: assert(T.is_full_dimensional())
	good = set([T.iso_sig(T.good_starts())])
	bad = set()
	while not current.empty():
		count += 1
		T = current.get()
		neighbours = [T.flip_edge(i, c) for i in T.flippable_edges() for c in [RED, BLUE]]
		for n in neighbours:
			s = n.iso_sig(n.good_starts())
			if s not in good and s not in bad:
				if skip_dimension or n.is_full_dimensional():
					good.add(s)
					current.put(n)
				else:
					bad.add(s)
		if count % 1 == 0:
			print('\r%d %d %d %0.3f               ' % (len(good), len(bad), current.qsize(), float(current.qsize()) / len(good)), end='')
			sys.stdout.flush()
	print('')
	print(len(good))

if __name__ == '__main__':
	# test()
	
	T = ColouredTriangulation.from_pA(flipper.load('S_1_1').mapping_class('aB'))  # [0]  # 2.
	T = ColouredTriangulation.from_pA(flipper.load('S_1_2').mapping_class('abC'))  # [1, 1, -1, -1] # 8797 in 1m47s.
	# T = ngon(6)  # [0, 0] # 18 in 1s.
	# T = ngon(8)  # [4] # 120 in 3s.
	T = ngon(10)  # [2, 2] # 2062 in 1m4s.
	# T = ngon(12)  # [8] # Was 59342 in 52m21s. Now 9116 in 17m8s.
	# T = ngon(14)  # [4, 4] # 
	
	build(T)

