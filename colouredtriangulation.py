
from __future__ import print_function

from Queue import Queue

from constants import *

try:
	from pyparma import Polyhedron
	from pyparma.utils import intize
except ImportError:
	print('PyParma unavailable.')
	pass

import flipper
norm = flipper.norm

def best_rotation(X):
	return min(X[i:] + X[:i] for i in range(len(X)))

class ColouredTriangulation(object):
	def __init__(self, triangulation, colouring, sanity=False):
		self.triangulation = triangulation
		self.colouring = colouring  # A dict : edge_indices --> {Red, Blue}
		self.zeta = self.triangulation.zeta
		if sanity:
			assert(all(colour in COLOURS for colour in self.colouring))
			assert(len(self.colouring) == self.zeta)
			assert(all(0 < [self.colouring[i] for i in t.indices].count(RED) < 3 for t in self.triangulation))
	
	@classmethod
	def from_pA(self, h):
		F = h.flat_structure()
		T = F.triangulation
		colouring = [RED if (F.edge_vectors[T.edge_lookup[i]].x > 0) == (F.edge_vectors[T.edge_lookup[i]].y > 0) else BLUE for i in range(T.zeta)]
		return ColouredTriangulation(T, colouring, sanity=True)
	
	@classmethod
	def from_QD(self, QD):
		return NotImplemented
	
	@classmethod
	def from_iso_sig(self, sig):
		# sig == (Edge colourings | triangles)
		# Hence len(sig) = zeta + F(S) = 5 / 3 zeta
		zeta = 3 * len(sig) // 5
		T = flipper.create_triangulation(sig[zeta:])
		colouring = sig[:zeta]
		return ColouredTriangulation(T, colouring, sanity=True)
	
	def __str__(self):
		return str(self.triangulation) + ', ' + str(self.colouring)
	def __repr__(self):
		return str(self)
	
	def is_abelian(self):
		if any(d % 2 == 1 for d in self.stratum()): return False
		# Perform BFT to check.
		
		oris = [None] * self.zeta
		oris[0] = True
		to_process = Queue()
		to_process.put(0)
		to_process.put(~0)
		while not to_process.empty():
			to_do = to_process.get()
			corner = self.triangulation.corner_of_edge(to_do)
			rev_oris = [corner.indices[i] != corner.labels[i] for i in range(3)]
			colours = [self.colouring[i] for i in corner.indices]
			
			if colours.count(RED) == 2:  # colours.count(RED) > colours.count(BLUE).
				new_ori1 = oris[corner.indices[0]] ^ (colours[0] == RED) ^ rev_oris[0] ^ rev_oris[1]
				new_ori2 = oris[corner.indices[0]] ^ (colours[2] == RED) ^ rev_oris[0] ^ rev_oris[2]
			else:  # #BLUE < # RED.
				new_ori1 = oris[corner.indices[0]] ^ (colours[1] == BLUE) ^ rev_oris[0] ^ rev_oris[1]
				new_ori2 = oris[corner.indices[0]] ^ (colours[0] == BLUE) ^ rev_oris[0] ^ rev_oris[2]
			
			for index, new_ori in enumerate([new_ori1, new_ori2], start=1):
				if oris[corner.indices[index]] is None:
					oris[corner.indices[index]] = new_ori
					to_process.put(~corner.labels[index])
				elif oris[corner.indices[index]] != new_ori:
					return False
		
		return True
	
	def stratum_dimension(self):
		return 2*self.triangulation.genus - 2 + self.triangulation.num_vertices + (1 if self.is_abelian() else 0)
	
	def colours_about_edge(self, i):
		return [self.colouring[e.index] for e in self.triangulation.square_about_edge(i)]
	
	def alternating_square(self, i):
		colours = self.colours_about_edge(i)
		return all(colours[j] != colours[(j+1) % 4] for j in range(4))
	
	def is_flippable(self, i):
		return self.triangulation.is_flippable(i) and self.alternating_square(i)
	
	def flippable_edges(self):
		return [i for i in self.triangulation.flippable_edges() if self.is_flippable(i)]
	
	def mostly_sloped_edges(self, slope):
		return [i for i in self.flippable_edges() if self.colouring[self.triangulation.corner_lookup[i].indices[1]] == (BLUE if slope == VERTICAL else RED)]
	
	def is_isomorphic_to(self, other):
		return self.iso_sig() == other.iso_sig()
	
	def isometries_to(self, other):
		return [isom for isom in self.triangulation.isometries_to(other.triangulation) if all(self.colouring[i] == other.colouring[isom.index_map[i]] for i in range(self.zeta))]
	
	def flip_edge(self, i, colour):
		assert(self.is_flippable(i))
		
		T = self.triangulation.flip_edge(i)
		colouring = list(self.colouring)
		colouring[norm(i)] = colour
		return ColouredTriangulation(T, colouring)
	
	def train_track_matrix(self, slope=VERTICAL):
		M = []
		for t in self.triangulation:
			corner = [corner for corner in t.corners if self.colouring[corner.indices[1]] == self.colouring[corner.indices[2]]][0]
			I = corner.indices
			corner_colour = self.colouring[I[1]]
			if (slope == VERTICAL) != (corner_colour == RED):
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
	
	def train_track_dimension(self, slope=VERTICAL):
		# Uses PyParma.
		M = self.train_track_matrix(slope)
		A = [[0] + [0] * i + [1] + [0] * (self.zeta - i - 1) for i in range(self.zeta)]
		for row in M:
			A.append([0] + [i for i in row])
			A.append([-0] + [-i for i in row])
		
		return Polyhedron(hrep=intize(A)).poly.affine_dimension()
	
	def geometric_matrix(self):
		G = []
		for i in self.flippable_edges():
			corner1 = self.triangulation.corner_lookup[i]
			corner2 = self.triangulation.corner_lookup[~i]
			
			if self.colouring[corner1.indices[1]] == BLUE:  # Mostly horizontal.
				row = [-1 if j == i else 0 for j in range(self.zeta)] + [1 if j == corner1.indices[1] or j == corner2.indices[2] else 0 for j in range(self.zeta)]
			else:  # Mostly vertical.
				row = [1 if j == corner1.indices[1] or j == corner2.indices[2] else 0  for j in range(self.zeta)] + [-1 if j == i else 0 for j in range(self.zeta)]
			G.append(row)
		return G
	
	def geometric_dimension(self):
		# Uses PyParma.
		V = self.train_track_matrix(VERTICAL)
		H = self.train_track_matrix(HORIZONTAL)
		G = self.geometric_matrix()
		
		A = [[0] + [0] * i + [1] + [0] * (2*self.zeta - i - 1) for i in range(2*self.zeta)]
		for row in V:
			A.append([0] + [i for i in row] + [0] * self.zeta)
			A.append([-0] + [-i for i in row] + [-0] * self.zeta)
		for row in V:
			A.append([0] + [0] * self.zeta + [i for i in row])
			A.append([-0] + [-0] * self.zeta + [-i for i in row])
		for row in G:
			A.append([0] + row)
		
		return Polyhedron(hrep=intize(A)).poly.affine_dimension()
	
	def sub_train_track_dimensions(self, slope=VERTICAL):
		# Uses PyParma.
		M = self.train_track_matrix(slope)
		for j in range(self.zeta):
			A = [[0] + [0] * i + [1] + [0] * (self.zeta - i - 1) for i in range(self.zeta)]
			for row in M:
				A.append([0] + [i for i in row])
				A.append([-0] + [-i for i in row])
			A.append([-0] + [0] * j + [-1] + [0] * (self.zeta - j - 1))
			yield Polyhedron(hrep=intize(A)).poly.affine_dimension()
	
	def vertex_data_dict(self):
		return dict((CC[0].vertex, (len([c for c in CC if self.colouring[c.indices[1]] == BLUE and self.colouring[c.indices[2]] == RED]) - 2, len(CC))) for CC in self.triangulation.corner_classes)
	
	def stratum(self):
		return sorted([d for d, v in self.vertex_data_dict().values()], reverse=True)
	
	def is_full_dimension(self):
		d = self.stratum_dimension()
		return all(self.train_track_dimension(slope) == d for slope in SLOPES)  # and self.geometric_dimension() >= 2*d - 1
	
	def good_starts(self):
		VD = self.vertex_data_dict()
		
		best = min(VD.values())
		def test(edge):
			if VD[edge.source_vertex] != best: return False
			if self.colouring[edge.index] != RED: return False
			if self.colouring[self.triangulation.corner_lookup[edge.label].indices[2]] != BLUE: return False
			return True
		
		return [edge.label for edge in self.triangulation.edges if test(edge)]
	
	def best_translation(self):
		best, best_translation = None, None
		start_edges = self.good_starts()
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
				best_translation = translate
		
		return best, best_translation
	
	def iso_sig(self):
		return self.best_translation()[0]
	
	def canonical(self):
		_, translate = self.best_translation()
		inv_translate = dict((v, k) for (k, v) in translate.items())
		
		T = flipper.create_triangulation([[translate[edge.label] if edge.label > 0 else ~translate[~edge.label] for edge in triangle] for triangle in self.triangulation])
		colouring = [self.colouring[norm(inv_translate[i])] for i in range(self.zeta)]
		return ColouredTriangulation(T, colouring)

def ngon(n):
	assert(n > 4 and n % 2 == 0)
	m = (n - 2) // 2
	T = flipper.create_triangulation( \
		[(i, 2*m+i, ~(i+1)) for i in range(m)] + \
		[(~0, ~(2*m), ~(m+1))] + \
		[(m+i, ~(2*m+i), ~(m+i+1)) for i in range(1, m-1)] + \
		[(2*m-1, ~(3*m-1), m)]
		)
	colouring = [RED] * (2*m) + [BLUE] * m
	return ColouredTriangulation(T, colouring)

def test():
	T = ColouredTriangulation(flipper.create_triangulation([(0,1,2), (~0,~1,~2)]), [RED, BLUE, RED])
	
	print(T)
	print(T.stratum())
	print(T.iso_sig())
	T2 = ColouredTriangulation.from_iso_sig(T.iso_sig())
	print(T2.iso_sig())
	print(T.geometric_dimension())
	
	print(T.flippable_edges())
	T2 = T.flip_edge(0, BLUE)
	print(T2)
	print(T2.stratum())
	print(T2.iso_sig())
	
	print(T2.canonical())

if __name__ == '__main__':
	test()

