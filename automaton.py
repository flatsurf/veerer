
from __future__ import print_function

import sys
from Queue import Queue
from numpy import matrix as Matrix
from math import log

from constants import *
from colouredtriangulation import ColouredTriangulation, ngon

import flipper
norm = flipper.norm

def compress_iso_sig(sig):
	d = len(CHARS)
	zeta = 3 * len(sig) // 5
	colours, tuples = sig[:zeta], sig[zeta:]
	
	if 2*zeta < d - 1:
		char_start = CHARS[zeta]
		digits = 1
	else:
		digits = int(log(2*zeta) / log(d)) + 1
		char_start = CHARS[-1] + CHARS[digits] + ''.join(CHARS[(2*zeta // d**i) % d] for i in range(digits))
	
	step = int(log(d) / log(2))
	colours = colours + tuple([BLUE] * step)  # Padding, just in case.
	char_colours = ''.join(CHARS[sum(2**i if colours[j+i] == RED else 0 for i in range(step))] for j in range(0, zeta, step))
	
	char_edges = ''.join(''.join(CHARS[(i+zeta // d**j) % d] for j in range(digits)) for t in tuples for i in t)
	
	return char_start + '.' + char_colours + '.' + char_edges

class Automaton(object):
	def __init__(self, graph, zeta):
		self.graph = graph
		self.zeta = zeta
		self.iso_sigs = sorted(graph)
		self.index = dict((sig, index) for index, sig in enumerate(self.iso_sigs))
	def __str__(self):
		return '\n'.join('%s -- {%s}' % (compress_iso_sig(g), ', '.join(compress_iso_sig(x) for x, _ in self.graph[g])) for g in self.graph)
	def __repr__(self):
		return str(self)
	def __len__(self):
		return len(self.iso_sigs)
	def __iter__(self):
		return iter(self.iso_sigs)
	def __contains__(self, item):
		return item in self.iso_sigs
	
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
		
		def is_full_dimension(t):
			if not all(t.train_track_dimension(slope) == d for slope in SLOPES): return False
			# if not t.geometric_dimension() >= 2*d - 1: return False
			if any(d2 == d for slope in SLOPES for d2 in t.sub_train_track_dimensions(slope)): return False
			return True
		
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
							current.put(neighbour.canonical())
						else:
							bad.add(s)
					if s in good:
						corner1 = T.triangulation.corner_lookup[i]
						corner2 = T.triangulation.corner_lookup[~i]
						M = Matrix([[1 if (j != i and k == j) or (j == i and k in [corner1.indices[1], corner2.indices[2]]) else 0 for k in range(zeta)] for j in range(zeta)])
						N = Matrix([[1 if j == norm(t[k]) else 0 for k in range(zeta)] for j in range(zeta)])
						
						neighbours.append((s, N * M))
			
			graph[T.iso_sig()] = sorted(neighbours, key=lambda (s, A): (s, A.tolist()))
			
			if verbose and count % 1 == 0:
				print('\r%d %d %d %0.3f               ' % (len(good), len(bad), current.qsize(), float(current.qsize()) / len(good)), end='')
				sys.stdout.flush()
		
		if verbose:
			print('\r%d %d %d %0.3f               ' % (len(good), len(bad), current.qsize(), float(current.qsize()) / len(good)))
		
		return Automaton(graph, zeta)

if __name__ == '__main__':
	T = ColouredTriangulation.from_pA(flipper.load('S_1_1').mapping_class('aB'))  # [0]  # 2.
	T = ColouredTriangulation.from_pA(flipper.load('S_1_2').mapping_class('abC'))  # [1, 1, -1, -1] # 8797 in 1m47s Now 1658.
	# T = ColouredTriangulation.from_pA(flipper.load('SB_4').mapping_class('s_0S_1'))  # [-1, -1, -1, -1] # 6 in 1.3s.
	# T = ColouredTriangulation(flipper.create_triangulation([(~11, ~10, ~8), (~9, ~4, 10), (~7, ~6, 11), (~5, 7, ~2), (~3, 8, 9), (~1, 5, 6), (~0, 3, 4), (0, 1, 2)]), [RED, BLUE, BLUE, BLUE, BLUE, RED, RED, RED, RED, RED, RED, BLUE])  # [3, 1] # 16.
	# T = ngon(6)  # [0, 0] # 16 in 1s.
	# T = ngon(8)  # [4] # 120 in 3s now 86 in 5s.
	# T = ngon(10)  # [2, 2] # 2062 in 1m4s now 876.
	# T = ColouredTriangulation.from_pA(flipper.load('S_2_1').mapping_class('a.a.a.a.B.c.d'))  # [2, 1, 1] # 16896.
	# T = ngon(12)  # [8] # Was 59342 in 52m21s. Now 9116 in 17m8s.
	# T = ngon(14)  # [4, 4] # 
	
	A = Automaton.from_triangulation(T, verbose=True)
	print(A)
	print(len(A))

