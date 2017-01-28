
from __future__ import print_function

from random import randint
from itertools import product

from colouredtriangulation import ColouredTriangulation, ngon
from automaton import Automaton
from constants import *

import numpy

import flipper

class Path(object):
	def __init__(self, zeta, start):
		self.zeta = zeta
		self.start = start
		self.matrices = [numpy.identity(self.zeta)]
		self.sigs = [self.start]
		self.sequence = []

class AutomataExplorer(object):
	def __init__(self, automaton, start, max_value, skip_earlier=True):
		self.automaton = automaton
		self.start = start
		self.max_value = max_value
		self.skip_earlier = skip_earlier
		
		self.zeta = self.automaton.zeta
		self.index = self.automaton.index[self.start]
	
	def all_paths(self, verbose=False):
		P = Path(self.zeta, self.start)
		while True:
			assert(len(P.sequence) == len(P.matrices) - 1)
			assert(len(P.sequence) == len(P.sigs) - 1)
			if verbose:
				if randint(0, 100) == 0:
					print('\r[' + ', '.join('%d/%d' % (i+1, len(self.automaton.graph[sig])) for i, sig in zip(P.sequence, P.sigs)) + ']' + ' '*10, end='')
				# print('[' + ', '.join('%d/%d' % (i, len(self.automaton.graph[sig])) for i, sig in zip(P.sequence, P.sigs)) + ']' + ' '*10)
				# print([self.automaton.index[sig] for sig in P.sigs])
			
			if self.automaton.index[P.sigs[-1]] >= self.index and numpy.sum(P.matrices[-1]) < self.max_value:  # Deeper.
				P.sequence.append(-1)
			else:
				P.matrices.pop()
				P.sigs.pop()
			
			new_sig = None
			while new_sig is None or self.automaton.index[new_sig] < self.index:
				while P.sequence[-1] == len(self.automaton.graph[P.sigs[-1]]) - 1:
					P.sequence.pop()
					P.matrices.pop()
					P.sigs.pop()
					if len(P.sequence) == 0:
						return
				
				P.sequence[-1] += 1
				new_sig, trans_matrix = self.automaton.graph[P.sigs[-1]][P.sequence[-1]]
			
			P.sigs.append(new_sig)
			P.matrices.append(trans_matrix * P.matrices[-1])
			yield P
		
		return
	
	def all_loops(self, verbose=False):
		for P in self.all_paths(verbose):
			if P.sigs[-1] == self.start:
				yield P
		return

if __name__ == '__main__':
	c = 1000
	# T = ColouredTriangulation.from_pA(flipper.load('SB_4').mapping_class('s_0S_1'))  # [-1, -1, -1, -1] # 6 in 1.3s.
	# T = ColouredTriangulation.from_pA(flipper.load('S_1_2').mapping_class('abC'))  # [1, 1, -1, -1] # 8797 in 1m47s Now 1658.
	# T = ngon(8)  # [4] # 86.
	for c in product(COLOURS, repeat=12):
		try:
			T = ColouredTriangulation(flipper.create_triangulation([(~11, 4, 7), (~10, ~7, 6), (~9, ~6, 5), (~8, ~5, ~0), (~4, 3, 11), (~3, 2, 10), (~2, 1, 9), (~1, 0, 8)]), dict(enumerate(c)), sanity=True)
			print(T.stratum())
			if T.stratum() == [3, 1]:
				print(T)
				break
		except AssertionError:
			pass
	# T = ColouredTriangulation.from_pA(flipper.load('S_2_1').mapping_class('a.a.a.a.B.c.d'))  # [2, 1, 1] # 16896.
	A = Automaton.from_triangulation(T, verbose=True)
	print(A)
	for i in reversed(range(len(A))):
		print('')
		print('='*60)
		print('Testing: %d' % i)
		print('='*60)
		AE = AutomataExplorer(A, A.iso_sigs[i], 40)
		for P in AE.all_loops(verbose=True):
			M = P.matrices[-1]
			pf_eigenvalue = max(numpy.linalg.eigvals(M))
			if pf_eigenvalue > 1 and pf_eigenvalue < c:
				c = pf_eigenvalue
				print('')
				print(c)
				print(M)
	
