from pyparma import Polyhedron
from pyparma.utils import fractionize, intize

# Build a matrix [b | A] of fractions.
A = [
	# Define positive orthant:
	[0,1,0],
	[0,0,1],
	# Define planes:
	[0,-1,0]
	]
P = Polyhedron(hrep=intize(A))  # Build Polyhedron corresponding to A.x + b >= 0.
print(P.poly.affine_dimension())  # Print dimension of Polyhedron.
