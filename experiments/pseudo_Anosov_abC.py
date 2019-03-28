r"""
The pseudo-Anosov abC on the torus.

The stratum is Q(1,1-,1,-1) but sadly we do not have the folded edges
by default. Also, this example has a symmetry that makes it live on
the sphere. We should figure out the quotient.
"""
import flipper
from veerer.layout import FlatVeeringTriangulationLayout

S = flipper.load("S_1_2")
f = S.mapping_class("abC")

V = FlatVeeringTriangulationLayout.from_pseudo_anosov(f)

V.flip(4)
V.flip(6)
V.flip(10)
V.flip(3)
V.flip(11)
V.flip(0)
