r"""
Veerer

A Python module to deal with Veering triangulation and L^infinity
Delaunay decomposition of surfaces.
"""
from __future__ import absolute_import

from .constants import RED, BLUE, HORIZONTAL, VERTICAL
from .triangulation import Triangulation
from .veering_triangulation import VeeringTriangulation
from .automaton import Automaton

from .env import sage
if sage is not None:
    from .layout import FlatVeeringTriangulationLayout

del sage, absolute_import
