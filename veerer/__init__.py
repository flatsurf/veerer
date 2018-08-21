r"""
Veerer

A Python module to deal with Veering triangulation and L^infinity
Delaunay decomposition of surfaces.
"""
from __future__ import absolute_import

from .constants import RED, BLUE, HORIZONTAL, VERTICAL
from .triangulation import Triangulation
from .coloured_triangulation import ColouredTriangulation
from .automaton import Automaton

try:
    import sage.all
except ImportError:
    HAS_SAGE = False
else:
    HAS_SAGE = True

if HAS_SAGE:
    from .layout import FlatTriangulation

del absolute_import
