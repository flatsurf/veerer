r"""
Veerer

A Python module to deal with Veering triangulation and L^infinity
Delaunay decomposition of surfaces.
"""
from __future__ import absolute_import

from .constants import RED, BLUE, HORIZONTAL, VERTICAL
from .triangulation import Triangulation
from .cover import TriangulationCover
from .veering_triangulation import VeeringTriangulation, VeeringTriangulations
from .measured_train_track import MeasuredTrainTrack
from .automaton import Automaton

from .env import sage
if sage is not None:
    from .layout import FlatVeeringTriangulationLayout

del sage, absolute_import
