r"""
Veerer

A Python module to deal with Veering triangulation and L^infinity
Delaunay decomposition of surfaces.
"""
from __future__ import absolute_import

from .constants import RED, BLUE, PURPLE, GREEN, HORIZONTAL, VERTICAL, RIGHT, LEFT, UP, DOWN
from .triangulation import Triangulation
from .cover import TriangulationCover
from .veering_triangulation import VeeringTriangulation, VeeringTriangulations
from .automaton import CoreAutomaton
from .flip_sequence import VeeringFlipSequence

from .env import sage
if sage is not None:
    from .flat_structure import FlatVeeringTriangulation
    from .layout import FlatVeeringTriangulationLayout
    from .measured_train_track import MeasuredTrainTrack

del sage, absolute_import
