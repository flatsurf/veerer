r"""
Veerer

A Python module to deal with Veering triangulation and L^infinity
Delaunay decomposition of surfaces.
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2023 Vincent Delecroix
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************

from .constants import RED, BLUE, PURPLE, GREEN, HORIZONTAL, VERTICAL, RIGHT, LEFT, UP, DOWN
from .triangulation import Triangulation
from .cover import TriangulationCover
from .veering_triangulation import VeeringTriangulation, VeeringTriangulations
from .linear_family import VeeringTriangulationLinearFamily
from .automaton import FlipGraph, CoreAutomaton, ReducedCoreAutomaton, GeometricAutomaton, GeometricAutomatonSubspace
from .flip_sequence import VeeringFlipSequence
from .flat_structure import FlatVeeringTriangulation
from .layout import FlatVeeringTriangulationLayout
from .measured_train_track import MeasuredTrainTrack
