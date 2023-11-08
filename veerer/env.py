r"""
Various environment modules.

TESTS::

    sage: from veerer.env import sage, flipper, surface_dynamics, ppl
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

CHECK = False

try:
    import sage.all
    import sage
except ImportError:
    sage = None
    import random
    op_LT = 0
    op_LE = 1
    op_EQ = 2
    op_NE = 3
    op_GT = 4
    op_GE = 5
    def rich_to_bool(op, c):
        if op == op_EQ:
            return c == 0
        elif op == op_NE:
            return c != 0
        elif op == op_LT:
            return c < 0
        elif op == op_LE:
            return c <= 0
        elif op == op_GE:
            return c >= 0
        elif op == op_GT:
            return c > 0
        else:
            raise RuntimeError
else:
    import sage.misc.prandom as random
    from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE, rich_to_bool


try:
    import surface_dynamics
except ImportError:
    surface_dynamics = None


try:
    import flipper
except ImportError:
    flipper = None

try:
    import curver
except ImportError:
    curver = None

try:
    import ppl
except ImportError:
    ppl = None


try:
    import PyNormaliz
except ImportError:
    PyNormaliz = None

error_msg = {
    'curver': 'the function {} can only be called when the package curver is installed.',

    'sage': 'the function {} can only be called when running inside of Sage. See https://www.sagemath.org/',

    'surface_dynamics': 'the function {} only works when the package surface_dynamics is installed. See https://pypi.org/project/surface_dynamics/ for instructions.',

    'flipper': 'the function {} only works when the package flipper is installed. See https://pypi.org/project/flipper/ for instructions',

    'ppl': 'the function {} only works when the package pplpy is installed. See https://pypi.org/project/pplpy/ for instructions.',

    'PyNormaliz': 'the function {} only works when the package PyNormaliz is installed.'
    }

missing_mods = {
    'curver': curver is None,
    'sage': sage is None,
    'flipper': flipper is None,
    'ppl': ppl is None,
    'surface_dynamics': surface_dynamics is None,
    'PyNormaliz': PyNormaliz is None
    }

# TODO: use the traceback to find out who called this function!
# https://docs.python.org/2/library/traceback.html#traceback-examples
def require_package(mod_name, caller):
    if missing_mods[mod_name]:
        raise ValueError(error_msg[mod_name].format(caller))
