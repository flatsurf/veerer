r"""
Various environment modules.

TESTS::

    sage: from veerer.env import sage, flipper, surface_dynamics, ppl
"""

try:
    import sage.all
    import sage
except ImportError:
    sage = None

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

error_msg = {
    'curver': 'the function {} can only be called when the package curver is installed.',

    'sage': 'the function {} can only be called when running inside of Sage. See http://www.sagemath.org/',

    'surface_dynamics': 'the function {} only works when the package surface_dynamics is installed. See https://pypi.org/project/surface_dynamics/ for instructions.',

    'flipper': 'the function {} only works when the package flipper is installed. See https://pypi.org/project/flipper/ for instructions',

    'pplpy': 'the function {} only works when the package pplpy is installed. See https://pypi.org/project/pplpy/ for instructions.'
    }

missing_mods = {
    'curver': curver is None,
    'sage': sage is None,
    'flipper': flipper is None,
    'ppl': ppl is None,
    'surface_dynamics': surface_dynamics is None
    }

# TODO: use the traceback to find out who called this function!
# https://docs.python.org/2/library/traceback.html#traceback-examples
def require_package(mod_name, caller):
    if missing_mods[mod_name]:
        raise ValueError(error_msg[mod_name].format(caller))
