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
    import ppl
except ImportError:
    ppl = None

msg = {
    'sage': 'the function {} can only be called when running inside of Sage. See http://www.sagemath.org/',

    'surface_dynamics': 'the function {} only works when the package surface_dynamics is installed. See https://pypi.org/project/surface_dynamics/ for instructions.',

    'flipper': 'the function {} only works when the package flipper is installed. See https://pypi.org/project/flipper/ for instructions',

    'pplpy': 'the function {} only works when the package pplpy is installed. See https://pypi.org/project/pplpy/ for instructions.'
    }

mods = {
    'sage': sage is not None,
    'flipper': flipper is not None,
    'ppl': ppl is not None,
    'surface_dynamics': surface_dynamics is not None
    }

def require_package(mod_name, caller):
    if mods[mod_name] is None:
        raise ValueError(error_msg.format(caller))
