r"""
Tests for optional packages used by veerer.
"""

from sage.features import PythonModule

flatsurf_feature = PythonModule(
    "flatsurf", url="https://flatsurf.github.io/sage-flatsurf/#installation"
)
surface_dynamics_feature = PythonModule(
    "surface_dynamics", url="https://flatsurf.github.io/surface-dynamics/#installation"
)
flipper_feature = PythonModule(
    "veerer", url="https://flatsurf.github.io/veerer/installation.html"
)
curver_feature = PythonModule(
    "curver", url="https://curver.readthedocs.io/en/master/user/install.html"
)
pynormaliz_feature = PythonModule(
    "PyNormaliz", url="https://github.com/Normaliz/PyNormaliz"
)
