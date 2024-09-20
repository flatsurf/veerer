r"""
Tests for optional packages used by veerer.
"""

from sage.features import PythonModule

sage_flatsurf_feature = PythonModule(
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
pyflatsurf_feature = PythonModule(
    "pyflatsurf", url="https://github.com/flatsurf/flatsurf"
)
pynormaliz_feature = PythonModule(
    "PyNormaliz", url="https://github.com/Normaliz/PyNormaliz"
)
