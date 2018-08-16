
try:
    from pyparma import Polyhedron
    from pyparma.utils import intize
except ImportError:
    print('PyParma unavailable.')
    pass

class Polyhedron(object):
    
