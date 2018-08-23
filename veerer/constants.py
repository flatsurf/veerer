
import string

RED = 'Red'
BLUE = 'Blue'
PURPLE = 'Purple'
GREEN = 'Green'

HORIZONTAL, VERTICAL = 'Horizontal', 'Vertical'

COLOURS = [RED, BLUE]
SLOPES = [HORIZONTAL, VERTICAL]

CORE = 0
GEOMETRIC = 1
CYLINDRICAL = 1 << 2

TYPE_COLOURS = {
    CORE: '#BFBFBF',          # light gray
    GEOMETRIC: '#50AA50',     # green
    CYLINDRICAL: '#FFA500',   # orange
    GEOMETRIC | CYLINDRICAL: '#FA8072' # salmon
    }

