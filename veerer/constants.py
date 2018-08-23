r"""
Various constants defined for veering triangulations.

EXAMPLES::

    sage: from veerer.constants import RED, BLUE, PURPLE, GREEN
"""

RED    = 1
BLUE   = 2
PURPLE = 4
GREEN  = 8

cRED = 'R'
cBLUE = 'B'
cPURPLE = 'P'
cGREEN = 'G'

COLOURS = [RED, BLUE]

def colour_from_char(c):
    r"""
    EXAMPLES::

        sage: from veerer.constants import colour_from_char, RED, BLUE, GREEN, PURPLE

        sage: colour_from_char('R') == RED
        True

        sage: colour_from_char('X')
        Traceback (most recent call last):
        ...
        ValueError: unknown color 'X'
    """
    if c == cRED:
        a[i] = RED
    elif c == cBLUE:
        a[i] = BLUE
    else:
        raise ValueError("unknown color '%s'" % c)

def colour_to_char(col):
    if col == RED:
        return 'R'
    elif col == BLUE:
        return 'B'
    else:
        raise ValueError("unknown color code '%d'" % d)

HORIZONTAL, VERTICAL = 'Horizontal', 'Vertical'

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

