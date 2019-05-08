r"""
Various constants defined for veering triangulations.

EXAMPLES::

    sage: from veerer.constants import RED, BLUE, PURPLE, GREEN
"""

NONE = 0

# edge colors
RED    = 1<<0   # positive slope
BLUE   = 1<<1   # negative slope
PURPLE = 1<<2   # horizontal
GREEN  = 1<<3   # vertical

cRED    = 'R'
cBLUE   = 'B'
cPURPLE = 'P'
cGREEN  = 'G'

COLOURS = [RED, BLUE, PURPLE, GREEN]

# additional veering triangulation properties
SQUARETILED    = 1 << 2
QUADRANGULABLE = 1 << 3
CYLINDRICAL    = RED | BLUE
GEOMETRIC      = 1 << 4

ST = SQUARETILED | GEOMETRIC

PROPERTIES_COLOURS = {
    NONE : '#FFFFFF',  # white
    RED  : '#8B0000',               # dark red
    RED | ST : '#FF0000', # red
    BLUE : '#00008B',              # dark blue
    BLUE | ST : '#0000FF',# blue
    GEOMETRIC       : '#00AA00',   # green
    GEOMETRIC | RED : '#FAAA00',
    GEOMETRIC | BLUE: '#00AAFF'
    }

def key_property(p):
    return -((1<<4) * bool(p & SQUARETILED) | \
           (1<<3) * bool(p & GEOMETRIC) | \
           (1<<2) * bool(p & QUADRANGULABLE) | \
           (1<<1) * bool(p & RED) | \
           (1<<0) * bool(p & BLUE))

def properties_to_string(p):
    r"""
    EXAMPLES::

        sage: from veerer import VeeringTriangulation
        sage: from veerer.constants import *
        sage: T = VeeringTriangulation("(0,1,8)(2,~7,~1)(3,~0,~2)(4,~5,~3)(5,6,~4)(7,~8,~6)", "BRRRRBRBR")
        sage: properties_to_string(T.properties_code())
        'red geometric'

        sage: properties_to_string(ST|BLUE)
        'blue square-tiled'
        sage: properties_to_string(ST|RED)
        'red square-tiled'
    """
    s = []

    if p & RED and p & BLUE:
        raise ValueError("bicolored!")

    if p & SQUARETILED:
        if p & ST != ST:
            raise ValueError("square-tiled implies GQHV")
        if p & RED:
            return 'red square-tiled'
        elif p & BLUE:
            return 'blue square-tiled'
        else:
            raise ValueError("square-tiled must be cylindrical")
    elif p & RED:
        s.append('red')
    elif p & BLUE:
        s.append('blue')
    elif p & QUADRANGULABLE:
        s.append('quadrangulable')

    if p & QUADRANGULABLE and p & CYLINDRICAL and not p & SQUARETILED:
        raise RuntimeError

    if p & GEOMETRIC:
        s.append('geometric')

    if s:
        return ' '.join(s)
    else:
        return 'none'

# slopes
HORIZONTAL = 1
VERTICAL   = 2

SLOPES = [HORIZONTAL, VERTICAL]

def colour_from_char(c):
    r"""
    EXAMPLES::

        sage: from veerer.constants import colour_from_char, RED, BLUE, GREEN, PURPLE

        sage: assert colour_from_char('R') == RED
        sage: assert colour_from_char('B') == BLUE
        sage: assert colour_from_char('G') == GREEN
        sage: assert colour_from_char('P') == PURPLE

        sage: colour_from_char('X')
        Traceback (most recent call last):
        ...
        ValueError: unknown color 'X'
    """
    if c == cRED:
        return RED
    elif c == cBLUE:
        return BLUE
    elif c == cPURPLE:
        return PURPLE
    elif c == cGREEN:
        return GREEN
    else:
        raise ValueError("unknown color '%s'" % c)

def colour_to_char(col):
    r"""
    EXAMPLES::

        sage: from veerer.constants import colour_to_char, RED, BLUE, GREEN, PURPLE

        sage: assert colour_to_char(RED) == 'R'
        sage: assert colour_to_char(BLUE) == 'B'
        sage: assert colour_to_char(GREEN) == 'G'
        sage: assert colour_to_char(PURPLE) == 'P'
    """
    if col == RED:
        return 'R'
    elif col == BLUE:
        return 'B'
    elif col == PURPLE:
        return 'P'
    elif col == GREEN:
        return 'G'
    else:
        raise ValueError("unknown color %s" % col)



