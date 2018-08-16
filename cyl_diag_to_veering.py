import flipper
from constants import *

def to_veering_triangulation(c):
    r"""
    Return a Veering triangulation from a cylinder decomposition.

    The input can be a cylinder diagram of either an orientable or
    non-orientable differential.

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: c = QuadraticCylinderDiagram('(0)-(0)')
        sage: to_veering_triangulation(c)
        ([(~2, 1, ~0), (~1, 0, 2)], 'RBB')

        sage: c = QuadraticCylinderDiagram('(0,0)-(1,1,2,2,3,3)')
        sage: to_veering_triangulation(c)
        ([(~11, 4, ~3), (~10, ~0, 11), (~9, 0, 10), (~8, 9, 1),
          (~7, 8, ~1), (~6, 7, 2), (~5, 6, ~2), (~4, 5, 3)],
          'RRRRBBBBBBBB')

    In order to build a cylinder diagram, you need the development
    version of ``surface_dynamics`` and do::

        sage: Q = QuadraticStratum(9,-1)
        sage: c1,c2 = Q.components()
        sage: c1.one_cylinder_diagram()
        (0,1,2,1,2)-(0,3,4,3,4,5,5)
        sage: c2.one_cylinder_diagram()
        (0,0,1)-(1,2,3,4,5,2,3,4,5)

    """
    unseen = [False] * c.nseps()
    triangles = []
    k = c.nseps()
    for bot,top in c.cylinders():
        # len(bot) + len(top)
        l = k + len(top) - 1
        for i in bot:
            assert isinstance(i,int)
            if unseen[i]:
                i = ~i
            else:
                unseen[i] = True
            triangles.append((i,l+1,~l))
            l += 1
        l = k + len(top) - 1
        for i in top:
            assert isinstance(i,int)
            if unseen[i]:
                i = ~i
            else:
                unseen[i] = True
            triangles.append((i,~(l-1),l))
            l -= 1
        # last one got wrong
        i, j, k = triangles.pop()
        triangles.append((i,~(k+len(bot)+len(top)-1),k))

    colors = [RED] * c.nseps() + [BLUE] * (2*c.nseps())

    print triangles
    import flipper
    return flipper.kernel.Triangulation.from_tuple(triangles), colors



