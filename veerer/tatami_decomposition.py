r"""
Quartic differentials built from rectangles.
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

from collections import defaultdict

from .constants import LEFT, RIGHT


def tatami_decomposition(rectangles, base_ring=None):
    r"""
    Return a sage-flatsurf surface built from the given ``rectangles``

    Each point is identified by its position on a separatrix. It is encoded in
    a triple ``(separatrix_label, side, distance_to_singularity)`` where
    - ``separatrix_label``: is a label identifying a separatrix
    - ``side``: specifies either left (``LEFT``) or right (``RIGHT``) of the separatrix
    - ``distance_to_singularity``: a positive real number

    The input ``rectangles`` must be a list of rectangles where each rectangle
    is given as a list of 8 points (encoded as above). These 8 points correspond
    to the position on the separatrices corresponding to the various sides as in
    the picture below::

          o----------------o
          |  p5         p4 |
          |p6            p3|
          |                |
          |p7            p2|
          | p0          p1 |
          o----------------o

    EXAMPLES::

        sage: from veerer.tatami_decomposition import tatami_decomposition  # random output due to deprecation warnings from realalg
        sage: from veerer.constants import LEFT, RIGHT
        sage: r0 = ((1, RIGHT, 1), (1, RIGHT, 0), (0, LEFT, 0), (0, LEFT, 1), (3, RIGHT, 2), (3, RIGHT, 1), (0, RIGHT, 2), (0, RIGHT, 1))
        sage: r1 = ((3, LEFT, 0), (3, LEFT, 2), (0, LEFT, 1), (0, LEFT, 2), (3, RIGHT, 1), (1, LEFT, 1), (0, RIGHT, 1), (0, RIGHT, 0))
        sage: tatami_decomposition([r0, r1])  # optional: sage_flatsurf
        Translation Surface built from a square and a rectangle

    TESTS::

        sage: tatami_decomposition([r0, r1], ZZ)  # optional: sage_flatsurf
        Translation Surface built from a square and a rectangle
        sage: tatami_decomposition([r0, r1], QQ)  # optional: sage_flatsurf
        Translation Surface built from a square and a rectangle
        sage: tatami_decomposition([r0, r1], AA)  # optional: sage_flatsurf
        Translation Surface built from a square and a rectangle

        sage: r0 = ((1, RIGHT, AA(1)), (1, RIGHT, 0), (0, LEFT, 0), (0, LEFT, 1), (3, RIGHT, 2), (3, RIGHT, 1), (0, RIGHT, 2), (0, RIGHT, 1))
        sage: r1 = ((3, LEFT, 0), (3, LEFT, 2), (0, LEFT, 1), (0, LEFT, 2), (3, RIGHT, 1), (1, LEFT, 1), (0, RIGHT, 1), (0, RIGHT, 0))
        sage: tatami_decomposition([r0, r1])  # optional: sage_flatsurf
        Translation Surface built from a square and a rectangle
    """
    if base_ring is None:
        coefficients = []
        for rectangle in rectangles:
            for sep, side, x in rectangle:
                coefficients.append(x)
        from sage.structure.sequence import Sequence
        base_ring = Sequence(coefficients).universe()

    # for each separatrix we build an iet
    zero = base_ring(0)
    intervals = defaultdict(lambda: defaultdict(list))
    for i, rectangle in enumerate(rectangles):
        if len(rectangle) != 8:
            raise ValueError('invalid data: each rectangle should be made of eight points on separatrices; got {}'.format(rectangle))
        lengths = []
        for i in range(0, 8, 2):
            (sep0, side0, x0) = rectangle[i]
            (sep1, side1, x1) = rectangle[i + 1]
            x0 = base_ring(x0)
            x1 = base_ring(x1)
            if x0 < 0:
                raise ValueError('invalid length data {}'.format(x0))
            if x1 < 0:
                raise ValueError('invalid length data {}'.format(x1))
            if sep0 == sep1:
                if side0 != side1:
                    raise ValueError('invalid data')
                if x1 < x0:
                    x0, x1 = x1, x0
                intervals[sep0][side0].append((x0, x1))
                lengths.append(x1 - x0)
            else:
                if side0 == side1:
                    raise ValueError('invalid data')
                intervals[sep0][side0].append((zero, x0))
                intervals[sep1][side1].append((zero, x1))
                lengths.append(x0 + x1)

        if lengths[0] == 0 or lengths[1] == 0 or lengths[0] != lengths[2] or lengths[1] != lengths[3]:
            raise ValueError('invalid data: each rectangle must have opposite sides of equal lengths. Got {} which have lengths {}'.format(rectangle, lengths))

    # refinement: for each interval (sep, side, x0, x1) we store additional cut points
    refinement = {}
    for sep in intervals:
        if len(intervals[sep]) != 2 or LEFT not in intervals[sep] or RIGHT not in intervals[sep]:
            raise ValueError('invalid data: for separatrix={} got intervals={}'.format(sep, intervals[sep]))
        left = intervals[sep][LEFT]
        right = intervals[sep][RIGHT]
        left.sort()
        right.sort()
        if (left[0][0] != 0 or right[0][0] != 0 or
            left[-1][1] != right[-1][1] or
             any(left[i][1] != left[i + 1][0] for i in range(len(left) - 1)) or
             any(right[i][1] != right[i + 1][0] for i in range(len(right) - 1))):
            raise ValueError('invalid data: for separatrix {} we got left={} right={}'.format(sep, left, right))

        ir = 0
        xr0, xr1 = right[0]
        for il in range(len(left)):
            xl0, xl1 = left[il]
            assert xr0 <= xl0 and xr1 >= xl0, (xl0, xl1, xr0, xr1)
            if xr1 == xl0:
                ir += 1
                xr0, xr1 = right[ir]
            cut = [xl0]
            while xl0 < xr1 < xl1:
                cut.append(xr1)
                ir += 1
                xr0, xr1 = right[ir]
            cut.append(xl1)
            assert all(cut[k] < cut[k+1] for k in range(len(cut) - 1))
            assert (sep, xl0, xl1) not in refinement
            refinement[sep, xl0, xl1] = cut

        il = 0
        xl0, xl1 = left[0]
        for ir in range(len(right)):
            xr0, xr1 = right[ir]
            assert xl0 <= xr0 and xl1 >= xr0, (xr0, xr1, xl0, xl1)
            cut = [xr0]
            if xl1 == xr0:
                il += 1
                xl0, xl1 = left[il]
            while xr0 < xl1 < xr1:
                cut.append(xl1)
                il += 1
                xl0, xl1 = left[il]
            cut.append(xr1)
            assert all(cut[k] < cut[k+1] for k in range(len(cut) - 1))
            if (sep, xr0, xr1) in refinement:
                assert refinement[sep, xr0, xr1] == [xr0, xr1]
            else:
                refinement[sep, xr0, xr1] = cut

    # rerun through rectangles and compute how each edge is cut into pieces
    # and turn this into proper sage-flatsurf rectangles recording edge indices
    from flatsurf.geometry.polygon import Polygon
    from flatsurf import MutableOrientedSimilaritySurface
    from sage.modules.free_module import FreeModule
    V = FreeModule(base_ring, 2)
    surface = MutableOrientedSimilaritySurface(base_ring)
    directions = [V((1, 0)), V((0, 1)), V((-1, 0)), V((0, -1))]
    interval_id = {}
    for i, rectangle in enumerate(rectangles):
        vertices = [V()]
        m = 0
        for j, direction in zip(range(0, 8, 2), directions):
            (sep0, side0, x0) = rectangle[j]
            (sep1, side1, x1) = rectangle[j + 1]
            x0 = base_ring(x0)
            x1 = base_ring(x1)
            if sep0 == sep1:
                if x1 < x0:
                    cuts = refinement[sep0, x1, x0]
                else:
                    cuts = refinement[sep0, x0, x1]
                lengths = [cuts[k + 1] - cuts[k] for k in range(len(cuts) - 1)]
                if x0 < x1:
                    for k in range(len(cuts) - 1):
                        interval_id[sep0, side0, cuts[k], cuts[k+1]] = (i, m)
                        m += 1
                else:
                    for k in range(len(cuts) - 2, -1, -1):
                        interval_id[sep0, side0, cuts[k], cuts[k + 1]] = (i, m)
                        m += 1
                    lengths.reverse()
            else:
                cuts = refinement[sep0, zero, x0]
                for k in range(len(cuts) - 2, -1, -1):
                    interval_id[sep0, side0, cuts[k], cuts[k + 1]] = (i, m)
                    m += 1
                lengths0 = [cuts[k + 1] - cuts[k] for k in range(len(cuts) - 1)]

                cuts = refinement[sep1, zero, x1]
                for k in range(len(cuts) - 1):
                    interval_id[sep1, side1, cuts[k], cuts[k + 1]] = (i, m)
                    m += 1
                lengths1 = [cuts[k + 1] - cuts[k] for k in range(len(cuts) - 1)]

                lengths = lengths0[::-1] + lengths1

            for l in lengths:
                vertices.append(vertices[-1] + l * direction)

        assert vertices[0] == vertices[-1], vertices
        vertices.pop()
        surface.add_polygon(Polygon(vertices=vertices, base_ring=base_ring))

    # gluing: we go through left sides
    for (sep, side, x, y) in interval_id:
        if side == LEFT:
            if (sep, RIGHT, x, y) not in interval_id:
                raise RuntimeError('ERROR: missing (sep={}, side=RIGHT, x={}, y={})'.format(sep, x, y))
            surface.glue(interval_id[sep, LEFT, x, y], interval_id[sep, RIGHT, x, y])
        elif side == RIGHT and (sep, LEFT, x, y) not in interval_id:
            raise RuntimeError('ERROR: missing (sep={}, side=LEFT, x={}, y={})'.format(sep, x, y))

    return surface
