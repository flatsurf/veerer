r"""
Flat veering triangulations layout in the plane

This module is mostly intended for plotting features. A layout is defined by a
choice of a rooted spanning forest together with a coordinate for each root.
Each tree in the spanning forest corresponds to triangulated polygon.

Note:

- when graphviz generates a svg it specifies a given size with the
  attributes "width" and "height". This would better be redefined
  to width="100%".
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2018 Mark Bell
#                     2018-2023 Vincent Delecroix
#                     2018 Saul Schleimer
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

import math
import itertools

from sage.categories.fields import Fields

from sage.rings.all import RDF

from sage.modules.free_module import VectorSpace

from sage.misc.prandom import shuffle
from sage.plot.graphics import Graphics
from sage.plot.line import line2d
from sage.plot.polygon import polygon2d
from sage.plot.text import text
from sage.plot.bezier_path import bezier_path
from sage.plot.point import point2d

from .constants import BLUE, RED, PURPLE, GREEN, HORIZONTAL, VERTICAL, RIGHT
from .permutation import perm_init, perm_check, perm_on_list
from .misc import flipper_edge, flipper_edge_perm, flipper_nf_to_sage, flipper_nf_element_to_sage, det2, flipper_face_edge_perms
from .triangulation import Triangulation
from .veering_triangulation import VeeringTriangulation
from .flat_structure import FlatVeeringTriangulation
from .flat_structure import vec_slope

_Fields = Fields()

EDGE_COLORS = {
    BLUE: 'blue',
    RED: 'red',
    PURPLE: 'purple',
    GREEN: 'green'}

TIKZ_EDGE_COLORS = {
    BLUE: 'veeringblue',
    RED: 'veeringred',
    PURPLE: 'veeringpurple',
    GREEN: 'veeringgreen'}

TIKZ_FACE_COLORS = {
    BLUE: 'veeringblue!20',
    RED: 'veeringred!20'}


class FlatVeeringTriangulationLayout(object):
    r"""
    A flat triangulation layout in the plane.

    EXAMPLES:

        sage: from veerer import *  # random output due to deprecation warnings from realalg

    A can be built from a list of triangles and vectors::

        sage: triangles = [(0, 1, 2), (-1, -2, -3)]
        sage: holonomies = [(1, 2), (-2, -1), (1, -1), (1, -1), (-2, -1), (1, 2)]
        sage: fvt = FlatVeeringTriangulation(triangles, holonomies)
        sage: layout = FlatVeeringTriangulationLayout(fvt)
        sage: layout  # random
        FlatTriangulationLayout(FlatVeeringTriangulation(Triangulation("(0,1,2)(~2,~0,~1)"), [(1, 2), (-2, -1), (1, -1), (-1, 1), (2, 1), (-1, -2)]))
        sage: layout.set_pos()
        sage: layout  # random
        FlatTriangulationLayout(FlatVeeringTriangulation(Triangulation("(0,1,2)(~2,~0,~1)"), [(1, 2), (-2, -1), (1, -1), (-1, 1), (2, 1), (-1, -2)]), [(1, 0), (2, 2), (0, 1), (3, 1), (1, 0), (2, 2)])

    Or from a preexisting flat surface::

        sage: T = VeeringTriangulation("(0,1,2)(~0,~1,3)", "BRBB")
        sage: F = FlatVeeringTriangulation.from_coloured_triangulation(T)
        sage: F.layout()  # random
        FlatTriangulationLayout(FlatVeeringTriangulation(Triangulation("(0,1,2)(3,~0,~1)"), [(1, -2), (1, 1), (-2, 1), (2, -1), (-1, -1), (-1, 2)]))
    """
    def __init__(self, flat_triangulation, pos=None):
        r"""
        INPUT:

        - flat_triangulation - a flat triangulation
        """

        self._triangulation = FlatVeeringTriangulation(flat_triangulation, mutable=True)

        # for actual display.
        self._pos = pos     # vertex positions (list of length n)

        self._check()

    def _check(self):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,3)", "BRRR")
            sage: assert T.is_core()
            sage: F = T.flat_structure_min().layout()
            sage: F.set_pos()
            sage: F._check()
        """
        self._triangulation._check()

        n = self._triangulation.num_half_edges()
        ep = self._triangulation.edge_permutation(copy=False)
        vectors = self._triangulation._holonomies

        for a in range(n):
            A = ep[a]
            u = vectors[a]
            v = vectors[A]
            if u != v and u != -v:
                raise ValueError('ep[%s] = %s but vec[%s] = %s' % (a, u, A, v))

        for a, b, c in self._triangulation.faces():
            va = vectors[a]
            vb = vectors[b]
            vc = vectors[c]
            if va + vb + vc:
                raise ValueError('vec[%s] = %s, vec[%s] = %s and vec[%s] = %s do not sum to zero' % (a, va, b, vb, c, vc))

            if det2(va, vb) <= 0 or det2(vb, vc) <= 0 or det2(vc, va) <= 0:
                raise ValueError('(%s, %s, %s) is a clockwise triangle' %
                        (a, b, c))

        pos = self._pos
        if pos is not None:
            for a, b, c in self._triangulation.faces():
                if (pos[a] is not None and pos[b] is not None and pos[a] + vectors[a] != pos[b]) or \
                    (pos[b] is not None and pos[c] is not None and pos[b] + vectors[b] != pos[c]) or \
                        (pos[c] is not None and pos[a] is not None and pos[c] + vectors[c] != pos[a]):
                    raise ValueError('pos[%s] = %s, pos[%s] = %s, pos[%s] = %s while vec[%s] = %s, vec[%s] = %s, vec[%s] = %s' % (
                        a, pos[a],
                        b, pos[b],
                        c, pos[c],
                        a, vectors[a],
                        b, vectors[b],
                        c, vectors[c]))

    def copy(self):
        res = FlatVeeringTriangulationLayout.__new__(FlatVeeringTriangulationLayout)
        res._triangulation = self._triangulation.copy()
        if self._pos is not None:
            res._pos = [v.__copy__() for v in self._pos]
        else:
            res._pos = None
        return res

    def xy_scaling(self, sx, sy):
        r"""
        Apply Teichmueller flow to actually see things.
        """
        for i, v in enumerate(self._triangulation._holonomies):
            self._triangulation._holonomies[i] = self._triangulation._V((sx * v[0], sy * v[1]))

    def _edge_slope(self, e):
        return vec_slope(self._triangulation._holonomies[e])

    def __repr__(self):
        if self._pos is None:
            return 'FlatTriangulationLayout({})'.format(
                self._triangulation)
        else:
            return 'FlatTriangulationLayout({}, {})'.format(
                self._triangulation, self._pos)

    def _edge_is_boundary(self, e):
        r"""
        Test whether the edge ``e`` is on the boundary of the display.
        """
        if self._pos is None:
            return False
        fp = self._triangulation.face_permutation(copy=False)
        ep = self._triangulation.edge_permutation(copy=False)
        pos = self._pos
        vectors = self._triangulation._holonomies
        E = ep[e]
        return pos[fp[e]] is None or pos[E] is None or pos[fp[e]] != pos[E] or vectors[e] != -vectors[E]

    # there is something wrong with edge gluing
    # sometimes we end up with non-valid positions...
    def glue_edge(self, e):
        r"""
        Glue the triangle across the edge ``e`` to ``E`` so that
        we have a quadrilateral around ``e``.

        TESTS::

            sage: from veerer import *
            sage: t = "(0,~10,9)(1,~14,~2)(2,13,~3)(3,14,~4)(4,~13,~5)(5,~11,~6)(6,10,~7)(7,~12,~8)(8,11,~9)(12,~1,~0)"
            sage: cols = "RBBBBBBRBBBRBRR"
            sage: T = VeeringTriangulation(t, cols)
            sage: F = T.flat_structure_middle().layout()
            sage: F.set_pos()
            sage: for _ in range(100):
            ....:     e = randrange(15)
            ....:     _ = F.glue_edge(e)
            ....:     F._check()
        """
        holonomies = self._triangulation._holonomies
        fp = self._triangulation.face_permutation(False)
        ep = self._triangulation.edge_permutation(False)
        a = ep[e]
        b = fp[a]
        c = fp[b]

        if a == e:
            return

        if a != e and holonomies[a] == holonomies[e]:
            # apply point-symmetry to the triangle (a,b,c)
            holonomies[a] *= -1
            holonomies[b] *= -1
            holonomies[c] *= -1

        if self._edge_is_boundary(e):
            pos = self._pos
            if pos is None or pos[e] is None:
                raise RuntimeError

            if a != e:
                pos[a] = pos[e] + holonomies[e]
            pos[b] = pos[a] + holonomies[a]
            pos[c] = pos[b] + holonomies[b]
            xmin = min(pos[a][0], pos[b][0], pos[c][0])
            xmax = max(pos[a][0], pos[b][0], pos[c][0])
            ymin = min(pos[a][1], pos[b][1], pos[c][1])
            ymax = max(pos[a][1], pos[b][1], pos[c][1])

            self._check()
            return xmin, xmax, ymin, ymax

    def sublayout(self, e):
        r"""
        Return the edges of the connected part of the layout reachable from e
        but not going through e

        EXAMPLES::

            sage: from veerer import *

        """
        if self._pos is None:
            return

        T = self._triangulation
        n = T._n
        ep = T._ep
        fp = T._fp
        holonomies = T._holonomies
        pos = self._pos

        triangle_list = []        # list of triangles
        seen = [False] * n        # visited edges (not needed unless we have 2pi angles)
        seen[e] = seen[ep[e]] = True
        todo = [e]
        while todo:
            # run around triangle
            e = todo.pop()
            a = fp[e]
            b = fp[a]
            triangle_list.append((e, a, b))
            # traverse a and b if needed
            if not seen[a]:
                A = ep[a]
                seen[a] = seen[A] = True
                if a != A and pos[A] == pos[a] + holonomies[a]:
                    todo.append(A)
            if not seen[b]:
                B = ep[b]
                seen[b] = seen[B] = True
                if b != B and pos[B] == pos[b] + holonomies[b]:
                    todo.append(B)

        return triangle_list

    def glue_subsurface(self, e):
        r"""
        Move the sublayout glued from `e` to its paired half edge.
        """
        if self._pos is None:
            return

        T = self._triangulation
        ep = T._ep
        holonomies = T._holonomies

        pos = self._pos

        E = ep[e]
        if e == E or pos[E] == pos[e] + holonomies[e]:
            # self-glued or already glued edge
            return

        if holonomies[e] == holonomies[E]:
            # reflection
            s = -1
            t = pos[E] + holonomies[E] + pos[e]
        else:
            # translation
            s = 1
            t = pos[E] + holonomies[E] - pos[e]

        for (a, b, c) in self.sublayout(e):
            if s == -1:
                T.triangle_upside_down(a)
            pos[a] = s * pos[a] + t
            pos[b] = s * pos[b] + t
            pos[c] = s * pos[c] + t

        self._check()

    def set_pos_outgoing_separatrices(self, x_space=1):
        r"""
        Set position according to ...
        """
        if not self._triangulation._translation:
            raise ValueError("only available for translation surface, not for half-translation")

        # determine outgoing separatrices
        T = self._triangulation
        V = T._V
        n = T.num_half_edges()
        pos = self._pos = [None] * n
        vp = T.vertex_permutation(copy=False)
        ep = T.edge_permutation(copy=False)
        fp = T.face_permutation(copy=False)
        hol = T._holonomies

        wedges = []
        bdries = [False] * n
        for e in range(n):
            f = vp[e]
            if T.edge_colour(e) == RED and \
               T.edge_colour(f) == BLUE and \
               hol[e][1] > 0:
                assert hol[f][1] > 0
                wedges.append((e, f))
                bdries[e] = bdries[f] = bdries[ep[e]] = bdries[ep[f]] = True

        x = 0   # current x coordinate
        for er, el in wedges:
            modif = []
            xmin = x
            xmax = x
            v = V((x, 0))
            El = ep[el]
            e = fp[er]
            pos[er] = v
            pos[El] = v + hol[el]
            pos[e] = v + hol[er]
            modif.extend((er, El, e))
            xmax = max(xmax, pos[e][0])
            xmin = min(xmin, pos[El][0])
            todo = []
            if not bdries[e]:
                todo.append(ep[e])
            while todo:
                a = todo.pop()
                b = fp[a]
                c = fp[b]
                A = ep[a]
                pos[a] = pos[A] + hol[A]
                pos[b] = pos[a] + hol[a]
                pos[c] = pos[b] + hol[b]
                modif.extend((a, b, c))
                xmax = max([xmax, pos[a][0], pos[b][0], pos[c][0]])
                xmin = min([xmin, pos[a][0], pos[b][0], pos[c][0]])
                if not bdries[b]:
                    todo.append(ep[b])
                if not bdries[c]:
                    todo.append(ep[c])
            v = V((x - xmin, 0))
            for e in modif:
                pos[e] += v
            x = x - xmin + xmax + x_space

        self._pos = pos

    def set_pos(self, cylinders=None, y_space=0.1):
        r"""
        Set position randomly.

        INPUT:

        - ``cylinders`` - an optional list of cylinders

        EXAMPLES::

            sage: from veerer import *

            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "BBR")
            sage: F = T.flat_structure_min(allow_degenerations=True) # not tested
            sage: F.set_pos(cylinders=T.cylinders(BLUE) + T.cylinders(RED)) # not tested
            sage: F.plot() # not tested
            Graphics object consisting of 13 graphics primitives

        TO be tested

        RBBRBRBRBRBR_791508a432b6_ba2987654310
        RBBBRBBRBBBR_64a5978301b2_ba9875643210
        """
        n = self._triangulation.num_half_edges()
        half_edge_to_face = [None] * n
        faces = self._triangulation.faces()

        for i, (e0, e1, e2) in enumerate(faces):
            half_edge_to_face[e0] = i
            half_edge_to_face[e1] = i
            half_edge_to_face[e2] = i

        nf = self._triangulation.num_faces()
        face_seen = [False] * nf
        fp = self._triangulation.face_permutation(copy=False)
        ep = self._triangulation.edge_permutation(copy=False)
        vectors = self._triangulation._holonomies
        pos = self._pos = [None] * (3 * nf)

        y = 0   # current height

        # set the cylinders
        if cylinders:
            for cyl, _, _, _ in cylinders:
                xmin = xmax = ymin = ymax = self._K.zero()
                a = cyl[0]
                pos[a] = self._triangulation._V((0, 0))

                if a == ep[a]:
                    cyl = cyl[:-1]

                for e in cyl:
                    x0min, x0max, y0min, y0max = self.glue_edge(e)
                    face_seen[half_edge_to_face[ep[e]]] = True
                    xmin = min(xmin, x0min)
                    xmax = max(xmax, x0max)
                    ymin = min(ymin, y0min)
                    ymax = max(ymax, y0max)

                # 2. translate positions according to bounding box
                #    and the current height position
                t = self._triangulation._V((-xmin, y - ymin))
                for e in cyl:
                    a = ep[e]
                    b = fp[a]
                    c = fp[b]
                    pos[a] += t
                    pos[b] += t
                    pos[c] += t

                # 3. set new starting height
                y += ymax - ymin + y_space

        n = self._triangulation.num_half_edges()
        for e in range(n):
            f = half_edge_to_face[e]
            if (pos[e] is None) != (face_seen[half_edge_to_face[e]] is False):
                a, b, c = faces[f]
                raise RuntimeError('cylinder badly set: pos[%s] = %s while its face (%s,%s,%s) is %s' % (e, self._pos[e], a, b, c, 'seen' if face_seen[half_edge_to_face[e]] else 'unseen'))

        # random forest
        while any(x is False for x in face_seen):
            # choose a starting face
            for start in range(nf):
                if face_seen[start] is False:
                    break

            # sets its position
            a, b, c = faces[start]
            pos[a] = self._triangulation._V.zero()
            pos[b] = pos[a] + vectors[a]
            pos[c] = pos[b] + vectors[b]
            xmin = xmax = ymin = ymax = 0
            xmin = min(xmin, pos[a][0], pos[b][0], pos[c][0])
            xmax = max(xmax, pos[a][0], pos[b][0], pos[c][0])
            ymin = min(ymin, pos[a][1], pos[b][1], pos[c][1])
            ymax = max(ymax, pos[a][1], pos[b][1], pos[c][1])

            # spanning tree
            edges = {t: [] for t in range(nf)}
            wait = [start]
            face_seen[start] = True
            q = [start]
            while wait:
                shuffle(wait)
                t = wait.pop()
                t_edges = list(faces[t])
                shuffle(t_edges)
                for e in t_edges:
                    a = ep[e]
                    assert pos[e] is not None
                    assert half_edge_to_face[e] == t
                    f = half_edge_to_face[a]
                    if face_seen[f] is False:
                        edges[t].append(e)
                        face_seen[f] = True
                        q.append(f)
                        wait.append(f)

                        x0min, x0max, y0min, y0max = self.glue_edge(e)
                        xmin = min(xmin, x0min)
                        xmax = max(xmax, x0max)
                        ymin = min(ymin, y0min)
                        ymax = max(ymax, y0max)

            # translate
            t = self._triangulation._V((-xmin, y - ymin))
            for f in q:
                a, b, c = faces[f]
                pos[a] += t
                pos[b] += t
                pos[c] += t

            # set new height
            y += ymax - ymin + y_space

        n = self._triangulation.num_half_edges()
        for e in range(n):
            if self._pos[e] is None:
                raise RuntimeError('pos[%s] not set properly' % e)

    def relabel(self, p):
        r"""
        Relabel according to the permutation ``p``.

        EXAMPLES::

            sage: import veerer
            sage: T = veerer.VeeringTriangulation("(0,1,2)(~0,~1,~2)", [1, 2, 2])
            sage: F = T.flat_structure_middle().layout()
            sage: F.set_pos()
            sage: F._check()
            sage: F.relabel('(0,1,2)')
            sage: F._check()
        """
        n = self._triangulation._n
        if not perm_check(p, n):
            p = perm_init(p, n, self._triangulation.edge_permutation(False))
            if not perm_check(p, n):
                raise ValueError('invalid relabeling permutation')

        self._triangulation.relabel(p)
        if self._pos is not None:
            perm_on_list(p, self._pos)

    def forward_flippable_edges(self):
        return self._triangulation.forward_flippable_edges()

    def backward_flippable_edges(self):
        return self._triangulation.backward_flippable_edges()

    def flip(self, e, folded_edge_convention=RIGHT):
        r"""
        Flip an edge of the layout of a flat triangulation.

        EXAMPLES::

            sage: from veerer import Triangulation, FlatVeeringTriangulation, LEFT, RIGHT
            sage: t = Triangulation("(0,~5,4)(1,2,~6)(3,5,6)")
            sage: holonomies = [(8, 1), (7, 2), (-4, 2), (9, 2), (-2, -3), (-6, 2), (-3, -4), (-3, -4), (-6, 2)]
            sage: l = FlatVeeringTriangulation(t, holonomies).layout()
            sage: l.set_pos()
            sage: l.flip(0, folded_edge_convention=RIGHT)

            sage: l = FlatVeeringTriangulation(t, holonomies).layout()
            sage: l.set_pos()
            sage: l.flip(0, folded_edge_convention=LEFT)
            sage: l._check()

            sage: t = Triangulation("(0,1,2)(~0,3,4)(~1,5,6)")
            sage: holonomies = [(-47,51), (-27, -67), (74, 16), (22, 79),
            ....:               (-69, -28), (-61, -31), (34, -36)]
            sage: l = FlatVeeringTriangulation(t, holonomies).layout()
            sage: l.set_pos()
            sage: l.flip(2)
            sage: l._check()
        """
        if not self._triangulation.is_forward_flippable(e):
            raise ValueError("not flippable")

        # be sure that e is in the middle of a quadrilateral
        # for the layout
        ep = self._triangulation.edge_permutation(copy=False)
        fp = self._triangulation.face_permutation(copy=False)
        pos = self._pos
        holonomies = self._triangulation._holonomies
        E = ep[e]

        if e != E and pos[e] + holonomies[e] != pos[E]:
            raise ValueError("can not flip an unglued edge")

        # -------------
        # \       e  /
        #  \a       /
        #   \      /
        #    \   b/
        #     \  /
        a = fp[e]
        A = ep[a]
        b = fp[a]
        B = ep[b]

        aglued = (a != A) and (pos[A] == pos[a] + holonomies[a])
        bglued = (b != B) and (pos[B] == pos[b] + holonomies[b])

        # Do the flip on the triangulation first
        self._triangulation.flip(e, folded_edge_convention=folded_edge_convention)

        if pos is not None:
            holonomies = self._triangulation._holonomies
            fp = self._triangulation.face_permutation(copy=False)
            if e == E:
                # folded edge
                #                          /  |      |  \
                #                         /b  |      |   \
                #                        /    |      |e  a\
                #   E                   /     |      |     \
                # -------------               |      |
                # \       e  /          \     |  or  |     /
                #  \a       /     ->     \a   |      |    /
                #   \      /              \  e|      |   /
                #    \   b/                \  |      | b/
                #     \  /                  \ |      | /
                #                        LEFT        RIGHT
                assert fp[e] == b and fp[b] == a
                if folded_edge_convention == RIGHT:
                    # we should move the a neighbors
                    A = ep[a]
                    pos[a] = pos[e]
                    pos[e] = pos[a] + holonomies[a]
                    if aglued:
                        self.glue_subsurface(A)
                else:
                    # we should move the b neighbors
                    B = ep[b]
                    pos[e] = pos[b]
                    pos[b] = pos[e] + holonomies[e]
                    if bglued:
                        self.glue_subsurface(B)

            else:
                # normal edge
                #    /  \             / | \
                #   /b   \           /b |  \
                #  /     a\         /   |E a\
                # /  e     \       /    |    \
                # ----------   ->       |
                # \     E  /       \   e|    /
                #  \c     /         \c  |   /
                #   \   d/           \  | d/
                #    \  /             \ | /
                b = fp[e]
                d = fp[E]

                pos[e] = pos[d]
                pos[E] = pos[b]

        self._check()

    ###################################################################
    # Plotting functions
    ###################################################################

    def _plot_edge(self, e, **opts):
        r"""
        Plot the edge ``e``.
        """
        assert self._pos is not None

        pos = self._pos
        fp = self._triangulation._fp
        ep = self._triangulation._ep

        oopts = {}
        E = ep[e]
        u = pos[e]
        v = pos[fp[e]]
        if not self._edge_is_boundary(e):
            # the edge is between adjacent faces
            oopts['alpha'] = 0.5
            oopts['linestyle'] = 'dotted'

        oopts['color'] = EDGE_COLORS[self._edge_slope(e)]

        oopts.update(opts)

        L = line2d([u, v], **oopts)

        if e == E:
            # folded edge
            L += point2d([(u + v) / 2], color='black', marker='x', pointsize=100)

        return L

    def _tikz_edge(self, output, e, edge_label, **opts):
        r"""
        Plot the edge ``e``.
        """
        assert self._pos is not None
        assert edge_label is True or edge_label is False

        pos = self._pos
        fp = self._triangulation._fp
        ep = self._triangulation._ep

        E = ep[e]
        u = pos[e]
        v = pos[fp[e]]

        edge_tikz_opts = TIKZ_EDGE_COLORS[self._edge_slope(e)]
        if not self._edge_is_boundary(e):
            # the edge is between adjacent faces
            edge_tikz_opts += ',dotted'

        if edge_label:
            if self._edge_is_boundary(e) and e > E:
                lab = "$\\sim%d$" % (E)
            else:
                lab = "$%d$" % e

            output.write('\\draw[%s] (%s,%s) -- node[sloped] {\\contour{white}{%s}} (%s,%s);\n' % (edge_tikz_opts, u[0], u[1], lab, v[0], v[1]))

        else:
            output.write('\\draw[%s] (%s,%s) -- (%s,%s);\n' % (edge_tikz_opts, u[0], u[1], v[0], v[1]))

        if e == E:
            # folded edge
            output.write('\\draw mid +(.1,.1) -- (-.1,-.1);\n')
            output.write('\\draw mid +(-.1,.1) -- (.1,-.1);\n')

    def _plot_edge_label(self, a, tilde=None, **opts):
        assert self._pos is not None

        fp = self._triangulation.face_permutation(copy=False)
        ep = self._triangulation.edge_permutation(copy=False)
        pos = self._pos
        vectors = self._triangulation._holonomies

        b = fp[a]
        c = fp[b]

        posa = pos[a].n()
        posb = pos[b].n()
        vc = vectors[c].n()
        vc /= vc.norm()
        relposc = posa - vc

        if a == ep[a]:
            # folded edge
            pos = (6.5 * posa + 6.5 * posb + relposc) / 14
        else:
            pos = (8 * posa + 5 * posb + relposc) / 14

        if tilde is None:
            tilde = self._edge_is_boundary(a)

        if tilde:
            if a > ep[a]:
                lab = "%s=~%s" % (a, ep[a])
            else:
                lab = str(a)
        else:
            lab = str(a)

        x, y = self._triangulation._holonomies[a]
        if y.is_zero():
            angle = 0
        elif x.is_zero():
            if y > 0:
                angle = 90
            else:
                angle = 270
        else:
            x = float(x)
            y = float(y)
            angle = math.acos(x / math.sqrt(x**2 + y**2))
            if y < 0:
                angle *= -1.0
            angle *= 180 / math.pi

        return text(lab, pos, rotation=angle, color='black')

    def _plot_face(self, a, edge_labels=True, fill=True):
        r"""
        Plot the face that contains the edge ``a``.
        """
        assert self._pos is not None

        fp = self._triangulation.face_permutation(copy=False)
        b = fp[a]
        c = fp[b]

        G = Graphics()

        if fill:
            # computing slopes in order to determine filling color
            pos = self._pos

            nred = nblue = 0
            for e in (a, b, c):
                slope = self._edge_slope(e)
                if slope == RED:
                    nred += 1
                elif slope == BLUE:
                    nblue += 1

            if nred == 2:
                color = (.75, 0.5, 0.5)
            elif nblue == 2:
                color = (.5, 0.5, .75)
            else:
                color = (.5, 0.5, 0.5)

            G += polygon2d([pos[a], pos[b], pos[c], pos[a]], alpha=0.41, color=color)

        if edge_labels:
            G += self._plot_edge_label(a)
            G += self._plot_edge_label(b)
            G += self._plot_edge_label(c)

        return G

    def _tikz_face(self, output, a,
            red='red!20', blue='blue!20', neutral='gray!20', tikz_face_options=None):
        assert self._pos is not None

        fp = self._triangulation.face_permutation(copy=False)
        b = fp[a]
        c = fp[b]
        pos = self._pos

        # computing slopes in order to determine filling color
        nred = nblue = 0
        for e in (a, b, c):
            slope = self._edge_slope(e)
            if slope == RED:
                nred += 1
            elif slope == BLUE:
                nblue += 1

        if nred == 2:
            color = TIKZ_FACE_COLORS[RED]
        elif nblue == 2:
            color = TIKZ_FACE_COLORS[BLUE]
        else:
            color = 'gray!20'

        output.write('\\fill[%s] (%s,%s) -- (%s,%s) -- (%s,%s) -- cycle;\n' % (
            color if tikz_face_options is None else color + ',' + tikz_face_options,
            pos[a][0], pos[a][1],
            pos[b][0], pos[b][1],
            pos[c][0], pos[c][1]))

    def _plot_train_track(self, slope):
        pos = self._pos
        vectors = self._triangulation._holonomies

        V2 = VectorSpace(RDF, 2)
        G = Graphics()

        if slope == HORIZONTAL:
            POS = RED
            NEG = BLUE
            color = 'purple'
        else:
            POS = BLUE
            NEG = RED
            color = 'green'

        for e0, e1, e2 in self._triangulation.faces():
            # determine the large edge
            v0 = vectors[e0]
            v1 = vectors[e1]
            v2 = vectors[e2]
            col0 = RED if v0[0] * v0[1] > 0 else BLUE
            col1 = RED if v1[0] * v1[1] > 0 else BLUE
            col2 = RED if v2[0] * v2[1] > 0 else BLUE
            if col0 == POS and col1 == NEG:
                # e2 is large
                l, s1, s2 = e2, e0, e1
            elif col1 == POS and col2 == NEG:
                # i is large
                l, s1, s2 = e0, e1, e2
            elif col2 == POS and col0 == NEG:
                # j is large
                l, s1, s2 = e1, e2, e0

            pl = V2(pos[l])
            ps1 = V2(pos[s1])
            ps2 = V2(pos[s2])

            cl = (pl + ps1) / 2
            vl = (ps1 - pl)
            cs1 = (ps1 + ps2) / 2
            vs1 = ps2 - ps1
            cs2 = (ps2 + pl) / 2
            vs2 = pl - ps2

            ol = V2((-vl[1], vl[0]))
            ol /= ol.norm()
            os1 = V2((-vs1[1], vs1[0]))
            os1 /= os1.norm()
            os2 = V2((-vs2[1], vs2[0]))
            os2 /= os2.norm()

            G += bezier_path([[cl, cl + 0.3 * ol,
                               cs1 + 0.3 * os1, cs1]], rgbcolor=color)
            G += bezier_path([[cl, cl + 0.3 * ol,
                               cs2 + 0.3 * os2, cs2]], rgbcolor=color)
        return G

    def _tikz_train_track(self, slope):
        raise NotImplementedError

    def plot(self, horizontal_train_track=False, vertical_train_track=False, edge_labels=True, fill=True):
        r"""
        Return a graphics.

        INPUT:

        - ``horizontal_train_track`` - boolean - whether to plot the horizontal
          train-track on the surface

        - ``vertical_train_track`` - boolean - whether to plot the vertical
          train-track on the surface

        EXAMPLES::

            sage: from veerer import *
            sage: faces = "(0, ~3, 4)(1, 2, ~7)(3, ~1, ~2)(5, ~8, ~4)(6, ~5, 8)(7, ~6, ~0)"
            sage: colours = 'RBRRBRBRB'
            sage: T = VeeringTriangulation(faces, colours)
            sage: FS = T.flat_structure_min().layout()
            sage: FS.plot(fill=False)
            Graphics object consisting of 31 graphics primitives
            sage: FS.plot(fill=True)
            Graphics object consisting of 37 graphics primitives
            sage: FS.plot(horizontal_train_track=True)
            Graphics object consisting of 49 graphics primitives
            sage: FS.plot(vertical_train_track=True)
            Graphics object consisting of 49 graphics primitives
        """
        if self._pos is None:
            self.set_pos()

        G = Graphics()
        n = self._triangulation.num_half_edges()
        fp = self._triangulation.face_permutation(copy=False)
        ep = self._triangulation.edge_permutation(copy=False)

        # 1. plot faces
        for e in range(n):
            if e < fp[e] and e < fp[fp[e]]:
                G += self._plot_face(e, edge_labels=edge_labels, fill=fill)

        # 2. plot edges
        for e in range(n):
            if self._edge_is_boundary(e) or e <= ep[e]:
                G += self._plot_edge(e)

        # 3. possibly plot train tracks
        if horizontal_train_track:
            G += self._plot_train_track(HORIZONTAL)
        if vertical_train_track:
            G += self._plot_train_track(VERTICAL)

        G.set_aspect_ratio(1)
        G.axes(False)
        return G

    def tikz(self, filename=None, horizontal_train_track=False, vertical_train_track=False, edge_labels=True):
        if filename is None:
            from sys import stdout as output
        else:
            output = open(filename, 'w')

        if self._pos is None:
            self.set_pos()

        n = self._triangulation.num_half_edges()
        fp = self._triangulation.face_permutation(copy=False)
        ep = self._triangulation.edge_permutation(copy=False)

        # 1. plot faces
        for e in range(n):
            if e < fp[e] and e < fp[fp[e]]:
                self._tikz_face(output, e)

        # 2. plot edges
        for e in range(n):
            if self._edge_is_boundary(e) or e <= ep[e]:
                self._tikz_edge(output, e, edge_labels)

        # 3. possibly plot train tracks
        if horizontal_train_track:
            self._tikz_train_track(output, HORIZONTAL)
        if vertical_train_track:
            self._tikz_train_track(output, VERTICAL)

        if filename is not None:
            output.close()

    def train_track(self, slope=VERTICAL):
        r"""
        Return the measured train-track corresponding to the vertical
        (or horizontal) direction.

        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,3)", "BRRR")
            sage: F = T.flat_structure_min().layout()
            sage: F.train_track()
            MeasuredTrainTrack(Triangulation("(0,1,2)(3,~0,~1)"), (1, 1, 2, 2, 1, 1))
            sage: F.train_track(HORIZONTAL)
            MeasuredTrainTrack(Triangulation("(0,1,2)(3,~0,~1)"), (1, 2, 1, 1, 2, 1))
        """
        from .measured_train_track import MeasuredTrainTrack
        if slope == VERTICAL:
            return MeasuredTrainTrack(self._triangulation,
                                      [x.abs() for x, y in self._triangulation._holonomies])
        elif slope == HORIZONTAL:
            return MeasuredTrainTrack(self._triangulation,
                                      [y.abs() for x, y in self._triangulation._holonomies])

    def plot_orbit(self, p, n, slope=VERTICAL, **kwds):
        r"""
        Plot a piece of orbit of the vertical flow.

        EXAMPLES::

            sage: from veerer import *
            sage: T = VeeringTriangulation("(0,1,2)(~0,~1,3)", "BRRR")
            sage: F = T.flat_structure_min().layout()
            sage: F.plot() + F.plot_orbit((0,1/4),4)  # not tested (warning in matplotlib)
            Graphics object consisting of 19 graphics primitives
        """
        G = Graphics()

        if slope == HORIZONTAL:
            raise NotImplementedError

        tt = self.train_track(slope)
        L = tt.lengths()
        V = self._triangulation._holonomies

        O = [q for q in itertools.islice(tt.orbit(p), 2 * n)]

        if self._pos is None:
            self.set_pos()

        for i in range(n):
            i1, x1 = O[2 * i]
            i2, x2 = O[2 * i + 1]

            p1 = self._pos[i1]
            p2 = self._pos[i2]

            G += line2d([p1 + x1 / L[i1] * V[i1], p2 + x2 / L[i2] * V[i2]], **kwds)

        return G
