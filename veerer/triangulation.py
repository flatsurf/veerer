r"""
Triangulation of surfaces.
"""
# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2018 Mark Bell
#                     2018-2024 Vincent Delecroix
#                     2024 Kai Fu
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

import numbers
from array import array

from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE, rich_to_bool

from .permutation import (perm_init, perm_check, perm_cycles, perm_dense_cycles,
                          perm_invert, perm_conjugate, perm_cycle_string, perm_cycles_lengths,
                          perm_cycles_to_string, perm_on_list,
                          perm_num_cycles, str_to_cycles, str_to_cycles_and_data, perm_compose, perm_from_base64_str,
                          uint_base64_str, uint_from_base64_str, perm_base64_str,
                          perms_are_transitive, triangulation_relabelling_from)
from .constellation import Constellation


def face_edge_perms_init(faces):
    r"""
    INPUT:: ``faces`` - a list or a string encoding a permutation

    EXAMPLES::

        sage: from veerer.triangulation import face_edge_perms_init  # random output due to deprecation warnings from realalg

        sage: face_edge_perms_init('(0,1,2)(~0,~1,~2)')
        (array('i', [1, 2, 0, 5, 3, 4]), array('i', [5, 4, 3, 2, 1, 0]))

        sage: face_edge_perms_init('(0,1,2)')
        (array('i', [1, 2, 0]), array('i', [0, 1, 2]))

        sage: face_edge_perms_init('(0,1,2)(~0)(~1)(~2)')
        (array('i', [1, 2, 0, 3, 4, 5]), array('i', [5, 4, 3, 2, 1, 0]))

    TESTS:

    Check that edge permutation do not depend on the details of faces::

        sage: f1 = "(0,~5,4)(3,5,6)(1,2,~6)"
        sage: f2 = "(0,6,5)(1,2,~6)(3,4,~5)"
        sage: f3 = "(6,4,3)(~6,~5,1)(5,0,2)"
        sage: assert face_edge_perms_init(f1)[1] == face_edge_perms_init(f2)[1] == face_edge_perms_init(f3)[1]
    """
    if isinstance(faces, str):
        l = str_to_cycles(faces)
    else:
        l = [[int(i) for i in c] for c in faces]

    pos = []
    neg = []
    for c in l:
        for e in c:
            if e < 0:
                neg.append(e)
            else:
                pos.append(e)

    for e in neg:
        if ~e not in pos:
            raise ValueError("inconsistent permutation data")

    pos.sort()
    neg.sort(reverse=True)
    for i in range(len(pos) - 1):
        if pos[i] == pos[i+1]:
            raise ValueError("repeated edge label {}".format(pos[i]))
        elif pos[i + 1] != pos[i] + 1:
            raise ValueError("missing edge label {} (pos={})".format(pos[i] + 1, pos))
    for i in range(len(neg) - 1):
        if neg[i] == neg[i+1]:
            raise ValueError("repeated edge label ~{}".format(~neg[i]))

    # number of half edges
    n = len(pos) + len(neg)

    # build the edge permutation ep
    ep = [-1] * n  # edge permutation
    m = n - 1   # label available at the back
    for e in neg:
        E = ~e
        e = m
        if ep[E] != -1:
            raise ValueError("inconsistent permutation data")
        m -= 1
        ep[E] = e
        ep[e] = E
    for i in range(n):
        if ep[i] == -1:
            ep[i] = i

    fp = [-1] * n  # face permutation
    for c in l:
        k = len(c)
        for i in range(k):
            e0 = c[i]
            e1 = c[(i + 1) % k]
            if e0 < 0:
                e0 = ep[~e0]
            if e1 < 0:
                e1 = ep[~e1]
            fp[e0] = e1

    return perm_init(fp), perm_init(ep)


def boundary_init(fp, ep, boundary):
    r"""
    Initialize boundary data given face and edge permutation.

    EXAMPLES:

    A test with folded edges::

        sage: from array import array
        sage: from veerer.triangulation import boundary_init
        sage: fp = array('i', [1, 2, 0, 4, 3])
        sage: ep = array('i', [0, 4, 3, 2, 1])
        sage: boundary_init(fp, ep, {0: 1})
        array('i', [1, 0, 0, 0, 0])
        sage: boundary_init(fp, ep, {1: 1})
        array('i', [0, 1, 0, 0, 0])
        sage: boundary_init(fp, ep, {2: 1})
        array('i', [0, 0, 1, 0, 0])
        sage: boundary_init(fp, ep, {-2: 1})
        array('i', [0, 0, 0, 0, 1])
        sage: boundary_init(fp, ep, {-3: 1})
        array('i', [0, 0, 0, 1, 0])
    """
    n = len(ep)

    if boundary is None:
        return array('i', [0] * n)
    elif isinstance(boundary, (array, tuple, list)):
        if len(boundary) != n:
            raise ValueError('invalid input argument')
        return array('i', boundary)
    elif isinstance(boundary, dict):
        edge_to_half_edge = {}
        for j, c in enumerate(perm_cycles(ep, n)):
            if len(c) == 1:
                edge_to_half_edge[j] = c[0]
            elif len(c) == 2:
                edge_to_half_edge[j] = c[0]
                edge_to_half_edge[~j] = c[1]
            else:
                raise ValueError

        output = array('i', [0] * n)
        for e, v in boundary.items():
            if isinstance(e, str):
                if not e:
                    raise ValueError('keys must be valid edges, got {!r}'.format(e))
                elif e[0] == '~':
                    e = ~int(e[1:])
                else:
                    e = int(e)
            elif isinstance(e, numbers.Integral):
                e = int(e)
            if e > n:
                raise ValueError('keys must be valid edges, got {!r}'.format(e))
            output[edge_to_half_edge[e]] = v
        return output
    else:
        raise TypeError('invalid boundary data')



# NOTE: we don't really care that we have a triangulation here. When
# there is no restriction on the cycle decomposition of _fp we got
# a coding for any cell decomposition (with possibly folded edges).
# flipping here still makes sense
#
#  +                 +                 +               _
#   \b             a/                   \            _/ /
#    \      e      /                     \      e __/  /
#     +-----------+       --------->      +    __/
#    /             \                     /  __/        \
#   /c             d\                   /__/            \
#  +                 +                 +/
#
# and we have operations that consist in removing/adding edges.
# For that purpose, it would be convenient to allow partial
# permutations of {0, 1, ..., n-1}.
#
# TODO: implement: is_removable(), remove(), add()
#
# NOTE: Much of the above is already done in the flatsurf package
# with a quite different encoding, see
#
#    https://github.com/flatsurf/sage-flatsurf


# TODO: surface_dynamics compatible datastructure
# {
#  _n: number of standard edges
#  _m: number of self-glued edges
#  _vp, _fp: array of size 2 _n + _m
#  _boundary: bitset
# }
# boundary requirement: boundary is necessarily glued to a non-boundary
# edges are 0, 1, ..., n+m-1 and half-edges are
# 0, ~0, 1, ~1, ..., (n-1), ~(n-1), n, n+1, ..., n+m-1

class Triangulation(Constellation):
    r"""
    A triangulation of an oriented surface

    with attributes:

    * _n  number of half-edges (an int)
    * _fp  face permutation (an array)
    * _ep  edge permutation (an array)
    * _vp  vertex permutation (an array)

    Each of fp, ep, vp is a permutation of oriented half-edges (aka
    "darts"), given as a function.  Our conventions for the
    permutations are set out in the following figure::

              ~b
            ----->
        w-----------*-----------v
         \             <-----  /
          \ \             b   / /
           \ \ c             / /~a
            \ \     F       / /
             \ v           / v
              *         ^ *
             ^ \       / /
              \ \   a / /
             ~c\ \   / /
                \ \   /
                   \ /
                    u

    Here the face permutation sends a to b, b to c, and c to a.  The
    vertex permutation send a to ~c, b to ~a, and c to ~b.  The edge
    permutation interchanges e and ~e for every edge; the edge e is
    folded if and only if e = ~e.  Thus folded edges are fixed by the
    edge permutation.

    EXAMPLES::

        sage: from veerer import *

        sage: T = Triangulation("(~2, 1, ~0)(~1, 0, 2)", mutable=True)
        sage: T.genus()
        1
        sage: T.num_triangles()
        2
        sage: T.num_vertices()
        1
        sage: T.flip(0)
        sage: T
        Triangulation("(0,~1,~2)(1,2,~0)")
        sage: T.flip(0)
        sage: T
        Triangulation("(0,~2,1)(2,~1,~0)")
        sage: T.flip(0)
        sage: T
        Triangulation("(0,1,2)(~2,~0,~1)")
        sage: T.flip(0)
        sage: T
        Triangulation("(0,2,~1)(1,~0,~2)")

        sage: T.set_immutable()
        sage: T.flip(0)
        Traceback (most recent call last):
        ...
        ValueError: immutable triangulation; use a mutable copy instead

    Non-connected surfaces are allowed::

        sage: Triangulation("(0,1,2)(3,4,5)")
        Triangulation("(0,1,2)(3,4,5)")

    Examples with boundaries::

        sage: Triangulation("(0,1,2)(~0)(~1)(~2)", boundary={"~0": 1, "~1": 1, "~2": 1})
        Triangulation("(0,1,2)", boundary="(~2:1)(~1:1)(~0:1)")
        sage: Triangulation("(0,1,2)", boundary="(~0:1)(~1:1,~2:1)")
        Triangulation("(0,1,2)", boundary="(~2:1,~1:1)(~0:1)")
        sage: Triangulation("(0,1,2)(~0,~1,~2)", boundary={"~0": 0})
        Triangulation("(0,1,2)(~2,~0,~1)")

    Example with boundary and folded edges::

        sage: Triangulation("(0,1,2)", boundary="(~1:1,~2:1)")
        Triangulation("(0,1,2)", boundary="(~2:1,~1:1)")

    Examples with invalid boundaries::

        sage: Triangulation("(0,1,2)(~0,~1,~2)", boundary={"~0": 1, "~1": 1, "~2": 0})
        Traceback (most recent call last):
        ...
        ValueError: invalid boundary data

    """
    __slots__ = ['_bdry']

    def __init__(self, triangles, boundary=None, mutable=False, check=True):
        if isinstance(triangles, Triangulation):
            fp = triangles.face_permutation(copy=True)
            ep = triangles.edge_permutation(copy=True)
            bdry = triangles.boundary_vector(copy=True)
        else:
            if boundary is not None and isinstance(boundary, str):
                boundary_cycles, boundary = str_to_cycles_and_data(boundary)
                if isinstance(triangles, str):
                    triangles = str_to_cycles(triangles)
                else:
                    triangles = list(triangles)
                triangles.extend(boundary_cycles)
            fp, ep = face_edge_perms_init(triangles)
            bdry = boundary_init(fp, ep, boundary)

        Constellation.__init__(self, len(fp), None, ep, fp, (bdry,), mutable, check)

    def _set_data_pointers(self):
        self._bdry = self._data[0]

    @staticmethod
    def from_face_edge_perms(fp, ep, vp=None, boundary=None, mutable=False, check=True):
        r"""
        INPUT:

        - ``fp``, ``ep``, ``vp`` -- the face, edge and vertex permutation

        - ``check`` - boolean (default: ``True``) - if set to ``False`` no
          check are performed

        EXAMPLES::

            sage: from veerer import Triangulation
            sage: from array import array

            sage: fp = array('i', [1, 2, 0, 4, 8, 6, 7, 5, 3])
            sage: ep = array('i', [8, 7, 2, 3, 4, 5, 6, 1, 0])
            sage: vp = array('i', [2, 8, 7, 0, 3, 1, 5, 6, 4])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            doctest:warning
            ...
            UserWarning: the method Triangulation.from_face_edge_perms is deprecated; use the classmethod from_permutations instead
            Triangulation("(0,1,2)(3,4,~0)(5,6,~1)")
        """
        import warnings
        warnings.warn('the method Triangulation.from_face_edge_perms is deprecated; use the classmethod from_permutations instead')

        n = len(fp)
        if vp is None:
            vp = array('i', [-1] * n)
            for i in range(n):
                vp[fp[ep[i]]] = i
        if boundary is None:
            bdry = array('i', [0] * n)
        else:
            bdry = array('i', boundary)

        return Triangulation.from_permutations(vp, ep, fp, (bdry,), mutable, check)

    def _check(self, error=RuntimeError):
        r"""
        TESTS::

            sage: from veerer import Triangulation

            sage: Triangulation("(0,1,3)")
            Traceback (most recent call last):
            ...
            ValueError: missing edge label 2 (pos=[0, 1, 3])

            sage: Triangulation("(0,1,~2)")
            Traceback (most recent call last):
            ...
            ValueError: inconsistent permutation data

            sage: Triangulation("(0)")
            Traceback (most recent call last):
            ...
            ValueError: non-trianglular internal face starting at half-edge i=0

            sage: from array import array
            sage: fp = array('i', [1,2,0])
            sage: ep = array('i', [0,1,2])
            sage: vp = array('i', [1,2,0])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            Traceback (most recent call last):
            ...
            ValueError: fev relation not satisfied at half-edge i=0

            sage: fp = array('i', [1,2,0])
            sage: ep = array('i', [1,2,0])
            sage: vp = array('i', [1,2,0])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            Traceback (most recent call last):
            ...
            ValueError: invalid edge permutation at half-edge i=0 (vp=array('i', [1, 2, 0]) ep=array('i', [1, 2, 0]) fp=array('i', [1, 2, 0]))
        """
        Constellation._check(self, error)
        n = self._n

        for face in self.faces():
            i = face[0]
            if self._bdry[i]:
                if not all(self._bdry[i] for i in face):
                    raise error('invalid boundary data')
            else:
                if any(self._bdry[i] for i in face):
                    raise error('invalid boundary data')
                if len(face) != 3:
                    raise error('non-trianglular internal face starting at half-edge i={}'.format(self._half_edge_string(i)))

    def to_flipper(self):
        r"""
        Return the corresponding flipper triangulation

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.to_flipper()  # optional - flipper
            [(~2, ~0, ~1), (0, 1, 2)]

        Conversion of veering triangulation forgot the colors::

            sage: V = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: V.to_flipper()  # optional - flipper
            [(~2, ~0, ~1), (0, 1, 2)]
        """
        from .features import flipper_feature
        flipper_feature.require()

        if any(self._data[0]):
            raise ValueError('triangulation has boundary')

        import flipper
        ep = self._ep
        F = []
        for f in self.faces():
            face = []
            for e in f:
                if ep[e] == e:
                    raise ValueError("flipper do not accept folded edges")
                if ep[e] < e:
                    face.append(~int(ep[e]))
                else:
                    face.append(int(e))
            F.append(tuple(face))

        return flipper.create_triangulation(F)

    def to_curver(self):
        r"""
        Return the corresponding curver triangulation

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.to_curver()  # optional - curver
            [(~2, ~0, ~1), (0, 1, 2)]

        Conversion of veering triangulation forgot the colors::

            sage: V = VeeringTriangulation("(0,1,2)(~0,~1,~2)", [RED, RED, BLUE])
            sage: V.to_curver()  # optional - curver
            [(~2, ~0, ~1), (0, 1, 2)]
        """
        from .features import curver_feature
        curver_feature.require()

        import curver
        ep = self._ep
        F = []
        for f in self.faces():
            face = []
            for e in f:
                if ep[e] == e:
                    raise ValueError("curver do not accept folded edges")
                if ep[e] < e:
                    face.append(~int(ep[e]))
                else:
                    face.append(int(e))
            F.append(tuple(face))

        return curver.create_triangulation(F)

    def num_triangles(self):
        r"""
        Return the number of triangles.
        """
        return sum(self._data[0][c[0]] == 0 for c in perm_cycles(self._fp))

    def triangles(self):
        r"""
        Return the list of internal faces as triples of half-edges

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.faces()
            [[0, 1, 2], [3, 4, 5], [6, 8, 7]]
        """
        return [c for c in perm_cycles(self._fp, True, self._n) if self._data[0][c[0]] == 0]

    def boundary_faces(self):
        r"""
        Return the list of boundaries as lists of half-edges.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.boundary_faces()
            []
            sage: T = Triangulation("(0,1,2)(~0)(~1,~2)", {"~0": 1, "~1": 1, "~2": 1})
            sage: T.boundary_faces()
            [[3, 4], [5]]
        """
        return [c for c in perm_cycles(self._fp, True, self._n) if self._data[0][c[0]]]

    def num_boundary_faces(self):
        r"""
        Return the number of boundary faces.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.num_boundary_faces()
            0
            sage: T = Triangulation("(0,1,2)(~0)(~1,~2)", {"~0": 1, "~1": 1, "~2": 1})
            sage: T.num_boundary_faces()
            2

            sage: fp = "(0,2,1)(~0,3,~1)"
            sage: bdry = "(~2:2,~3:2)"
            sage: Triangulation(fp, bdry).num_boundary_faces()
            1
        """
        return sum(bool(self._data[0][c[0]]) for c in perm_cycles(self._fp))

    def euler_characteristic(self):
        r"""
        Return the Euler characteristic of this triangulation.

        EXAMPLES::

            sage: from veerer import Triangulation

        A sphere::

            sage: T = Triangulation("(0,1,2)")
            sage: T.euler_characteristic()
            2

        A torus::

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.euler_characteristic()
            0

        A genus 2 surface::

            sage: T = Triangulation("(0,1,2)(~2,3,4)(~4,5,6)(~6,~0,7)(~7,~1,8)(~8,~3,~5)")
            sage: T.euler_characteristic()
            -2

        A cylinder::

            sage: T = Triangulation("(0,1,2)(~0,3,4)(~1,~2)(~3,~4)", {"~1": 1, "~2": 1, "~3": 1, "~4": 1})
            sage: T.euler_characteristic()
            0

        A pair of pants::

            sage: T = Triangulation("(0,1,2)(~0)(~1)(~2)", {"~0": 1, "~1": 1, "~2": 1})
            sage: T.euler_characteristic()
            -1
        """
        return self.num_triangles() - self.num_edges() + (self.num_vertices() + self.num_folded_edges())

    def genus(self):
        r"""
        Return the genus of this triangulation.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,~0)(2,~1,~2)")
            sage: T.genus()
            0

            sage: T = Triangulation([(-3,1,-1), (-2,0,2)])
            sage: T.genus()
            1

        A cylinder::

            sage: T = Triangulation("(0,1,2)(~0,3,4)(~1,~2)(~3,~4)", {"~1": 1, "~2": 1, "~3": 1, "~4": 1})
            sage: T.genus()
            0

        A pair of pants::

            sage: T = Triangulation("(0,1,2)(~0)(~1)(~2)", {"~0": 1, "~1": 1, "~2": 1})
            sage: T.genus()
            0
        """
        if not self.is_connected():
            raise NotImplementedError

        # chi = 2 - 2g - n
        return (2 - self.euler_characteristic() - self.num_boundary_faces()) // 2

    def __str__(self):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: str(T)
            'Triangulation("(0,1,2)(~2,~0,~1)")'
        """
        cycles = perm_cycles(self._fp, n=self._n)
        face_cycles = perm_cycles_to_string([c for c in cycles if not self._data[0][c[0]]], involution=self._ep)
        bdry_cycles = perm_cycles_to_string([c for c in cycles if self._data[0][c[0]]], involution=self._ep, data=self._data[0])
        if bdry_cycles:
            return 'Triangulation("%s", boundary="%s")' % (face_cycles, bdry_cycles)
        else:
            return 'Triangulation("%s")' % face_cycles

    def __repr__(self):
        return str(self)

    def _check_homology_matrix(self, m):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: m = matrix([[1,1],[-1,0],[0,-1]])
            sage: T._check_homology_matrix(m)
        """
        ne = self.num_edges()
        ep = self._ep

        assert m.nrows() == ne
        try:
            V = m.row_ambient_module()
        except AttributeError:
            V = m._row_ambient_module()

        # boundary condition
        for F in self.faces():
            v = V.zero()
            for e in F:
                E = ep[e]
                if E > e:
                    v += m[e]
                elif E < e:
                    v -= m[E]
                else:
                    # folded edge condition
                    # NOTE: it is stupid to keep all these zeros in
                    # the matrix! We should label the folded edges
                    # after the non-folded ones...
                    # though we can allow non-zero stuff assuming
                    # that the half edge is oriented toward the pole
                    assert m[e].is_zero()
            assert v.is_zero()

    # TODO: also compute the intersection pairing!
    def homology_matrix(self):
        r"""
        Return a basis of homology as a matrix.

        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)")
            sage: T.homology_matrix().transpose().echelon_form()
            [ 1  0  1  0  1  1]
            [ 0  1  0  0 -1  0]
            [ 0  0  0  1  0 -1]

            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,6,0)")
            sage: T.homology_matrix().transpose().echelon_form()
            [ 1  0  1  0  1  1  0]
            [ 0  1  0  0 -1  0  0]
        """
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ

        ep = self._ep
        nf = self.num_faces()
        ne = self.num_edges()
        nfe = self.num_folded_edges()
        m = matrix(ZZ, nf + nfe, ne)

        # face equations
        for i, f in enumerate(self.faces()):
            for e in f:
                if ep[e] == e:
                    continue
                elif ep[e] < e:
                    m[i, ep[e]] -= 1
                else:
                    m[i, e] += 1

        # force the folded edge to have coefficient zero
        for e in range(ne):
            if ep[e] == e:
                m[i, e] = 1
                i += 1

        # compute and check
        h = m.right_kernel_matrix().transpose()
        self._check_homology_matrix(h)
        return h

    def flip_homological_action(self, e, m, twist=False):
        r"""
        Multiply the matrix ``m`` on the left by the homology action of
        the ``e``-flip.

        The matrix ``m`` must have ``ne`` rows and each column represents a
        vector in cohomology (possibly twisted for quadratic differentials).
        That is to say, each column has to satisfy the following conditions:

        - for each face `F` the evaluation on its boundary is zero (up to twist).

        - for each folded edge, the corresponding entry is zero (up to twist).

        INPUT:

        - ``e`` - a half edge

        - ``m`` - matrix

        - ``twist`` - whether we consider the twisted homology (mainly useful for
          veering triangulations)

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)", mutable=True)
            sage: A = matrix([[1,1],[-1,0],[0,-1]])
            sage: B = copy(A)
            sage: for e in [0,1,0,1]:
            ....:     T.flip_homological_action(e, B)
            ....:     T.flip(e)
            sage: B
            [ 1 -3]
            [-1  4]
            [ 0 -1]

        In the above example we are back to the initial triangulation and
        one can recover the homology matrix using pseudo inverses::

            sage: A.pseudoinverse() * B
            [ 1 -4]
            [ 0  1]

        Another torus example (a single Dehn twist along the curve
        w)::

            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)", mutable=True)
            sage: u = [0,1,0,0,-1,0]
            sage: v = [0,0,0,1,0,-1]
            sage: w = [1,1,1,1,0,0]
            sage: A = matrix(ZZ, 3, 6, [u,v,w]).transpose()
            sage: B = copy(A)
            sage: for i in (0,2,1,3):
            ....:     T.flip_homological_action(i, B)
            ....:     T.flip(i)
            ....:     T._check_homology_matrix(B)
            sage: p = "(0,2)(~0,~2)(1,3)(~1,~3)"
            sage: T.relabel_homological_action(p, B)
            sage: T.relabel(p)
            sage: T._check_homology_matrix(B)
            sage: A.pseudoinverse() * B
            [1 0 0]
            [0 1 0]
            [1 1 1]
        """
        ne = self.num_edges()
        assert m.nrows() == ne
        ep = self._ep

        if not twist and ep[e] == e:
            return
        elif ep[e] < e:
            e = ep[e]

        a, b, c, d = self.square_about_edge(e)
        # v_e use to be v_c + v_d and becomes v_d + v_a
        # v<----------u     v<----------u
        # |     a    ^^     |^    a     ^
        # |         / |     | \         |
        # |  F     /  |     |  \     G  |
        # |       /   |     |   \       |
        # |b    e/   d| --> |b   \e    d|
        # |     /     |     |     \     |
        # |    /      |     |      \    |
        # |   /       |     |       \   |
        # |  /     G  |     | F      \  |
        # | /         |     |         \ |
        # v/    c     |     v     c    \|
        # w---------->x     w---------->x

        # NOTE: We might want to use optimized row operations such as
        #   swap_rows(i, j)
        #   add_multiple_of_row(i, j, s)

        A = ep[a]
        D = ep[d]
        if twist:
            m[e] = (m[a] if a < A else m[A]) + (m[d] if d < D else m[D])
        else:
            if a == A and d == D:
                try:
                    m[e] = m.row_ambient_module().zero()
                except AttributeError:
                    m[e] = m._row_ambient_module().zero()
            elif a == A:
                m[e] = m[d] if d < D else -m[D]
            elif d == D:
                m[e] = m[a] if a < A else -m[A]
            else:
                m[e] = (m[a] if a < A else -m[A]) + (m[d] if d < D else -m[D])

    def relabel_homological_action(self, p, m, twist=False, check=True):
        r"""
        Acts on the homological vectors ``m`` by the relabeling permutation ``p``.

        The matrix ``m`` must have ``ne`` rows and each column represents a
        vector in cohomology. That is to say, each column has to satisfy the
        following conditions:

        - for each face `F` the evaluation on its boundary is zero.

        - for each folded edge, the corresponding entry is zero.

        INPUT:

        - ``p`` - automorphism

        - ``m`` - matrix

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)", mutable=True)
            sage: u = [0,1,0,0,-1,0]
            sage: v = [0,0,0,1,0,-1]
            sage: w = [1,1,1,1,0,0]
            sage: A = matrix(ZZ, 3, 6, [u,v,w]).transpose()
            sage: T._check_homology_matrix(A)
            sage: perms = ["(0,1)(~0,~1)", "(3,~3)", "(3,~5)(~3,5)",
            ....:          "(0,1,2,~0,~1,~2)",
            ....:          "(0,~1,4,~0,1,~4)(2,3)(~2,~3)"]
            sage: for p in perms * 5:
            ....:     T.relabel_homological_action(p, A)
            ....:     T.relabel(p)
            ....:     T._check_homology_matrix(A)
        """
        n = self._n
        ne = self.num_edges()
        if check and not perm_check(p, n):
            p = perm_init(p, self._n, self._ep)
            if not perm_check(p, n):
                raise ValueError('invalid relabeling permutation')

        q = perm_invert(p)

        ep = self._ep
        seen = [False] * ne

        for e0 in range(ne):
            if seen[e0]:
                continue

            seen[e0] = True

            e = q[e0]
            E = ep[e]
            if E < e:
                is_e0_neg = True
                e, E = E, e
            else:
                is_e0_neg = False
            while not seen[e]:
                assert e < ne
                seen[e] = True

                ee = q[e]
                EE = ep[ee]

                if EE < ee:
                    is_neg = True
                    ee, EE = EE, ee
                else:
                    is_neg = False
                m.swap_rows(e, ee)
                if is_neg and not twist:
                    m[e] *= -1

                e = ee

            # one more sign change?
            assert e == e0
            if is_e0_neg and not twist:
                m[e] *= -1

    def is_flippable(self, e, check=True):
        r"""
        Check whether the half-edge e is flippable.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~2,4)(~1,3,~3)")
            sage: T.is_flippable(0)
            True
            sage: T.is_flippable(1)
            True
            sage: T.is_flippable(3)
            False
            sage: T.is_flippable(4)
            True

        A torus with boundary::

            sage: t = Triangulation("(0,2,1)(3,~1,~0)", boundary="(~3:1,~2:1)")
            sage: t.is_flippable(0)
            True
            sage: t.is_flippable(1)
            True
            sage: t.is_flippable(2)
            False
            sage: t.is_flippable(3)
            False
        """
        if check:
            e = self._check_half_edge(e)
        E = self._ep[e]
        a = self._fp[e]
        b = self._fp[a]
        return not self._data[0][e] and not self._data[0][E] and a != E and b != E and self._data[0][e] == 0 and self._data[0][E] == 0

    def flippable_edges(self):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.flippable_edges()
            [0, 1, 2]
            sage: V = VeeringTriangulation(T, [RED, RED, BLUE])
            sage: V.flippable_edges()
            [0, 1]

            sage: T = Triangulation("(0,1,2)(~0,3,4)(~1,~2)(~3,~4)", {"~1": 1, "~2": 1, "~3": 1, "~4": 1})
            sage: T.flippable_edges()
            [0]
        """
        n = self._n
        ep = self._ep
        return [e for e in range(n) if e <= ep[e] and self.is_flippable(e)]

    def square_about_edge(self, e, check=True):
        r"""
        Return the four edges that makes ``e`` the diagonal of a quadrilateral.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.square_about_edge(0)
            (1, 2, 4, 3)

            sage: T = Triangulation("(0,1,2)")
            sage: T.square_about_edge(0)
            (1, 2, 1, 2)
        """
        # x<----------x
        # |     a    ^^
        # |         / |
        # |        /  |
        # |       /   |
        # |b    e/   d|
        # |     /     |
        # |    /      |
        # |   /       |
        # |  /        |
        # | /         |
        # v/    c     |
        # x---------->x

        if check:
            e = self._check_half_edge(e)

        E = self._ep[e]
        if check and (self._data[0][e] or self._data[0][E]):
            raise ValueError('non internal edge')

        a = self._fp[e]
        b = self._fp[a]
        c = self._fp[E]
        d = self._fp[c]

        return a, b, c, d

    def flip(self, e, check=True):
        r"""
        Flip the edge ``e``.

        EXAMPLES::

            sage: from veerer import Triangulation

        Flipping a folded edge is an involution::

            sage: T = Triangulation("(0,1,2)", mutable=True)
            sage: T
            Triangulation("(0,1,2)")
            sage: T.flip(0)
            sage: T
            Triangulation("(0,2,1)")
            sage: T.flip(0)
            sage: T
            Triangulation("(0,1,2)")

            sage: T == Triangulation("(0,1,2)")
            True

        An example with boundaries::

            sage: t = Triangulation("(0,1,2)(~0,3,4)", boundary="(~4:1,~3:1,~2:1,~1:1)", mutable=True)
            sage: t.flip(0)
            sage: t
            Triangulation("(0,2,3)(1,~0,4)", boundary="(~4:1,~3:1,~2:1,~1:1)")
            sage: t.flip(2)
            Traceback (most recent call last):
            ...
            ValueError: can not flip non internal edge 2
        """
        # v<----------u     v<----------u
        # |     a    ^^     |^    a     ^
        # |         / |     | \         |
        # |  F     /  |     |  \     G  |
        # |       /   |     |   \       |
        # |b    e/   d| --> |b   \     d|
        # |     /     |     |     \     |
        # |    /      |     |     e\    |
        # |   /       |     |       \   |
        # |  /     G  |     | F      \  |
        # | /         |     |         \ |
        # v/    c     |     v     c    \|
        # w---------->x     w---------->x

        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        if check:
            e = self._check_half_edge(e)

        E = self._ep[e]
        if self._data[0][e] or self._data[0][E]:
            raise ValueError('can not flip non internal edge %s' % self._norm(e))

        a = self._fp[e]
        b = self._fp[a]
        if a == E or b == E:
            raise ValueError('edge %s is not flippable' % self._norm(e))
        c = self._fp[E]
        d = self._fp[c]

        A = self._ep[a]
        B = self._ep[b]
        C = self._ep[c]
        D = self._ep[d]

        # Disabled for now
        # F = self._fl[e]
        # G = self._fl[E]

        # v = self._vl[b]
        # x = self._vl[d]

        # fix face perm and cycles
        self._fp[e] = b
        self._fp[b] = c
        self._fp[c] = e
        self._fp[a] = E
        self._fp[E] = d
        self._fp[d] = a

        # Face labels
        # self._fl[a] = G
        # self._fl[c] = F

        # fix vertex perm
        self._vp[a] = D
        self._vp[b] = E
        self._vp[E] = A
        self._vp[c] = B
        self._vp[d] = e
        self._vp[e] = C

        # Vertex labels
        # self._vl[e] = x
        # self._vl[E] = v

    def flip_back(self, e, check=True):
        r"""
        Flip back the edge ``e``.

        EXAMPLES::

            sage: from veerer import *
            sage: T0 = Triangulation([(0,1,2),(-1,-2,-3)], mutable=True)
            sage: T = T0.copy()
            sage: T.flip(0)
            sage: T.flip_back(0)
            sage: T == T0
            True

            sage: T.flip(1)
            sage: T.flip(2)
            sage: T.flip_back(2)
            sage: T.flip_back(1)
            sage: T == T0
            True

            sage: T.set_immutable()
            sage: T.flip(0)
            Traceback (most recent call last):
            ...
            ValueError: immutable triangulation; use a mutable copy instead
        """
        # Use the following for reference:
        # v<----------u     v<----------u
        # |     a    ^^     |\    a     ^
        # |         / |     | \         |
        # |  F     /  |     |  \     F  |
        # |       /   |     |   \       |
        # |b    e/   d| --> |b   \e    d|
        # |     /     |     |     \     |
        # |    /      |     |      \    |
        # |   /       |     |       \   |
        # |  /     G  |     | G      \  |
        # | /         |     |         \ |
        # v/    c     |     v     c    v|
        # w---------->x     w---------->x

        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        if check:
            e = self._check_half_edge(e)

        E = self._ep[e]

        if self._data[0][e] or self._data[0][E]:
            raise ValueError('can not flip non internal edge %s' % self._norm(e))

        a = self._fp[e]
        b = self._fp[a]
        if a == E or b == E:
            raise ValueError('edge %s is not flippable' % self._norm(e))
        c = self._fp[E]
        d = self._fp[c]

        A = self._ep[a]
        B = self._ep[b]
        C = self._ep[c]
        D = self._ep[d]

        # Disabled for now
        # F = self._fl[e]
        # G = self._fl[E]

        # v = self._vl[b]
        # x = self._vl[d]

        # fix face perm and cycles
        self._fp[e] = d
        self._fp[d] = a
        self._fp[a] = e
        self._fp[b] = c
        self._fp[c] = E
        self._fp[E] = b

        # Face labels
        # Disabled for now
        # self._fl[a] = G
        # self._fl[c] = F

        # fix vertex perm
        self._vp[a] = D
        self._vp[b] = e
        self._vp[e] = A
        self._vp[c] = B
        self._vp[d] = E
        self._vp[E] = C

        # Vertex labels
        # Disabled for now
        # self._vl[e] = x
        # self._vl[E] = v

    def conjugate(self):
        r"""
        Conjugate this triangulation.

        The face permutation is replaced by its inverse, conjugated by
        the edge permutation. The vertex permutation is replaced by
        its inverse.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)", mutable=True)
            sage: T.conjugate()
            sage: T
            Triangulation("(0,2,4)(1,3,5)(~5,~4,~3)(~2,~1,~0)")

            sage: T = Triangulation("(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)", mutable=False)
            sage: T.conjugate()
            Traceback (most recent call last):
            ...
            ValueError: immutable triangulation; use a mutable copy instead
        """
        # for reference
        #
        # x<----------x     x
        # |     a    ^^     ^\
        # |         /       | \
        # |        /        |  \
        # |       /         |   \
        # |b    c/      --> |    \
        # |     /           |     \
        # |    /            |B    C\
        # |   /             |       \
        # |  /              |        \
        # | /               |         \
        # v/                |     A    v
        # x                 x<----------x
        #
        # (a, b, c)     -->  (C, B, A)

        if not self._mutable:
            raise ValueError('immutable triangulation; use a mutable copy instead')

        self._fp = perm_conjugate(perm_invert(self._fp), self._ep)
        self._vp = perm_invert(self._vp)

    # TODO: deprecate
    is_isomorphic_to = Constellation.is_isomorphic

    def cover(self, c, mutable=False, check=True):
        from .cover import TriangulationCover
        return TriangulationCover(self, c, mutable=mutable, check=check)

    def _check_xy(self, x, y):
        if len(x) == self.num_edges():
            x = [x[self._norm(i)] for i in range(self._n)]
        elif len(x) != self._n or any(a < 0 for a in x):
            raise ValueError('invalid argument x')
        if len(y) == self.num_edges():
            y = [y[self._norm(i)] for i in range(self._n)]
        elif len(y) != self._n or any(a < 0 for a in y):
            raise ValueError('invalid argument y')
        return (x, y)

    def colouring_from_xy(self, x, y, check=True):
        r"""
        Return the veering colouring associated with the holonomy data ``x`` and ``y``.

        EXAMPLES::

            sage: from veerer import Triangulation
            sage: t = "(0,1,2)(~0,~1,~2)"
            sage: t = Triangulation(t)
            sage: x = [1, 2, 1]
            sage: y = [1, 1, 2]
            sage: t.colouring_from_xy(x, y)
            array('i', [1, 2, 2, 2, 2, 1])
            sage: t = "(0,1,2)(~0,~1,4)(~2,5,3)(~3,~4,~5)"
            sage: t = Triangulation(t)
            sage: x = [1, 2, 1, 2, 1, 1]
            sage: y = [1, 1, 2, 1, 2, 3]
            sage: t.colouring_from_xy(x, y)
            array('i', [1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1])
        """
        from .constants import BLUE, RED, PURPLE, GREEN

        if check:
            x, y = self._check_xy(x, y)

        def check_set_colouring(colouring, e, ep, col):
            if colouring[e] is None:
                colouring[e] = colouring[ep[e]] = col
            elif colouring[e] != col:
                raise ValueError('inconsistent colouring between x and y for edge e={}'.format(e))

        colouring = [None] * self._n
        faces = perm_cycles(self._fp, True, self._n)
        ep = self._ep
        for face in faces:
            a, b, c = face
            x_degenerate = (x[a] == 0) + (x[b] == 0) + (x[c] == 0)
            if x_degenerate == 0:
                if x[a] == x[b] + x[c]:
                    xlarge = 0
                elif x[b] == x[c] + x[a]:
                    xlarge = 1
                elif x[c] == x[a] + x[b]:
                    xlarge = 2
                else:
                    raise ValueError('inconsistent x data for triangle {}'.format(face))
                check_set_colouring(colouring, face[(xlarge + 1) % 3], ep, BLUE)
                check_set_colouring(colouring, face[(xlarge + 2) % 3], ep, RED)
            elif x_degenerate == 1:
                if x[a] == 0:
                    xvert = 0
                elif x[b] == 0:
                    xvert = 1
                elif x[c] == 0:
                    xvert = 2
                check_set_colouring(colouring, face[xvert], ep, GREEN)
                check_set_colouring(colouring, face[(xvert + 1) % 3], ep, RED)
                check_set_colouring(colouring, face[(xvert + 2) % 3], ep, BLUE)
            else:
                raise ValueError('inconsistent x data for triangle {}'.format(face))

            y_degenerate = (y[a] == 0) + (y[b] == 0) + (y[c] == 0)
            if y_degenerate == 0:
                if y[a] == y[b] + y[c]:
                    ylarge = 0
                elif y[b] == y[c] + y[a]:
                    ylarge = 1
                elif y[c] == y[a] + y[b]:
                    ylarge = 2
                else:
                    raise ValueError('inconsistent y data for triangle {}'.format(face))
                check_set_colouring(colouring, face[(ylarge + 1) % 3], ep, RED)
                check_set_colouring(colouring, face[(ylarge + 2) % 3], ep, BLUE)
            elif y_degenerate == 1:
                if y[a] == 0:
                    yhor = 0
                elif y[b] == 0:
                    yhor = 1
                elif y[c] == 0:
                    yhor = 2
                check_set_colouring(colouring, face[yhor], ep, PURPLE)
                check_set_colouring(colouring, face[(yhor + 1) % 3], ep, BLUE)
                check_set_colouring(colouring, face[(yhor + 2) % 3], ep, RED)
            else:
                raise ValueError('inconsistent y data for triangle {}'.format(face))

        return array('i', colouring)
