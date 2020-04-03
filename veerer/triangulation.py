r"""
Triangulation of surfaces.
"""
from __future__ import absolute_import, print_function
from six.moves import range, map, filter, zip

from array import array
from .permutation import *
from .env import require_package, flipper, curver

def face_edge_perms_init(data):
    r"""
    EXAMPLES::

        sage: from veerer.triangulation import face_edge_perms_init

        sage: face_edge_perms_init('(0,1,2)(~0,~1,~2)')
        (array('l', [1, 2, 0, 5, 3, 4]), array('l', [5, 4, 3, 2, 1, 0]))

        sage: face_edge_perms_init('(0,1,2)')
        (array('l', [1, 2, 0]), array('l', [0, 1, 2]))

    TESTS:

    Check that edge permutation do not depend on the details of faces::

        sage: f1 = "(0,~5,4)(3,5,6)(1,2,~6)"
        sage: f2 = "(0,6,5)(1,2,~6)(3,4,~5)"
        sage: f3 = "(6,4,3)(~6,~5,1)(5,0,2)"
        sage: assert face_edge_perms_init(f1)[1] == face_edge_perms_init(f2)[1] == face_edge_perms_init(f3)[1]
    """
    if isinstance(data, str):
        l = str_to_cycles(data)
    else:
        l = [[int(i) for i in c] for c in data]

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
    if pos != list(range(len(pos))):
        raise ValueError("inconsistent permutation data")

    # number of half edges
    n = len(pos) + len(neg)

    # build the edge permutation ep
    ep = [-1] * n # edge permutation
    m = n - 1 # label available at the back
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

    fp = [-1] * n # face permutation
    for c in l:
        k = len(c)
        for i in range(k):
            e0 = c[i]
            e1 = c[(i+1) % k]
            if e0 < 0:
                e0 = ep[~e0]
            if e1 < 0:
                e1 = ep[~e1]
            fp[e0] = e1

    return perm_init(fp), perm_init(ep)

# NOTE: we don't really care that we have a triangulation here. When
# there is no restriction on the cycle decomposition of _fp we got
# a coding for any cell decomposition (with possibly folded edges).
# flipping here still makes sense
#
#  +                 +                 +               _ +
#   \b             a/                   \            _/ /
#    \      e      /                     \      e __/  /
#     +-----------+       --------->      +    __/    +
#    /             \                     /  __/        \
#   /c             d\                   /__/            \
#  +                 +                 +/                +
#
# and we have operations that consist in removing/adding edges.
# For that purpose, it would be convenient to allow partial
# permutations of {0, 1, ..., n-1}.
#
# TODO: implement: is_removable(), remove(), add()
#
# NOTE: Many of the above is already done in the flatsurf package
# with a quite different encoding
#
#    https://github.com/videlec/sage-flatsurf

class Triangulation(object):
    r"""
    A triangulation of an oriented surface

    with attributes:

    * _n  number of half-edges (an int)
    * _fp  face permutation (an array)
    * _ep  edge permutation (an array)
    * _vp  vertex permutation (an array)

    Each of fp, ep, vp is a permutation of oriented half-edges (aka
    "darts"), given as a function.  Our conventions for the
    permutations are set out in the following figure:

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

        sage: T = Triangulation("(~2, 1, ~0)(~1, 0, 2)")
        sage: T.genus()
        1
        sage: T.num_faces()
        2
        sage: T.num_vertices()
        1
        sage: T.flip(0); T
        Triangulation("(0,~1,~2)(1,2,~0)")
        sage: T.flip(0); T
        Triangulation("(0,~2,1)(2,~1,~0)")
        sage: T.flip(0); T
        Triangulation("(0,1,2)(~2,~0,~1)")
        sage: T.flip(0); T
        Triangulation("(0,2,~1)(1,~0,~2)")

    The surface must be connected::

        sage: Triangulation("(0,1,2)(3,4,5)")
        Traceback (most recent call last):
        ...
        ValueError: (fp, ep, vp) do not generate a transitive group
    """
    __slots__ = ['_n', '_fp', '_ep', '_vp']

    def __init__(self, triangles, check=True):
        if isinstance(triangles, Triangulation):
            self._fp = triangles.face_permutation(copy=True)
            self._ep = triangles.edge_permutation(copy=True)
        else:
            self._fp, self._ep = face_edge_perms_init(triangles)

        fp = self._fp
        ep = self._ep
        n = self._n = len(fp)

        # TODO: face labels are disabled for now. We should actually
        # make it a *partition* of the faces (two faces are allowed to
        # have the same label) The trivial partition {0,1,2,...,F-1}
        # would correspond to no labels and the atomic
        # {{0},{1},...,{F-1}} would correspond to everybody being
        # labelled
        # fl = self._fl = [None] * F # face labels

        # TODO: edge labels are disabled for now.
        # el = self._fl = [None] * E # edge labels

        vp = self._vp = array('l', [-1] * n)
        for i in range(n):
            vp[fp[ep[i]]] = i

        # TODO: vertex labels are disabled for now
        # vl = self._vl = perm_cycles(vp)[1]....

        if check:
            self._check(ValueError)

    @staticmethod
    def from_face_edge_perms(fp, ep, vp=None, check=True):
        r"""
        INPUT:

        - ``fp``, ``ep``, ``vp`` -- the face, edge and vertex permutation

        - ``check`` - boolean (default: ``True``) - if set to ``False`` no
          check are performed

        EXAMPLES::

            sage: from veerer import Triangulation
            sage: from array import array

            sage: fp = array('l', [1, 2, 0, 4, 8, 6, 7, 5, 3])
            sage: ep = array('l', [8, 7, 2, 3, 4, 5, 6, 1, 0])
            sage: vp = array('l', [2, 8, 7, 0, 3, 1, 5, 6, 4])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            Triangulation("(0,1,2)(3,4,~0)(5,6,~1)")
        """
        T = Triangulation.__new__(Triangulation)
        n = T._n = len(fp)
        T._fp = fp
        T._ep = ep
        if vp is None:
            fp = T._fp
            ep = T._ep
            vp = array('l', [-1] * n)
            for i in range(n):
                vp[fp[ep[i]]] = i
        T._vp = vp

        if check:
            T._check(ValueError)

        return T

    @staticmethod
    def from_string(s, check=True):
        r"""
        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: Triangulation.from_string(T.to_string()) == T
            True
        """
        n,fp,ep = s.split('_')
        n = uint_from_base64_str(n)
        fp = perm_from_base64_str(fp, n)
        ep = perm_from_base64_str(ep, n)
        return Triangulation.from_face_edge_perms(fp, ep, check=check)

    def _check(self, error=RuntimeError):
        r"""
        TESTS::

            sage: from veerer import Triangulation

            sage: Triangulation("(0,1,3)")
            Traceback (most recent call last):
            ...
            ValueError: inconsistent permutation data

            sage: Triangulation("(0,1,~2)")
            Traceback (most recent call last):
            ...
            ValueError: inconsistent permutation data

            sage: Triangulation("(0)")
            Traceback (most recent call last):
            ...
            ValueError: broken face permutation

            sage: from array import array
            sage: fp = array('l', [1,2,0])
            sage: ep = array('l', [0,1,2])
            sage: vp = array('l', [1,2,0])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            Traceback (most recent call last):
            ...
            ValueError: fev relation not satisfied

            sage: fp = array('l', [1,2,0])
            sage: ep = array('l', [1,2,0])
            sage: vp = array('l', [1,2,0])
            sage: Triangulation.from_face_edge_perms(fp, ep, vp)
            Traceback (most recent call last):
            ...
            ValueError: broken edge permutation
        """
        n = self._n

        if not perm_check(self._fp, n):
            raise error('fp is not a permutation')
        if not perm_check(self._ep, n):
            raise error('ep is not permutation')
        if not perm_check(self._vp, n):
            raise error('vp is not a permutation')
        if not perms_are_transitive([self._fp, self._ep, self._vp]):
            raise error('(fp, ep, vp) do not generate a transitive group')

        # The face, edge, and vertex permutations fp, ep, and vp must
        # satisfy the relations fp^3 = ep^2 = fp.ep.vp = Id.  Note
        # that this is a representation of PSL(2, \ZZ) \isom S^2(2, 3,
        # \infty) into the symmetric group on the half edges.

        for i in range(n):
            # NOTE: this first test is relevant only if we enforce triangulations
            # when we generalize to CW complex we can simply remove it
            if self._fp[i] == i or self._fp[self._fp[self._fp[i]]] != i:
                raise error('broken face permutation')
            if self._ep[self._ep[i]] != i:
                raise error('broken edge permutation')
            if self._fp[self._ep[self._vp[i]]] != i:
                raise error('fev relation not satisfied')

    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n == other._n and self._fp == other._fp and self._ep == other._ep

    def __ne__(self, other):
        if type(self) != type(other):
            raise TypeError
        return self._n != other._n or self._fp != other._fp or self._ep != other._ep

    def copy(self):
        r"""
        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation([[0,1,2],[-1,-2,-3]])
            sage: U = T.copy()
            sage: T == U
            True
            sage: T.flip(0)
            sage: T == U
            False
        """
        T = Triangulation.__new__(Triangulation)
        T._n = self._n
        T._fp = self._fp[:]
        T._ep = self._ep[:]
        T._vp = self._vp[:]
        return T

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
        require_package('flipper', 'to_flipper')
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
        require_package('curver', 'to_curver')
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

    def face_permutation(self, copy=True):
        if copy:
            return self._fp[:]
        else:
            return self._fp

    def edge_permutation(self, copy=True):
        if copy:
            return self._ep[:]
        else:
            return self._ep

    def vertex_permutation(self, copy=True):
        if copy:
            return self._vp[:]
        else:
            return self._vp

    def num_half_edges(self):
        return self._n

    def folded_edges(self):
        n = self._n
        ep = self._ep
        return [i for i in range(n) if ep[i]==i]

    def num_folded_edges(self):
        n = self._n
        ep = self._ep
        return sum(ep[i] == i for i in range(n))

    def num_edges(self):
        return (self._n + self.num_folded_edges()) // 2

    def _edge_rep(self, e):
        f = self._ep[e]
        if f < e: return '~%d' % f
        else: return str(e)

    def _norm(self, e):
        f = self._ep[e]
        return f if f < e else e

    def num_faces(self):
        return perm_num_cycles(self._fp)

    def faces(self):
        r"""
        Return the list of faces as triples of half-edges

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.faces()
            [[0, 1, 2], [3, 4, 5], [6, 8, 7]]
        """
        return perm_cycles(self._fp)

    def edges(self):
        r"""
        Return the list of faces as tuples of half-edges

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.edges()
            [[0, 8], [1], [2], [3, 7], [4], [5], [6]]
        """
        return perm_cycles(self._ep)

    def vertices(self):
        r"""
        Return the list of faces as tuples of half-edges

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: T.vertices()
            [[0, 2, 1, 8, 6, 3, 5, 4, 7]]
        """
        return perm_cycles(self._vp)

    def num_vertices(self):
        return perm_num_cycles(self._vp)

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
        """
        return self.num_faces() - self.num_edges() + (self.num_vertices() + self.num_folded_edges())

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
        """
        # chi = 2 - 2g
        return (2 - self.euler_characteristic()) // 2

    def __str__(self):
        r"""
        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: str(T)
            'Triangulation("(0,1,2)(~2,~0,~1)")'
        """
        return 'Triangulation("%s")' % perm_cycle_string(self._fp, n=self._n, involution=self._ep)

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
            sage: T.homology_matrix()
            [ 1  0  0]
            [ 0  1  0]
            [ 1  0  0]
            [ 0  0  1]
            [ 1 -1  0]
            [ 1  0 -1]

            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,6,0)")
            sage: T.homology_matrix()
            [ 1  0]
            [ 0  1]
            [ 1  0]
            [ 0  0]
            [ 1 -1]
            [ 1  0]
            [ 0  0]
        """
        require_package('sage', 'homology_matrix')
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ

        ep = self._ep
        nf = self.num_faces()
        ne = self.num_edges()
        nfe = self.num_folded_edges()
        m = matrix(ZZ, nf + nfe, ne)

        # face equations
        for i,f in enumerate(self.faces()):
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

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
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

            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)")
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

        a,b,c,d = self.square_about_edge(e)
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
                m[e] = V = m._row_ambient_module().zero()
            elif a == A:
                m[e] = m[d] if d < D else -m[D]
            elif d == D:
                m[e] = m[a] if a < A else -m[A]
            else:
                m[e] = (m[a] if a < A else -m[A]) + (m[d] if d < D else -m[D])

    def relabel_homological_action(self, p, m, twist=False):
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

            sage: T = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)")
            sage: u = [0,1,0,0,-1,0]
            sage: v = [0,0,0,1,0,-1]
            sage: w = [1,1,1,1,0,0]
            sage: A = matrix(ZZ, 3, 6, [u,v,w]).transpose()
            sage: T._check_homology_matrix(A)
            sage: perms = ["(0,1)(~0,~1)", "(3,~3)", "(3,~5)(~3,5)",
            ....:          "(0,1,2,~0,~1,~2)",
            ....:          "(0,~1,4,~0,1,~4)(2,3)(~2,~3)"]
            sage: for p in perms*5:
            ....:     T.relabel_homological_action(p, A)
            ....:     T.relabel(p)
            ....:     T._check_homology_matrix(A)
        """
        n = self._n
        ne = self.num_edges()
        if not perm_check(p, n):
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
                e,E = E,e
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

    def is_flippable(self, e):
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
        """
        e = int(e)
        E = self._ep[e]
        a = self._fp[e]
        b = self._fp[a]
        return a != E and b != E

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
        """
        n = self._n
        ep = self._ep
        return [e for e in range(n) if e <= ep[e] and self.is_flippable(e)]

    def square_about_edge(self, e):
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

        e = int(e)
        E = self._ep[e]

        a = self._fp[e]
        b = self._fp[a]
        c = self._fp[E]
        d = self._fp[c]

        return a,b,c,d

    def swap(self, e):
        r"""
        Change the orientation of the edge ``e``.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.swap(0)
            sage: T
            Triangulation("(0,~1,~2)(1,2,~0)")
            sage: T.swap(1)
            sage: T
            Triangulation("(0,1,~2)(2,~0,~1)")
            sage: T.swap(2)
            sage: T
            Triangulation("(0,1,2)(~2,~0,~1)")

            sage: T = Triangulation("(0,~5,4)(3,5,6)(1,2,~6)")
            sage: T.swap(0)
            sage: T
            Triangulation("(0,~5,4)(1,2,~6)(3,5,6)")
            sage: T.swap(5)
            sage: T
            Triangulation("(0,5,4)(1,2,~6)(3,~5,6)")

        Also works for veering triangulations::

            sage: from veerer import VeeringTriangulation

            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: V = VeeringTriangulation(fp, cols)
            sage: V.swap(0)
            sage: V.swap(10)
            sage: V
            VeeringTriangulation("(0,1,~3)(2,~0,~1)(3,4,~5)(5,~9,~7)(6,~2,~4)(7,~6,8)(9,~10,~11)(10,11,~8)", "BRBBBRRBBBBR")

        One can alternatively use ``relabel``::

            sage: T = Triangulation("(0,~5,4)(3,5,6)(1,2,~6)")
            sage: T1 = T.copy()
            sage: T1.swap(3)
            sage: T1.swap(6)
            sage: T2 = T.copy()
            sage: T2.relabel("(3,~3)(6,~6)")
            sage: T1 == T2
            True
        """
        vp = self._vp
        ep = self._ep
        fp = self._fp
        E = ep[e]

        if e == E:
            return

        # images/preimages by vp
        e_vp = vp[e]
        E_vp = vp[E]
        e_vp_inv = fp[E]
        E_vp_inv = fp[e]
        assert vp[e_vp_inv] == e
        assert vp[E_vp_inv] == E

        # images/preimages by fp
        e_fp = fp[e]
        E_fp = fp[E]
        e_fp_inv = ep[e_vp]
        E_fp_inv = ep[E_vp]
        assert fp[e_fp_inv] == e
        assert fp[E_fp_inv] == E

        fp[e_fp_inv] = E
        fp[E_fp_inv] = e
        vp[e_vp_inv] = E
        vp[E_vp_inv] = e
        fp[e] = E_fp
        fp[E] = e_fp
        vp[e] = E_vp
        vp[E] = e_vp

        # TODO: remove check
        self._check()

    def relabel(self, p):
        r"""
        Relabel this triangulation inplace according to the permutation ``p``.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.relabel("(0,~0)"); T
            Triangulation("(0,~1,~2)(1,2,~0)")
            sage: T.relabel("(0,1,~2)"); T
            Triangulation("(0,1,2)(~2,~0,~1)")

        An example of a flip sequence which forms a loop after non-trivial relabelling::

            sage: T0 = Triangulation("(1,~0,4)(2,~4,~1)(3,~2,5)(~5,~3,0)")
            sage: T = T0.copy()
            sage: T.flip_back(1)
            sage: T.flip_back(3)
            sage: T.flip_back(0)
            sage: T.flip_back(2)
            sage: T.relabel("(0,2)(1,3)")
            sage: T == T0
            True
        """
        n = self._n
        if not perm_check(p, n):
            # if the input is not a valid permutation, we assume that half-edges
            # are not separated
            p = perm_init(p, self._n, self._ep)
            if not perm_check(p, n, self._ep):
                raise ValueError('invalid relabeling permutation')

        self._vp = perm_conjugate(self._vp, p)
        self._ep = perm_conjugate(self._ep, p)
        self._fp = perm_conjugate(self._fp, p)


    def flip(self, e):
        r"""
        Flip the edge ``e``.

        EXAMPLES::

            sage: from veerer import Triangulation

        Flipping a folded edge is an involution::

            sage: T = Triangulation("(0,1,2)")
            sage: T
            Triangulation("(0,1,2)")
            sage: T.flip(0); T
            Triangulation("(0,2,1)")
            sage: T.flip(0); T
            Triangulation("(0,1,2)")

            sage: T == Triangulation("(0,1,2)")
            True
        """
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

        e = int(e)
        E = self._ep[e]

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

    def flip_back(self, e):
        r"""
        Flip back the edge ``e``.

        EXAMPLES::

            sage: from veerer import *
            sage: T0 = Triangulation([(0,1,2),(-1,-2,-3)])
            sage: T = T0.copy()
            sage: T.flip(0)
            sage: T.flip_back(0)
            sage: T == T0
            True

            sage: T.flip(1); T.flip(2)
            sage: T.flip_back(2); T.flip_back(1)
            sage: T == T0
            True
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

        e = int(e)
        E = self._ep[e]

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

    def to_string(self):
        r"""
        Serialize this triangulation as a string.

        EXAMPLES::

            sage: from veerer import *

            sage: T = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: T.to_string()
            '6_120534_543210'
        """
        return uint_base64_str(self._n) + '_' + perm_base64_str(self._fp) + '_' + perm_base64_str(self._ep)

    def conjugate(self):
        r"""
        Conjugate this triangulation.

        The face permutation is replaced by its inverse, conjugated by
        the edge permutation. The vertex permutation is replaced by
        its inverse.

        EXAMPLES::

            sage: from veerer import Triangulation

            sage: T = Triangulation("(0,1,2)(~0,~4,~2)(3,4,5)(~3,~1,~5)")
            sage: T.conjugate()
            sage: T
            Triangulation("(0,2,4)(1,3,5)(~5,~4,~3)(~2,~1,~0)")
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

        self._fp = perm_conjugate(perm_invert(self._fp), self._ep)
        self._vp = perm_invert(self._vp)

    def _relabelling_from(self, start_edge):
        r"""
        Return a canonical relabelling map obtained from walking
        along the triangulation starting at ``start_edge``.

        The returned relabelling array maps the current edge to the new
        labelling.

        EXAMPLES::

            sage: from veerer import *
            sage: from array import array

        The torus example (6 symmetries)::

            sage: fp = array('l', [1, 5, 4, 2, 3, 0])
            sage: ep = array('l', [4, 3, 5, 1, 0, 2])
            sage: vp = array('l', [2, 4, 1, 0, 5, 3])
            sage: T = Triangulation.from_face_edge_perms(fp, ep, vp)
            sage: T._relabelling_from(3)
            array('l', [4, 5, 3, 0, 1, 2])

            sage: p = T._relabelling_from(0)
            sage: T.relabel(p)
            sage: for i in range(6):
            ....:     p = T._relabelling_from(i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     assert S == T

        The sphere example (3 symmetries)::

            sage: fp = array('l', [1, 2, 0])
            sage: ep = array('l', [0, 1, 2])
            sage: vp = array('l', [2, 0, 1])
            sage: T = Triangulation.from_face_edge_perms(fp, ep, vp)
            sage: T._relabelling_from(1)
            array('l', [1, 0, 2])
            sage: p = T._relabelling_from(0)
            sage: T.relabel(p)
            sage: for i in range(3):
            ....:     p = T._relabelling_from(i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     assert S == T

        An example with no automorphism::

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)")
            sage: p = T._relabelling_from(0)
            sage: T.relabel(p)
            sage: for i in range(1, 9):
            ....:     p = T._relabelling_from(i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     S._check()
            ....:     assert S != T
        """
        n = self._n
        ep = self._ep
        vp = self._vp

        k = 0     # current available label at the front.
        m = n - 1 # current available label at the back.
        relabelling = array('l', [-1] * n)
        relabelling[start_edge] = 0
        k = k + 1

        if ep[start_edge] != start_edge:
            relabelling[ep[start_edge]] = m
            m = m - 1

        to_process = [start_edge]
        if ep[start_edge] != start_edge:
            to_process.append(ep[start_edge])

        while to_process:
            e0 = to_process.pop()
            e = vp[e0]
            while e != e0:
                if relabelling[e] == -1:
                    relabelling[e] = k
                    k = k + 1
                    if ep[e] != e:
                        relabelling[ep[e]] = m
                        m = m - 1
                        to_process.append(ep[e])
                e = vp[e]

        # check that everybody has a name!
        assert k == m + 1

        return relabelling

    def _automorphism_good_starts(self):
        # we discriminate based on the lengths of cycles in
        #   vp, ep, fp and vp[ep[fp]]
        n = self._n

        lv = [-1] * n
        for c in perm_cycles(self._vp):
            m = len(c)
            for i in c:
                lv[i] = m

        le = [-1] * n
        for c in perm_cycles(self._ep):
            m = len(c)
            for i in c:
                le[i] = m

        lf = [-1] * n
        for c in perm_cycles(self._fp):
            m = len(c)
            for i in c:
                lf[i] = m

        lvef = [-1] * n
        vef = perm_compose(perm_compose(self._fp, self._ep), self._vp)
        for c in perm_cycles(vef):
            m = len(c)
            for i in c:
                lvef[i] = m

        d = {}
        for i in range(n):
            s = (lv[i], le[i], lf[i], lvef[i])
            if s in d:
                d[s].append(i)
            else:
                d[s] = [i]
        m = min(len(x) for x in d.values())
        candidates = [s for s,v in d.items() if len(v) == m]
        winner = min(candidates)

        return d[winner]

    def automorphisms(self):
        r"""
        Return the list of automorphisms of this triangulation.

        The output is a list of arrays that are permutations acting on the set
        of half edges.

        EXAMPLES::

            sage: from veerer import *

        An example with 4 symmetries in genus 2::

            sage: T = Triangulation("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)")
            sage: A = T.automorphisms()
            sage: len(A)
            4

        And the "sphere octagon" has 8::

            sage: s  = "(0,8,~7)(1,9,~0)(2,10,~1)(3,11,~2)(4,12,~3)(5,13,~4)(6,14,~5)(7,15,~6)"
            sage: len(Triangulation(s).automorphisms())
            8
        """
        n = self._n
        fp = self._fp
        ep = self._ep

        best = None
        best_relabellings = []

        for start_edge in self._automorphism_good_starts():
            relabelling = self._relabelling_from(start_edge)

            fp_new = perm_conjugate(fp, relabelling)
            ep_new = perm_conjugate(ep, relabelling)

            T = (fp_new, ep_new)
            if best is None or T == best:
                best_relabellings.append(relabelling)
                best = T
            elif T < best:
                del best_relabellings[:]
                best_relabellings.append(relabelling)
                best = T

        p0 = perm_invert(best_relabellings[0])
        return [perm_compose(p, p0) for p in best_relabellings]

    def best_relabelling(self):
        n = self._n
        fp = self._fp
        ep = self._ep

        best = None

        for start_edge in self._automorphism_good_starts():
            relabelling = self._relabelling_from(start_edge)

            fp_new = perm_conjugate(fp, relabelling)
            ep_new = perm_conjugate(ep, relabelling)

            T = (fp_new, ep_new)
            if best is None or T < best:
                best_relabelling = relabelling
                best = T

        return best_relabelling, best

    def iso_sig(self):
        r"""
        Return a canonical signature for this triangulation.

        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(0,3,1)(~0,4,2)(~1,~2,~4)")
            sage: T.iso_sig()
            '9_345876210_087654321'
            sage: TT = Triangulation.from_string(T.iso_sig())
            sage: TT
            Triangulation("(0,3,~1)(1,4,~2)(2,~4,~3)")
            sage: TT.iso_sig() == T.iso_sig()
            True

            sage: T = Triangulation("(0,10,~6)(1,12,~2)(2,14,~3)(3,16,~4)(4,~13,~5)(5,~1,~0)(6,~17,~7)(7,~14,~8)(8,13,~9)(9,~11,~10)(11,~15,~12)(15,17,~16)")
            sage: T.iso_sig()
            'A_lkp9ijfmcvegq0osnyutxdrbaw87654321zh_zyxwvutsrqponmlkjihgfedcba9876543210'
            sage: Triangulation.from_string(T.iso_sig())
            Triangulation("(0,~14,13)(1,~15,~2)(2,~10,~3)(3,9,~4)(4,~17,~5)(5,~16,~6)(6,15,~7)(7,~13,~8)(8,12,~9)(10,14,~11)(11,16,~12)(17,~1,~0)")
        """
        n = self._n
        _, (fp, ep) = self.best_relabelling()

        fp = perm_base64_str(fp)
        ep = perm_base64_str(ep)

        return uint_base64_str(n) + '_' + fp + '_' + ep

    def is_isomorphic_to(self, other, certificate=False):
        r"""
        Check wheter ``self`` is isomorphic to ``other``.

        TESTS::

            sage: from veerer import Triangulation
            sage: from veerer.permutation import perm_random_centralizer
            sage: T = Triangulation("(0,5,1)(~0,4,2)(~1,~2,~4)(3,6,~5)")
            sage: TT = T.copy()
            sage: for _ in range(10):
            ....:     rel = perm_random_centralizer(TT.edge_permutation(False))
            ....:     TT.relabel(rel)
            ....:     assert T.is_isomorphic_to(TT)
        """
        if type(self) is not type(other):
            raise TypeError("can only check isomorphisms between two triangulations")

        if self._n != other._n or self.num_folded_edges() != other.num_folded_edges():
            return (False, None) if certificate else False

        r1, data1 = self.best_relabelling()
        r2, data2 = other.best_relabelling()

        if data1 != data2:
            return (False, None) if certificate else False
        elif certificate:
            return (True, perm_compose(r1, perm_invert(r2)))
        else:
            return True

    def cover(self, c):
        from .cover import TriangulationCover
        return TriangulationCover(self, c)
