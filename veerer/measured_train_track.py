# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2018-2023 Vincent Delecroix
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

from sage.categories.fields import Fields

from sage.structure.sequence import Sequence

from .triangulation import Triangulation

_Fields = Fields()

class MeasuredTrainTrack(object):
    r"""
    Train-track endowed with a transverse measure.
    """
    def __init__(self, triangulation, lengths, base_ring=None):
        self._triangulation = Triangulation(triangulation)

        n = self._triangulation.num_half_edges()
        m = self._triangulation.num_edges()
        ep = self._triangulation.edge_permutation(copy=False)
        if len(lengths) == m:
            for e in range(m):
                E = ep[e]
                if e != E and E < m:
                    raise ValueError("edge perm not in standard form")
            lengths = list(lengths) + [lengths[ep[e]] for e in range(m,n)]
        if len(lengths) != n:
            raise ValueError('wrong number of vectors')


        if base_ring is None:
            lengths = Sequence(lengths)
            base_ring = lengths.universe()
            lengths = list(lengths)
        else:
            lengths = [base_ring.coerce(l) for l in lengths]

        self._base_ring = base_ring
        self._lengths = tuple(lengths)

        fp = self._triangulation.face_permutation(copy=False)

        for e in range(n):
            E = ep[e]
            if self._lengths[e] != self._lengths[E]:
                raise ValueError("non-compatible length data")

        # we record below whether the given size is (horizontally) big or not
        self._edge_type = [None] * n

        for a in range(n):
            b = fp[a]
            c = fp[b]
            if lengths[a] >= lengths[b] and lengths[a] >= lengths[c]:
                if lengths[a] != lengths[b] + lengths[c]:
                    raise ValueError("non-compatible length data")
                self._edge_type[a] = 0
                self._edge_type[b] = 1
                self._edge_type[c] = 2

    def lengths(self):
        return self._lengths

    def __repr__(self):
        return "MeasuredTrainTrack({}, {})".format(
                self._triangulation, self._lengths)

    def __call__(self, p, iterations=1):
        r"""
        p = (i,x) where
        i = half-edge
        x = value

        EXAMPLES::

            sage: from veerer import *  # random output due to deprecation warnings from realalg
            sage: v0 = vector((1, 0, 1, 1))
            sage: v1 = vector((0, 1, 1, 1))
            sage: t = Triangulation("(0,1,2)(~0,~1,3)")
            sage: tt = MeasuredTrainTrack(t, 2*v0 + 3*v1)
            sage: tt((2,3/2))
            (4, 3/2)
            sage: tt((2,3/2),2) == tt((4,3/2))
            True
            sage: tt((2,3/2),5) == tt((4,3/2),4)
            True
        """
        it = self.orbit(p)
        for _ in range(2*iterations):
            next(it)
        return next(it)

    def orbit(self, p):
        r"""
        Return an iterator through the (infinite) orbit of p.

        (intermediate steps are yielded as well)
        """
        n = self._triangulation.num_half_edges()
        ep = self._triangulation.edge_permutation(copy=False)
        fp = self._triangulation.face_permutation(copy=False)
        L = self._lengths
        i, x = p

        while True:
            assert 0 <= i < n
            assert 0 < x < self._lengths[i]
            yield (i,x)

            # 1. cross the tripode (first involution)
            if self._edge_type[i] == 0:
                # big edge
                b = i
                s1 = fp[b]
                s2 = fp[s1]
                if x < L[s2]:
                    i = s2
                    x = L[s2] - x
                else:
                    i = s1
                    x = L[s1] + L[s2] - x
            elif self._edge_type[i] == 1:
                # first small edge
                s1 = i
                s2 = fp[s1]
                b = fp[s2]
                i = b
                x = L[s1] + L[s2] - x
            elif self._edge_type[i] == 2:
                # second small edge
                s2 = i
                b = fp[s2]
                s1 = fp[b]
                i = b
                x = L[s2] - x

            yield (i,x)

            # 2. pass through the edges (second involution)
            i = ep[i]
            x = L[i] - x
