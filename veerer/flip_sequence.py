r"""
Dynamical forward flip sequences (and relabeling) in veering triangulations.
"""
######################################################################
# This file is part of veering.
#
#       Copyright (C) 2020 Vincent Delecroix
#
# veerer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# veerer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with veerer. If not, see <https://www.gnu.org/licenses/>.
######################################################################

from __future__ import absolute_import

from .constants import colour_from_char, colour_to_char
from .permutation import perm_init, perm_check, perm_id, perm_is_one, perm_preimage, perm_invert, perm_cycle_string, perm_compose
from .veering_triangulation import VeeringTriangulation

def flip_sequence_to_string(sequence):
    return " ".join("%d%s" % (e, colour_to_char(col)) for e,col in sequence)

def flip_sequence_from_string(s):
    return [(int(f[:-1]), colour_from_char(f[-1])) for f in s.split()]

class VeeringFlipSequence(object):
    r"""
    A dynamical sequence of forward flips followed by a relabelling.

    EXAMPLES::

        sage: from veerer import VeeringTriangulation, VeeringFlipSequence, BLUE, RED
        sage: T = VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB")
        sage: F = VeeringFlipSequence(T)
        sage: F
        VeeringFlipSequence(VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB"), (0)(1)(2)(~2)(~1)(~0))
        sage: F.flip(1, RED)
        sage: F.flip(0, RED)
        sage: F
        VeeringFlipSequence(VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB"), "1R 0R", (0)(1)(2)(~2)(~1)(~0))

    The flips can also be specified in the input as a string or as a list::

        sage: VeeringFlipSequence(T, "1R 0R")
        VeeringFlipSequence(VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB"), "1R 0R", (0)(1)(2)(~2)(~1)(~0))
        sage: VeeringFlipSequence(T, [(1, RED), (0, RED)])
        VeeringFlipSequence(VeeringTriangulation("(0,1,2)(~2,~0,~1)", "RRB"), "1R 0R", (0)(1)(2)(~2)(~1)(~0))
    """
    def __init__(self, start, sequence=None, relabelling=None, reduced=False):
        if not isinstance(start, VeeringTriangulation):
            raise TypeError("'start' must be a VeeringTriangulation")
        self._start = start
        self._end = start.copy()
        self._relabelling = perm_id(self._start._n)
        self._flips = []   # list of triples (e, col_after, col_before)
        self._reduced = reduced
        if reduced:
            self._start.forgot_forward_flippable_colour()
        if sequence is not None:
            if isinstance(sequence, str):
                sequence = flip_sequence_from_string(sequence)
            for e, col in sequence:
                self.flip(e, col)
        if relabelling is not None:
            self.relabel(relabelling)

    def _check(self):
        T = VeeringFlipSequence(self._start, [f[:2] for f in self._flips], self._relabelling, self._reduced)
        if T != self:
            raise ValueError

    def __repr__(self):
        r"""
        TESTS::

            sage: from veerer import VeeringTriangulation, VeeringFlipSequence

            sage: T = VeeringTriangulation("(0,1,2)(~1,~2,~0)", "RRB")
            sage: VeeringFlipSequence(T, "1R 0R")
            VeeringFlipSequence(VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RRB"), "1R 0R", (0)(1)(2)(~0)(~2)(~1))
            sage: VeeringFlipSequence(T, "1R 0R", reduced=True)
            VeeringFlipSequence(VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RPB"), "1R 0R", (0)(1)(2)(~0)(~2)(~1), reduced=True)
            sage: VeeringFlipSequence(T, "1R 0R", relabelling="(0,5)(1,3)")
            VeeringFlipSequence(VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RPB"), "1R 0R", (0,~1)(1,~0)(2)(~2))
            sage: VeeringFlipSequence(T, "1R 0R", relabelling="(0,5)(1,3)", reduced=True)
            VeeringFlipSequence(VeeringTriangulation("(0,1,2)(~0,~1,~2)", "RPB"), "1R 0R", (0,~1)(1,~0)(2)(~2), reduced=True)
        """
        args = [repr(self._start)]
        if self._flips:
            args.append("\"%s\"" % flip_sequence_to_string(self.flips()))
        if self._relabelling is not None:
            args.append(perm_cycle_string(self._relabelling, self._end._n, involution=self._end._ep))
        if self._reduced:
            args.append("reduced=True")
        return "VeeringFlipSequence({})".format(", ".join(args))

    # properties
    def __eq__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._start == other._start and \
               self._end == other._end and \
               self._flips == other._flips and \
               self._reduced == other._reduced and \
               self._relabelling == other._relabelling

    def __ne__(self, other):
        return not (self == other)

    def copy(self):
        F = VeeringFlipSequence.__new__(VeeringFlipSequence)
        F._start = self._start.copy()
        F._end = self._end.copy()
        F._flips = self._flips[:]
        F._reduced = self._reduced
        F._relabelling = self._relabelling
        return F

    def matrix(self, twist=True):
        from sage.matrix.constructor import matrix
        m = matrix(ZZ, self._start._n)
        V = self._start.copy()
        for e, col, _ in self._flips:
            V.flip_homological_action(e, m, twist)
            V.flip(e, col)
        if self._relabelling is not None:
            V.relabel_homological_action(self._relabelling, m, twist)
        return m

    def flips(self):
        return [flip[:2] for flip in self._flips]

    def is_closed(self):
        return self._start == self._end

    def is_pseudo_anosov(self):
        r"""
        Test whether the flip sequence defines a pseudo-Anosov mapping class.

        EXAMPLES::

            sage: from veerer import VeeringTriangulation, VeeringFlipSequence
            sage: F2 = VeeringFlipSequence(VeeringTriangulation("(0,3,4)(1,~3,5)(2,6,~4)", "PPPBRRB"), "0B 1B", reduced=True)
            sage: F3 = VeeringFlipSequence(VeeringTriangulation("(0,4,3)(1,5,~3)(2,6,~4)", "BBPPRRB"), "3B", "(0,1)", reduced=True)
            sage: F4 = VeeringFlipSequence(VeeringTriangulation("(0,4,3)(1,5,~3)(2,6,~4)", "BBPPRRB"), "2R 3R", "(0,6,1)(2,5)(3,4)", reduced=True)
            sage: F6 = VeeringFlipSequence(VeeringTriangulation("(0,4,3)(1,5,~3)(2,6,~4)", "BBPPRRB"), "2B", "(2,6)", reduced=True)

        Pseudo-Anosov examples::

            sage: assert (F2 * F4 * F3).is_pseudo_anosov()
            sage: assert (F4 * F6).is_pseudo_anosov()
            sage: assert (F4 * F4 * F6).is_pseudo_anosov()

        Non pseudo-Anosov::

            sage: assert not F2.is_pseudo_anosov()
            sage: assert not (F2 * F3).is_pseudo_anosov()
            sage: assert not (F3 * F2).is_pseudo_anosov()
            sage: assert not (F4 * F4 * F4 * F6).is_pseudo_anosov()
        """
        if self._start != self._end:
            return False
        n = self._start.num_edges()
        flipped = set(x[0] for x in self._flips)
        if len(flipped) == n:
            return True

        new = list(flipped)
        very_new = []
        modified = True
        r = self._relabelling
        while len(flipped) < n and modified:
            modified = False
            for e in new:
                e = r[e]
                if e not in flipped:
                    modified = True
                    flipped.add(e)
                    very_new.append(e)
            new, very_new = very_new, new
            del very_new[:]

        return len(flipped) == n

    def self_similar_surface(self):
        pass

    # change
    def flip(self, e, col):
        oldcol = self._end._colouring[e]
        E = self._end._ep[e]
        if E < e:
            e = E
        self._end.flip(e, col, reduced=self._reduced)
        if self._relabelling is not None:
            # push the flip to the left of relabelling
            e = perm_preimage(self._relabelling, e)
        self._flips.append((e, col, oldcol))

    def relabel(self, r):
        end = self._end
        if not perm_check(r, end._n):
            r = perm_init(r, end._n, end._ep)
            if not perm_check(r, end._n):
                raise ValueError('invalid relabelling permutation')

        end.relabel(r)
        self._relabelling = perm_compose(self._relabelling, r)

    def __imul__(self, other):
        r"""
        EXAMPLES::

            sage: from veerer import VeeringTriangulation, VeeringFlipSequence

            sage: V0 = VeeringTriangulation("(0,3,4)(1,~3,5)(2,6,~4)", "PPPBRRB")
            sage: F = VeeringFlipSequence(V0, "2B", "(2,6)", reduced=True)
            sage: assert F.is_closed()       
            sage: F *= F
            sage: F
            VeeringFlipSequence(VeeringTriangulation("(0,3,4)(1,~3,5)(2,6,~4)", "PPPBRRB"), "2B 6B", (0)(1)(2)(3)(4)(5)(6)(~4)(~3), reduced=True)
        """
        if type(self) != type(other):
            raise TypeError
        if self._end != other._start or self._reduced != other._reduced:
            raise ValueError("composition undefined")

        r = perm_invert(self._relabelling, self._start._n)
        self._flips.extend([(r[e], col, oldcol) for e, col, oldcol in other._flips])
        self._end = other._end.copy()
        self._relabelling = perm_compose(self._relabelling, other._relabelling)

        # TODO: remove this expensive check!!
        self._check()
        return self

    def __mul__(self, other):
        r"""
        TESTS::

            sage: from veerer import VeeringTriangulation, VeeringFlipSequence
            sage: V2 = VeeringTriangulation("(0,3,4)(1,~3,5)(2,6,~4)", "PPPBRRB")
            sage: V3 = VeeringTriangulation("(0,4,3)(1,5,~3)(2,6,~4)", "BBPPRRB")
            sage: F2 = VeeringFlipSequence(V2, "0B 1B", reduced=True)
            sage: F3 = VeeringFlipSequence(V3, "3B", "(0,1)", reduced=True)
            sage: (F2 * F3) * (F2 * F3) == (F2 * (F3 * F2)) * F3 == F2 * ((F3 * F2) * F3)
            True
        """
        res = self.copy()
        res *= other
        return res
