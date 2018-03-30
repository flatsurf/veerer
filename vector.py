
import math
from constants import *

EPSILON = 10**-8
CUTOFF = 1.0
CUTOFF2 = CUTOFF**2

class Vector2(object):
    # Warning: This needs to be updated if the interals of this class ever change.
    __slots__ = ['x', 'y']
    def __init__(self, x, y):
        self.x, self.y = x, y
    def __str__(self):
        return "(%s, %s)" % (self.x, self.y)
    def __repr__(self):
        return str(self)
    def __reduce__(self):
        # Having __slots__ means we need to pickle manually.
        return (self.__class__, (self.x, self.y))
    def approx(self, other, epsilon=EPSILON):
        return (self - other).norm2() < epsilon
    def __eq__(self, other):
        raise TypeError('Susceptible to FPE.')
        return self.x == other.x and self.y == other.y
    def __ne__(self, other):
        return not (self == other)
    def __neg__(self):
        return Vector2(-self.x, -self.y)
    def __add__(self, other):
        return Vector2(self.x + other.x, self.y + other.y)
    def __radd__(self, other):
        return self + other
    def __sub__(self, other):
        return self + -other
    def __mul__(self, other):
        return Vector2(self.x * other, self.y * other)
    def __rmul__(self, other):
        return self * other
    def dot(self, other):
        return self.x * other.x + self.y * other.y
    def norm2(self):
        return self.dot(self)
    def norm(self):
        return math.sqrt(self.norm2())
    def colour(self):
        return RED if self.x * self.y > 0 else BLUE
    def slope(self):
        return HORIZONTAL if self.x**2 > self.y**2 else VERTICAL
    def flow(self, time):
        return Vector2(self.x * math.exp(time / 2), self.y * math.exp(-time / 2))
    def short_times(self, clock):
        disc = CUTOFF**4 - 4 * self.x**2 * self.y**2
        if disc < 0:
            return float('inf'), float('-inf')
        else:
            return math.log(CUTOFF**2 - math.sqrt(disc)) - math.log(2 * self.x**2) + clock, math.log(CUTOFF**2 + math.sqrt(disc)) - math.log(2 * self.x**2) + clock
    def min_length2(self):
        return 2 * abs(self.x) * abs(self.y)
    def min_time(self, clock):
        return math.log(abs(self.y)) - math.log(abs(self.x)) + clock  # = log(aspect ratio).

