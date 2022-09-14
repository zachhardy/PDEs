import copy
import numpy as np


class CartesianVector:
    """
    Implementation of a 3D spatial coordinate or vector.
    """
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self) -> str:
        return f"({self.x}, {self.y}, {self.z})"

    def __repr__(self):
        return f"CartesianVector{self.__str__()}"

    def __iadd__(self, other):
        if not isinstance(other, type(self)):
            raise TypeError

        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __add__(self, other):
        return copy.deepcopy(self).__iadd__(other)

    def __isub__(self, other):
        if not isinstance(other, type(self)):
            raise TypeError

        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __sub__(self, other):
        return copy.deepcopy(self).__isub__(other)

    def __imul__(self, other):
        if not isinstance(other, (int, float, type(self))):
            raise TypeError

        if isinstance(other, (int, float)):
            self.x *= other
            self.y *= other
            self.z *= other
        elif isinstance(other, type(self)):
            self.x *= other.x
            self.y *= other.y
            self.z *= other.z
        return self

    def __mul__(self, other):
        return copy.deepcopy(self).__imul__(other)

    def __itruediv__(self, other):
        if not isinstance(other, (int, float, type(self))):
            raise TypeError

        if isinstance(other, (int, float)):
            self.x /= other
            self.y /= other
            self.z /= other
        elif isinstance(other, type(self)):
            self.x /= other.x
            self.y /= other.y
            self.z /= other.z
        return self

    def __truediv__(self, other):
        return copy.deepcopy(self).__itruediv__(other)

    def __ipow__(self, power, modulo=None):
        if not isinstance(power, (int, float)):
            raise TypeError

        self.x **= power
        self.y **= power
        self.z **= power
        return self

    def __pow__(self, power, modulo=None):
        return copy.deepcopy(self).__ipow__(power, modulo)

    def __abs__(self):
        return CartesianVector(abs(self.x), abs(self.y), abs(self.z))

    def __neg__(self):
        return CartesianVector(-self.x, -self.y, -self.z)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            raise TypeError

        return (self.x == other.x and
                self.y == other.y and
                self.z == other.z)

    def __ne__(self, other):
        return not self == other

    def dot(self, other):
        if not isinstance(other, type(self)):
            raise TypeError

        dot = self.x * other.x
        dot += self.y * other.y
        dot += self.z * other.z
        return dot

    def cross(self, other):
        if not isinstance(other, type(self)):
            raise TypeError

        x = self.y * other.z - self.z * other.y
        y = self.z * other.x - self.x * other.z
        z = self.x * other.y - self.y * other.x
        return CartesianVector(x, y, z)

    def norm(self):
        return np.sqrt(self.dot(self))

    def direction(self):
        norm = self.norm()
        if norm == 0:
            norm = 1
        return self.__itruediv__(norm)
