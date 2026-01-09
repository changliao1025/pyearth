import numpy as np


class pynvector(object):
    """
    The vector class

    Args:
        object (_type_): _description_

    Returns:
        _type_: _description_
    """

    def __init__(self, aParameter):
        self.dX = float(aParameter.get("x", 0.0))
        self.dY = float(aParameter.get("y", 0.0))
        self.dZ = float(aParameter.get("z", 0.0))

        length = self.length()
        if length > 0:
            self.dX /= length
            self.dY /= length
            self.dZ /= length
        else:
            self.dX = self.dY = self.dZ = 0.0  # Normalize zero vector to zero

    def length(self):
        return np.sqrt(self.dX * self.dX + self.dY * self.dY + self.dZ * self.dZ)

    def dot(self, other):
        """Calculate dot product with another nvector."""
        if not isinstance(other, pynvector):
            raise TypeError(f"Cannot compute dot product with {type(other).__name__}")
        return self.dX * other.dX + self.dY * other.dY + self.dZ * other.dZ

    def __mul__(self, scalar):
        """Scalar multiplication of nvector."""
        return pynvector(
            {"x": self.dX * scalar, "y": self.dY * scalar, "z": self.dZ * scalar}
        )

    def __rmul__(self, scalar):
        """Right scalar multiplication of nvector."""
        return self.__mul__(scalar)

    def normalize(self):
        """Normalize the nvector to unit length."""
        length = self.length()
        if length > 0:
            return pynvector(
                {"x": self.dX / length, "y": self.dY / length, "z": self.dZ / length}
            )
        else:
            return pynvector({"x": 0.0, "y": 0.0, "z": 0.0})

    def __add__(self, other):
        if not isinstance(other, pynvector):
            return NotImplemented
        return pynvector(
            {"x": self.dX + other.dX, "y": self.dY + other.dY, "z": self.dZ + other.dZ}
        )

    def __repr__(self):
        return f"pynvector(x={self.dX}, y={self.dY}, z={self.dZ})"

    def __str__(self):
        return f"[{self.dX:.4f}, {self.dY:.4f}, {self.dZ:.4f}]"

    def toLatLon(self):
        from pyearth.toolbox.mesh.point import pypoint

        lat_rad = np.arctan2(self.dZ, np.sqrt(self.dX**2 + self.dY**2))
        lon_rad = np.arctan2(self.dY, self.dX)
        return pypoint(
            {
                "dLongitude_degree": np.degrees(lon_rad),
                "dLatitude_degree": np.degrees(lat_rad),
            }
        )
