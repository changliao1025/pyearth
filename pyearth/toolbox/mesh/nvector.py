import numpy as np


class pynvector(object):
    """
    N-vector representation for geographic points on a sphere.

    Supports high-precision mode using numpy.float128 for improved
    accuracy in spherical interpolation operations.

    Args:
        aParameter (dict): Dictionary with 'x', 'y', 'z' coordinates
        use_high_precision (bool): Use float128 for calculations (default: False)

    Returns:
        pynvector: Normalized n-vector
    """

    def __init__(self, aParameter, use_high_precision=False):
        self.use_high_precision = use_high_precision
        dtype = np.float128 if use_high_precision else np.float64

        # Extract and convert to appropriate precision
        self.dX = dtype(aParameter.get("x", 0.0))
        self.dY = dtype(aParameter.get("y", 0.0))
        self.dZ = dtype(aParameter.get("z", 0.0))

        # Normalize to unit length
        length = self.length()
        if length > 0:
            self.dX /= length
            self.dY /= length
            self.dZ /= length
        else:
            self.dX = self.dY = self.dZ = dtype(0.0)  # Normalize zero vector to zero

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
            {"x": self.dX * scalar, "y": self.dY * scalar, "z": self.dZ * scalar},
            use_high_precision=self.use_high_precision
        )

    def __rmul__(self, scalar):
        """Right scalar multiplication of nvector."""
        return self.__mul__(scalar)

    def normalize(self):
        """Normalize the nvector to unit length."""
        length = self.length()
        if length > 0:
            return pynvector(
                {"x": self.dX / length, "y": self.dY / length, "z": self.dZ / length},
                use_high_precision=self.use_high_precision
            )
        else:
            dtype = np.float128 if self.use_high_precision else np.float64
            return pynvector({"x": dtype(0.0), "y": dtype(0.0), "z": dtype(0.0)},
                           use_high_precision=self.use_high_precision)

    def __add__(self, other):
        if not isinstance(other, pynvector):
            return NotImplemented
        return pynvector(
            {"x": self.dX + other.dX, "y": self.dY + other.dY, "z": self.dZ + other.dZ},
            use_high_precision=self.use_high_precision
        )

    def __repr__(self):
        return f"pynvector(x={self.dX}, y={self.dY}, z={self.dZ})"

    def __str__(self):
        return f"[{self.dX:.4f}, {self.dY:.4f}, {self.dZ:.4f}]"

    def toLatLon(self):
        """
        Convert n-vector back to latitude/longitude coordinates.

        Uses high precision (float128) for intermediate calculations if enabled,
        but always returns float64 coordinates for pypoint storage.

        Returns:
            pypoint: Point with latitude and longitude in degrees
        """
        from pyearth.toolbox.mesh.point import pypoint

        # Perform calculations in high precision if enabled
        if self.use_high_precision:
            # Ensure we're working with float128
            x = np.float128(self.dX)
            y = np.float128(self.dY)
            z = np.float128(self.dZ)
            lat_rad = np.arctan2(z, np.sqrt(x**2 + y**2))
            lon_rad = np.arctan2(y, x)
        else:
            lat_rad = np.arctan2(self.dZ, np.sqrt(self.dX**2 + self.dY**2))
            lon_rad = np.arctan2(self.dY, self.dX)

        # Always return as float64 for storage
        return pypoint(
            {
                "dLongitude_degree": float(np.degrees(lon_rad)),
                "dLatitude_degree": float(np.degrees(lat_rad)),
            }
        )
