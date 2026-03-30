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

    def __init__(self, aParameter, use_high_precision=False, normalize=True):
        self.use_high_precision = use_high_precision
        dtype = np.float128 if use_high_precision else np.float64

        # Extract and convert to appropriate precision
        self.dX = dtype(aParameter.get("x", 0.0))
        self.dY = dtype(aParameter.get("y", 0.0))
        self.dZ = dtype(aParameter.get("z", 0.0))

        # Normalize to unit length unless this is an intermediate linear-combination value.
        if normalize:
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
            use_high_precision=self.use_high_precision,
            normalize=False,
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
                use_high_precision=self.use_high_precision,
                normalize=False,
            )
        else:
            dtype = np.float128 if self.use_high_precision else np.float64
            return pynvector({"x": dtype(0.0), "y": dtype(0.0), "z": dtype(0.0)},
                           use_high_precision=self.use_high_precision,
                           normalize=False)

    def __add__(self, other):
        if not isinstance(other, pynvector):
            return NotImplemented
        return pynvector(
            {"x": self.dX + other.dX, "y": self.dY + other.dY, "z": self.dZ + other.dZ},
            use_high_precision=self.use_high_precision,
            normalize=False,
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

        dtype = np.float128 if self.use_high_precision else np.float64

        # Re-normalize in the target precision to reduce drift before conversion.
        x = dtype(self.dX)
        y = dtype(self.dY)
        z = dtype(self.dZ)
        norm = np.sqrt(x * x + y * y + z * z)

        if norm > dtype(0.0):
            x = x / norm
            y = y / norm
            z = z / norm

        lat_rad = np.arctan2(z, np.sqrt(x * x + y * y))
        lon_rad = np.arctan2(y, x)

        # Always return as float64 for storage
        return pypoint(
            {
                "dLongitude_degree": float(np.degrees(lon_rad)),
                "dLatitude_degree": float(np.degrees(lat_rad)),
            }
        )
