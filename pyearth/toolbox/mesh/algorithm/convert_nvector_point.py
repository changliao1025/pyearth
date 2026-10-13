"""
Adapter functions to convert between pypoint and pynvector representations.
"""

import numpy as np
from pyearth.toolbox.mesh.point import pypoint
from pyearth.toolbox.mesh.nvector import pynvector


_HIGH_PRECISION_DTYPE = getattr(np, "float128", np.longdouble)


def point_to_nvector(pPoint: pypoint, use_high_precision: bool = False) -> pynvector:
    """
    Convert a pypoint object to an n-vector representation.

    Uses high precision for trigonometric calculations when enabled
    to maintain accuracy during spherical interpolation operations.

    Args:
        pPoint (pypoint): Point with geographic coordinates in degrees.
        use_high_precision (bool): Use extended precision for calculations (default: False).

    Returns:
        pynvector: An n-vector object with unit length.
    """

    if use_high_precision:
        dtype = _HIGH_PRECISION_DTYPE
        a = dtype(pPoint.dLatitude_radian)
        b = dtype(pPoint.dLongitude_radian)
    else:
        a = pPoint.dLatitude_radian
        b = pPoint.dLongitude_radian

    c = np.sin(a)
    e = np.cos(a)
    d = np.sin(b)
    f = np.cos(b)

    # Right-handed vector: x -> 0°E,0°N; y -> 90°E,0°N, z -> 90°N
    x = e * f
    y = e * d
    z = c

    point_dict = {"x": x, "y": y, "z": z}
    return pynvector(point_dict, use_high_precision=use_high_precision)


def nvector_to_point(pNvector: pynvector) -> pypoint:
    """
    Convert an n-vector object back to a pypoint object.

    Uses high precision (float128) for intermediate calculations if enabled,
    but returns float64 coordinates for pypoint storage.

    Args:
        pNvector (pynvector): N-vector on a unit sphere.

    Returns:
        pypoint: Point with latitude and longitude in degrees.
    """

    dtype = _HIGH_PRECISION_DTYPE if pNvector.use_high_precision else np.float64

    # Re-normalize in target precision to reduce drift before conversion
    x = dtype(pNvector.dX)
    y = dtype(pNvector.dY)
    z = dtype(pNvector.dZ)
    norm = np.sqrt(x * x + y * y + z * z)

    if norm > dtype(0.0):
        x = x / norm
        y = y / norm
        z = z / norm

    lat_rad = np.arctan2(z, np.sqrt(x * x + y * y))
    lon_rad = np.arctan2(y, x)

    return pypoint(
        {
            "dLongitude_degree": float(np.degrees(lon_rad)),
            "dLatitude_degree": float(np.degrees(lat_rad)),
        }
    )
