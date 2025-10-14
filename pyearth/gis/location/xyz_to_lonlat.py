"""
Convert 3D Cartesian coordinates to longitude/latitude on a sphere.

This module provides utilities for converting 3D Cartesian (x, y, z) coordinates
to geographic coordinates (longitude, latitude) on a unit sphere or any spherical surface.
"""
from typing import Tuple, Union
import numpy as np


def xyz_to_lonlat(
    x: Union[float, np.ndarray],
    y: Union[float, np.ndarray],
    z: Union[float, np.ndarray]
) -> Tuple[float, float]:
    """Convert 3D Cartesian coordinates to longitude/latitude coordinates.

    Transforms 3D Cartesian coordinates (x, y, z) to geographic coordinates
    (longitude, latitude) on a sphere. The input coordinates are automatically
    normalized to a unit sphere before conversion.

    Parameters
    ----------
    x : float or np.ndarray
        X coordinate in 3D Cartesian space.
        - Points to (0°N, 0°E) - Equator at Prime Meridian
        - For unit sphere: x = cos(lat) * cos(lon)
    y : float or np.ndarray
        Y coordinate in 3D Cartesian space.
        - Points to (0°N, 90°E) - Equator at 90° East
        - For unit sphere: y = cos(lat) * sin(lon)
    z : float or np.ndarray
        Z coordinate in 3D Cartesian space.
        - Points to North Pole (90°N)
        - For unit sphere: z = sin(lat)

    Returns
    -------
    Tuple[float, float]
        Geographic coordinates (longitude, latitude) in decimal degrees:
        - longitude: Range [-180, 180] degrees (negative = West, positive = East)
        - latitude: Range [-90, 90] degrees (negative = South, positive = North)

    Raises
    ------
    ValueError
        If all coordinates are zero (degenerate case, no unique point).
        If inputs contain NaN or infinite values.
    TypeError
        If inputs cannot be converted to numeric values.

    Notes
    -----
    - Input coordinates are normalized before conversion (divided by their magnitude)
    - Works for any non-zero radius; automatically projects to unit sphere
    - Uses atan2 for longitude to handle all quadrants correctly
    - Uses asin for latitude (assumes normalized coordinates)
    - This is the inverse of convert_longitude_latitude_to_sphere_3d

    Coordinate System
    -----------------
    The 3D coordinate system follows standard geographic conventions:
    - X-axis: Points to (0°N, 0°E) - Equator at Prime Meridian
    - Y-axis: Points to (0°N, 90°E) - Equator at 90° East
    - Z-axis: Points to North Pole (90°N)
    - Origin: Earth's center

    Mathematical Formulas
    ---------------------
    After normalization to unit sphere:
        longitude = atan2(y, x)
        latitude = asin(z)

    Examples
    --------
    >>> # North Pole
    >>> lon, lat = xyz_to_lonlat(0, 0, 1)
    >>> np.allclose([lon, lat], [0, 90])
    True

    >>> # Equator at Prime Meridian
    >>> lon, lat = xyz_to_lonlat(1, 0, 0)
    >>> np.allclose([lon, lat], [0, 0])
    True

    >>> # Equator at 90° East
    >>> lon, lat = xyz_to_lonlat(0, 1, 0)
    >>> np.allclose([lon, lat], [90, 0])
    True

    >>> # Equator at 180° (or -180°)
    >>> lon, lat = xyz_to_lonlat(-1, 0, 0)
    >>> np.allclose(abs(lon), 180) and np.allclose(lat, 0)
    True

    >>> # Point at 45°N, 45°E (on unit sphere)
    >>> x = np.cos(np.radians(45)) * np.cos(np.radians(45))
    >>> y = np.cos(np.radians(45)) * np.sin(np.radians(45))
    >>> z = np.sin(np.radians(45))
    >>> lon, lat = xyz_to_lonlat(x, y, z)
    >>> np.allclose([lon, lat], [45, 45], atol=1e-10)
    True

    >>> # Works with any radius (auto-normalized)
    >>> lon, lat = xyz_to_lonlat(2, 0, 0)  # 2x unit sphere radius
    >>> np.allclose([lon, lat], [0, 0])
    True

    See Also
    --------
    convert_longitude_latitude_to_sphere_3d : Inverse transformation (lon/lat to xyz)
    numpy.arctan2 : 2-argument arctangent for correct quadrant
    numpy.arcsin : Arcsine function

    References
    ----------
    - Spherical coordinate system: https://en.wikipedia.org/wiki/Spherical_coordinate_system
    - Geographic coordinate system: https://en.wikipedia.org/wiki/Geographic_coordinate_system
    """
    # Validate and convert inputs
    try:
        x_val = float(x) if np.isscalar(x) else np.asarray(x, dtype=float)
        y_val = float(y) if np.isscalar(y) else np.asarray(y, dtype=float)
        z_val = float(z) if np.isscalar(z) else np.asarray(z, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(
            f"Coordinates must be numeric values. "
            f"Got x={x}, y={y}, z={z}. Error: {e}"
        )

    # Check for NaN or infinite values
    if np.any(np.isnan([x_val, y_val, z_val])):
        raise ValueError(
            "Coordinates cannot contain NaN values. "
            f"Got x={x_val}, y={y_val}, z={z_val}"
        )
    if np.any(np.isinf([x_val, y_val, z_val])):
        raise ValueError(
            "Coordinates cannot contain infinite values. "
            f"Got x={x_val}, y={y_val}, z={z_val}"
        )

    # Calculate magnitude (norm)
    norm = np.sqrt(x_val**2 + y_val**2 + z_val**2)

    # Check for degenerate case (all zeros)
    if norm == 0 or np.isclose(norm, 0, atol=1e-15):
        raise ValueError(
            "Cannot convert origin (0, 0, 0) to longitude/latitude. "
            "All coordinates are zero (or extremely close to zero)."
        )

    # Normalize coordinates to unit sphere
    x_norm = x_val / norm
    y_norm = y_val / norm
    z_norm = z_val / norm

    # Convert to spherical coordinates
    # Longitude: use atan2 for correct quadrant handling
    lon_rad = np.arctan2(y_norm, x_norm)

    # Latitude: use asin (z is already normalized)
    # Clamp z_norm to [-1, 1] to handle floating-point errors
    z_clamped = np.clip(z_norm, -1.0, 1.0)
    lat_rad = np.arcsin(z_clamped)

    # Convert radians to degrees
    lon_deg = np.degrees(lon_rad)
    lat_deg = np.degrees(lat_rad)

    return float(lon_deg), float(lat_deg)
