"""
Convert geographic coordinates to 3D Cartesian coordinates on a unit sphere.

This module provides utilities for converting longitude/latitude coordinates
to 3D Cartesian (x, y, z) coordinates on a sphere with radius 1.
"""
import numpy as np
from typing import Union


def convert_longitude_latitude_to_sphere_3d(
    dLongitude_in: float,
    dLatitude_in: float,
    iFlag_radian: bool = False
) -> np.ndarray:
    """Convert longitude/latitude to 3D Cartesian coordinates on a unit sphere.

    Transforms geographic coordinates (longitude, latitude) to 3D Cartesian
    coordinates (x, y, z) on a sphere with radius 1. This is useful for
    spherical geometry calculations, great circle distances, and 3D visualizations.

    Parameters
    ----------
    dLongitude_in : float
        Longitude coordinate. In degrees by default, or radians if iFlag_radian=True.
        Valid range: [-180, 180] degrees or [-π, π] radians.
    dLatitude_in : float
        Latitude coordinate. In degrees by default, or radians if iFlag_radian=True.
        Valid range: [-90, 90] degrees or [-π/2, π/2] radians.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians. If False (default),
        input coordinates are in degrees.

    Returns
    -------
    np.ndarray
        3D Cartesian coordinates [x, y, z] on a unit sphere (radius = 1).
        - x = cos(lat) * cos(lon)
        - y = cos(lat) * sin(lon)
        - z = sin(lat)

    Raises
    ------
    ValueError
        If latitude is outside valid range [-90, 90] degrees or [-π/2, π/2] radians.
        If longitude is outside valid range [-180, 180] degrees or [-π, π] radians.
    TypeError
        If input coordinates cannot be converted to float.

    Notes
    -----
    - The returned vector has magnitude (norm) equal to 1
    - This uses the standard geographic coordinate system:
      * Longitude: 0° at Prime Meridian, positive East, negative West
      * Latitude: 0° at Equator, positive North, negative South
    - The 3D coordinate system:
      * x-axis: points to (0°N, 0°E) - Equator at Prime Meridian
      * y-axis: points to (0°N, 90°E) - Equator at 90° East
      * z-axis: points to North Pole (90°N)

    Examples
    --------
    >>> # North Pole
    >>> coords = convert_longitude_latitude_to_sphere_3d(0, 90)
    >>> np.allclose(coords, [0, 0, 1])
    True

    >>> # Equator at Prime Meridian
    >>> coords = convert_longitude_latitude_to_sphere_3d(0, 0)
    >>> np.allclose(coords, [1, 0, 0])
    True

    >>> # Equator at 90° East
    >>> coords = convert_longitude_latitude_to_sphere_3d(90, 0)
    >>> np.allclose(coords, [0, 1, 0])
    True

    >>> # Using radians
    >>> coords = convert_longitude_latitude_to_sphere_3d(np.pi/2, 0, iFlag_radian=True)
    >>> np.allclose(coords, [0, 1, 0])
    True

    >>> # Verify unit sphere (magnitude = 1)
    >>> coords = convert_longitude_latitude_to_sphere_3d(45, 45)
    >>> np.allclose(np.linalg.norm(coords), 1.0)
    True
    """
    # Validate and convert inputs to float
    try:
        longitude = float(dLongitude_in)
        latitude = float(dLatitude_in)
    except (TypeError, ValueError) as e:
        raise TypeError(
            f"Coordinates must be numeric values. "
            f"Got longitude={dLongitude_in}, latitude={dLatitude_in}. Error: {e}"
        )

    # Convert to radians if needed and validate ranges
    if not iFlag_radian:
        # Input in degrees - validate and convert
        if not -90 <= latitude <= 90:
            raise ValueError(
                f"Latitude must be in range [-90, 90] degrees. Got {latitude}°"
            )
        if not -180 <= longitude <= 180:
            raise ValueError(
                f"Longitude must be in range [-180, 180] degrees. Got {longitude}°"
            )
        longitude_rad, latitude_rad = np.radians([longitude, latitude])
    else:
        # Input in radians - validate
        if not -np.pi/2 <= latitude <= np.pi/2:
            raise ValueError(
                f"Latitude must be in range [-π/2, π/2] radians. Got {latitude} rad"
            )
        if not -np.pi <= longitude <= np.pi:
            raise ValueError(
                f"Longitude must be in range [-π, π] radians. Got {longitude} rad"
            )
        longitude_rad = longitude
        latitude_rad = latitude

    # Convert to 3D Cartesian coordinates on unit sphere
    # x = cos(lat) * cos(lon)
    # y = cos(lat) * sin(lon)
    # z = sin(lat)
    cos_lat = np.cos(latitude_rad)

    x = cos_lat * np.cos(longitude_rad)
    y = cos_lat * np.sin(longitude_rad)
    z = np.sin(latitude_rad)

    return np.array([x, y, z])
