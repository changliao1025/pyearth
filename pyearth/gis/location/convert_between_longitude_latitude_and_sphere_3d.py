"""
Bidirectional conversion between geographic and 3D Cartesian coordinates.

This module provides utilities for converting between longitude/latitude
coordinates and 3D Cartesian (x, y, z) coordinates on a sphere. These
conversions are fundamental for spherical geometry calculations, great
circle computations, and 3D visualizations of geographic data.

Main Functions
--------------
convert_longitude_latitude_to_sphere_3d : Convert lon/lat to 3D Cartesian (x, y, z)
convert_sphere_3d_to_longitude_latitude : Convert 3D Cartesian to lon/lat (inverse)

Coordinate Systems
------------------
Geographic (Longitude/Latitude):
    - Longitude: [-180, 180] degrees (0° = Prime Meridian, + East, - West)
    - Latitude: [-90, 90] degrees (0° = Equator, + North, - South)

3D Cartesian (on unit sphere):
    - X-axis: Points to (0°N, 0°E) - Equator at Prime Meridian
    - Y-axis: Points to (0°N, 90°E) - Equator at 90° East
    - Z-axis: Points to North Pole (90°N)
    - Origin: Center of sphere (Earth's center)

Mathematical Relationships
--------------------------
Forward (lon/lat → xyz):
    x = cos(lat) * cos(lon)
    y = cos(lat) * sin(lon)
    z = sin(lat)

Inverse (xyz → lon/lat):
    lon = atan2(y, x)
    lat = asin(z / ||xyz||)

Use Cases
---------
- Great circle distance and bearing calculations
- Spherical interpolation (slerp)
- 3D visualization of geographic data
- Vector operations on spherical surfaces
- Mesh generation for spherical domains

See Also
--------
calculate_distance_based_on_longitude_latitude : Great circle distance
calculate_intersect_on_great_circle : Great circle intersection calculations

References
----------
.. [1] Spherical coordinate system:
       https://en.wikipedia.org/wiki/Spherical_coordinate_system
.. [2] Geographic coordinate system:
       https://en.wikipedia.org/wiki/Geographic_coordinate_system
"""

import numpy as np
from typing import Union, Tuple


def convert_longitude_latitude_to_sphere_3d(
    dLongitude_in: float, dLatitude_in: float, iFlag_radian: bool = False
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
        if not -np.pi / 2 <= latitude <= np.pi / 2:
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


def convert_sphere_3d_to_longitude_latitude(
    x: Union[float, np.ndarray],
    y: Union[float, np.ndarray],
    z: Union[float, np.ndarray],
) -> Tuple[float, float]:
    """Convert 3D Cartesian coordinates to longitude/latitude coordinates.

    Transforms 3D Cartesian coordinates (x, y, z) to geographic coordinates
    (longitude, latitude) on a sphere. The input coordinates are automatically
    normalized to a unit sphere before conversion.

    This is the inverse operation of convert_longitude_latitude_to_sphere_3d.

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
    - Z-coordinate is clamped to [-1, 1] to handle floating-point precision errors

    Mathematical Formulas
    ---------------------
    After normalization to unit sphere:
        longitude = atan2(y, x)
        latitude = asin(z)

    Examples
    --------
    >>> # North Pole
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(0, 0, 1)
    >>> np.allclose([lon, lat], [0, 90])
    True

    >>> # Equator at Prime Meridian
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(1, 0, 0)
    >>> np.allclose([lon, lat], [0, 0])
    True

    >>> # Equator at 90° East
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(0, 1, 0)
    >>> np.allclose([lon, lat], [90, 0])
    True

    >>> # Equator at 180° (or -180°)
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(-1, 0, 0)
    >>> np.allclose(abs(lon), 180) and np.allclose(lat, 0)
    True

    >>> # Point at 45°N, 45°E (on unit sphere)
    >>> x = np.cos(np.radians(45)) * np.cos(np.radians(45))
    >>> y = np.cos(np.radians(45)) * np.sin(np.radians(45))
    >>> z = np.sin(np.radians(45))
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(x, y, z)
    >>> np.allclose([lon, lat], [45, 45], atol=1e-10)
    True

    >>> # Works with any radius (auto-normalized)
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(2, 0, 0)
    >>> np.allclose([lon, lat], [0, 0])
    True

    >>> # Verify round-trip conversion
    >>> xyz = convert_longitude_latitude_to_sphere_3d(45, 30)
    >>> lon, lat = convert_sphere_3d_to_longitude_latitude(*xyz)
    >>> np.allclose([lon, lat], [45, 30])
    True

    See Also
    --------
    convert_longitude_latitude_to_sphere_3d : Inverse transformation (lon/lat to xyz)
    numpy.arctan2 : 2-argument arctangent for correct quadrant
    numpy.arcsin : Arcsine function

    References
    ----------
    .. [1] Spherical coordinate system:
           https://en.wikipedia.org/wiki/Spherical_coordinate_system
    .. [2] Geographic coordinate system:
           https://en.wikipedia.org/wiki/Geographic_coordinate_system
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
