"""Bidirectional conversion between geographic and 3D Cartesian coordinates."""

from typing import Tuple, Union
import numpy as np

__all__ = [
    "convert_longitude_latitude_to_sphere_3d",
    "convert_sphere_3d_to_longitude_latitude",
]


def convert_longitude_latitude_to_sphere_3d(
    dLongitude_in: Union[float, np.ndarray],
    dLatitude_in: Union[float, np.ndarray],
    iFlag_radian: bool = False,
    radius: float = 1.0,
) -> np.ndarray:
    """Convert longitude/latitude to 3D Cartesian coordinates on a sphere.

    Transforms geographic coordinates (longitude, latitude) to 3D Cartesian
    coordinates (x, y, z) on a sphere with the specified radius (default unit sphere).

    Parameters
    ----------
    dLongitude_in : float or np.ndarray
        Longitude coordinate(s). In degrees by default, or radians if iFlag_radian=True.
        Valid range: [-180, 180] degrees or [-π, π] radians.
    dLatitude_in : float or np.ndarray
        Latitude coordinate(s). In degrees by default, or radians if iFlag_radian=True.
        Valid range: [-90, 90] degrees or [-π/2, π/2] radians.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians. If False (default), in degrees.
    radius : float, optional
        Sphere radius in desired distance units. Default is 1.0 (unit sphere).

    Returns
    -------
    np.ndarray
        3D Cartesian coordinates [x, y, z] on sphere. Shape (3,) for scalar input,
        or (N, 3) for 1D array input.

    Raises
    ------
    ValueError
        If latitude or longitude is outside valid range, or radius <= 0.
    TypeError
        If input coordinates cannot be converted to numeric arrays.
    """
    if radius <= 0:
        raise ValueError(f"Sphere radius must be positive. Got {radius}")

    try:
        lon_arr = np.asarray(dLongitude_in, dtype=float)
        lat_arr = np.asarray(dLatitude_in, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(f"Coordinates must be numeric values. Error: {e}")

    if not iFlag_radian:
        if np.any((lat_arr < -90.0) | (lat_arr > 90.0)):
            raise ValueError("Latitude must be in range [-90, 90] degrees.")
        if np.any((lon_arr < -180.0) | (lon_arr > 180.0)):
            raise ValueError("Longitude must be in range [-180, 180] degrees.")
        lon_rad = np.radians(lon_arr)
        lat_rad = np.radians(lat_arr)
    else:
        if np.any((lat_arr < -np.pi / 2.0) | (lat_arr > np.pi / 2.0)):
            raise ValueError("Latitude must be in range [-π/2, π/2] radians.")
        if np.any((lon_arr < -np.pi) | (lon_arr > np.pi)):
            raise ValueError("Longitude must be in range [-π, π] radians.")
        lon_rad = lon_arr
        lat_rad = lat_arr

    cos_lat = np.cos(lat_rad)
    x = radius * cos_lat * np.cos(lon_rad)
    y = radius * cos_lat * np.sin(lon_rad)
    z = radius * np.sin(lat_rad)

    if lon_arr.ndim == 0:
        return np.array([float(x), float(y), float(z)])
    return np.column_stack([x, y, z])


def convert_sphere_3d_to_longitude_latitude(
    x: Union[float, np.ndarray],
    y: Union[float, np.ndarray],
    z: Union[float, np.ndarray],
) -> Tuple[float, float]:
    """Convert 3D Cartesian coordinates to longitude/latitude coordinates.

    Transforms 3D Cartesian coordinates (x, y, z) to geographic coordinates
    (longitude, latitude) on a sphere. Works with any sphere radius.

    Parameters
    ----------
    x : float or np.ndarray
        X coordinate in 3D Cartesian space.
    y : float or np.ndarray
        Y coordinate in 3D Cartesian space.
    z : float or np.ndarray
        Z coordinate in 3D Cartesian space.

    Returns
    -------
    Tuple[float, float]
        Geographic coordinates (longitude, latitude) in decimal degrees:
        - longitude: Range [-180, 180] degrees
        - latitude: Range [-90, 90] degrees

    Raises
    ------
    ValueError
        If all coordinates are zero or contain NaN/inf.
    TypeError
        If inputs are not numeric.
    """
    try:
        x_val = np.asarray(x, dtype=float)
        y_val = np.asarray(y, dtype=float)
        z_val = np.asarray(z, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(f"Coordinates must be numeric values. Error: {e}")

    if np.any(np.isnan([x_val, y_val, z_val])):
        raise ValueError("Coordinates cannot contain NaN values.")
    if np.any(np.isinf([x_val, y_val, z_val])):
        raise ValueError("Coordinates cannot contain infinite values.")

    norm = np.sqrt(x_val**2 + y_val**2 + z_val**2)
    if np.all(norm == 0) or np.any(np.isclose(norm, 0, atol=1e-15)):
        raise ValueError("Cannot convert origin (0, 0, 0) to longitude/latitude.")

    x_norm = x_val / norm
    y_norm = y_val / norm
    z_norm = z_val / norm

    lon_rad = np.arctan2(y_norm, x_norm)
    z_clamped = np.clip(z_norm, -1.0, 1.0)
    lat_rad = np.arcsin(z_clamped)

    lon_deg = np.degrees(lon_rad)
    lat_deg = np.degrees(lat_rad)

    if x_val.ndim == 0:
        return float(lon_deg), float(lat_deg)
    return lon_deg, lat_deg
