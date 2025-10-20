

cimport cython
from libc.math cimport M_PI, sin, cos
import numpy as np

"""
Geospatial 3D conversion utilities for PyEarth.

Authors: Chang Liao
"""

"""
"""
# Authors: Chang Liao

#constant

cdef double dRadius = 6378137.0

@cython.boundscheck(False)
cpdef np.ndarray convert_longitude_latitude_to_sphere_3d(double dLongitude_in, double dLatitude_in, bint bFlag_radian=False):
    """
    Convert geographic coordinates (longitude, latitude) to 3D Cartesian coordinates on a unit sphere.

    Parameters
    ----------
    dLongitude_in : double
        Longitude coordinate. In degrees by default, or radians if bFlag_radian=True.
        Valid range: [-180, 180] degrees or [-π, π] radians.
    dLatitude_in : double
        Latitude coordinate. In degrees by default, or radians if bFlag_radian=True.
        Valid range: [-90, 90] degrees or [-π/2, π/2] radians.
    bFlag_radian : bint, optional
        If True, input coordinates are in radians. If False (default), input coordinates are in degrees.

    Returns
    -------
    np.ndarray
        3D Cartesian coordinates [x, y, z] on a unit sphere (radius = 1).

    Raises
    ------
    ValueError
        If latitude or longitude is outside valid range for the selected unit.

    Notes
    -----
    - The returned vector has magnitude (norm) equal to 1
    - Longitude: 0° at Prime Meridian, positive East, negative West
    - Latitude: 0° at Equator, positive North, negative South
    """
    cdef double longitude, latitude
    cdef double longitude_rad, latitude_rad
    cdef double x, y, z

    longitude = dLongitude_in
    latitude = dLatitude_in

    if not bFlag_radian:
        if latitude < -90.0 or latitude > 90.0:
            raise ValueError("Latitude must be in range [-90, 90] degrees.")
        if longitude < -180.0 or longitude > 180.0:
            raise ValueError("Longitude must be in range [-180, 180] degrees.")
        longitude_rad = longitude / 180.0 * M_PI
        latitude_rad = latitude / 180.0 * M_PI
    else:
        if latitude < -M_PI/2 or latitude > M_PI/2:
            raise ValueError("Latitude must be in range [-π/2, π/2] radians.")
        if longitude < -M_PI or longitude > M_PI:
            raise ValueError("Longitude must be in range [-π, π] radians.")
        longitude_rad = longitude
        latitude_rad = latitude

    x = cos(latitude_rad) * cos(longitude_rad)
    y = cos(latitude_rad) * sin(longitude_rad)
    z = sin(latitude_rad)

    return np.array([x, y, z], dtype=np.float64)

