"""Calculate perpendicular distance from a point to a great-circle plane."""

import math
from typing import Optional

import numpy as np

from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
)
from pyearth.system.define_global_variables import earth_radius

__all__ = ["calculate_distance_to_plane"]


def calculate_distance_to_plane(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
    radius: Optional[float] = None,
) -> float:
    """Calculate perpendicular distance from query point 2 to the great-circle plane of points 1 and 3.

    The great-circle plane passes through Earth's center (origin) and points 1 and 3.
    Point 2 is the query point.

    Parameters
    ----------
    dLongitude1_in : float
        Longitude of first plane-defining point.
    dLatitude1_in : float
        Latitude of first plane-defining point.
    dLongitude2_in : float
        Longitude of query point.
    dLatitude2_in : float
        Latitude of query point.
    dLongitude3_in : float
        Longitude of second plane-defining point.
    dLatitude3_in : float
        Latitude of second plane-defining point.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians; if False (default), in degrees.
    radius : float, optional
        If provided, returns distance in meters (or units of radius).
        If None (default), returns dimensionless distance on unit sphere.

    Returns
    -------
    float
        Perpendicular distance from point 2 to the plane through (0,0,0), point 1, and point 3.
    """
    if iFlag_radian:
        lon1_rad, lat1_rad = dLongitude1_in, dLatitude1_in
        lon2_rad, lat2_rad = dLongitude2_in, dLatitude2_in
        lon3_rad, lat3_rad = dLongitude3_in, dLatitude3_in
    else:
        lon1_rad, lat1_rad = np.radians([dLongitude1_in, dLatitude1_in])
        lon2_rad, lat2_rad = np.radians([dLongitude2_in, dLatitude2_in])
        lon3_rad, lat3_rad = np.radians([dLongitude3_in, dLatitude3_in])

    p1 = convert_longitude_latitude_to_sphere_3d(lon1_rad, lat1_rad, iFlag_radian=True)
    p2 = convert_longitude_latitude_to_sphere_3d(lon2_rad, lat2_rad, iFlag_radian=True)
    p3 = convert_longitude_latitude_to_sphere_3d(lon3_rad, lat3_rad, iFlag_radian=True)

    # Great-circle plane through origin, p1, and p3 has normal n = p1 x p3
    normal = np.cross(p1, p3)
    norm_mag = np.linalg.norm(normal)

    if norm_mag < 1e-12:
        # Collinear or antipodal points 1 and 3 do not define a unique great-circle plane
        return 0.0

    # Distance to plane through origin: |n . p2| / ||n||
    dist_unit = abs(float(np.dot(normal, p2))) / float(norm_mag)

    if radius is not None:
        return float(dist_unit * radius)
    return float(dist_unit)
