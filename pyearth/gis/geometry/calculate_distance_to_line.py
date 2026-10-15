"""Calculate distance from a point to a line segment on Earth's surface."""

from typing import Optional
import numpy as np

from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
)
from pyearth.system.define_global_variables import earth_radius

__all__ = ["calculate_distance_to_line"]


def calculate_distance_to_line(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
    radius: Optional[float] = None,
) -> float:
    """Calculate the shortest distance from a point to a line segment on Earth's surface.

    Point 2 is the query point. Points 1 and 3 define the segment endpoints.

    Parameters
    ----------
    dLongitude1_in : float
        Longitude of line segment start point.
    dLatitude1_in : float
        Latitude of line segment start point.
    dLongitude2_in : float
        Longitude of query point.
    dLatitude2_in : float
        Latitude of query point.
    dLongitude3_in : float
        Longitude of line segment end point.
    dLatitude3_in : float
        Latitude of line segment end point.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians; if False (default), in degrees.
    radius : float, optional
        Sphere radius in meters. Defaults to earth_radius.

    Returns
    -------
    float
        Shortest distance in meters.
    """
    if iFlag_radian:
        lon1_rad, lat1_rad = dLongitude1_in, dLatitude1_in
        lon2_rad, lat2_rad = dLongitude2_in, dLatitude2_in
        lon3_rad, lat3_rad = dLongitude3_in, dLatitude3_in
    else:
        lon1_rad, lat1_rad = np.radians([dLongitude1_in, dLatitude1_in])
        lon2_rad, lat2_rad = np.radians([dLongitude2_in, dLatitude2_in])
        lon3_rad, lat3_rad = np.radians([dLongitude3_in, dLatitude3_in])

    x1, y1, z1 = convert_longitude_latitude_to_sphere_3d(
        lon1_rad, lat1_rad, iFlag_radian=True
    )
    x2, y2, z2 = convert_longitude_latitude_to_sphere_3d(
        lon2_rad, lat2_rad, iFlag_radian=True
    )
    x3, y3, z3 = convert_longitude_latitude_to_sphere_3d(
        lon3_rad, lat3_rad, iFlag_radian=True
    )

    p1 = np.array([x1, y1, z1])
    p2 = np.array([x2, y2, z2])
    p3 = np.array([x3, y3, z3])

    line_vector = p3 - p1
    point_vector = p2 - p1
    line_length = np.linalg.norm(line_vector)

    if line_length < 1e-10:
        distance_on_sphere = np.linalg.norm(point_vector)
    else:
        line_unit_vector = line_vector / line_length
        projection_length = np.dot(point_vector, line_unit_vector)
        projection_length_clamped = np.clip(projection_length, 0.0, line_length)
        projection_point = p1 + projection_length_clamped * line_unit_vector
        distance_on_sphere = np.linalg.norm(p2 - projection_point)

    r = radius if radius is not None else earth_radius
    return float(distance_on_sphere * r)
