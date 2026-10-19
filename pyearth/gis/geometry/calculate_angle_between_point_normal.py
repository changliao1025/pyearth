"""Calculate angle between point normal for geographic coordinates."""

import numpy as np

from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
)

__all__ = ["calculate_angle_between_point_normal"]


def calculate_angle_between_point_normal(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> float:
    """Calculate the angle at point 2 formed by three points on a sphere.

    Computes the oriented angle measured at the middle point (point 2) from
    point 1 to point 3 around the local surface normal, normalized to [0, 360) degrees.

    Parameters
    ----------
    dLongitude1_in : float
        Longitude of first point (degrees or radians).
    dLatitude1_in : float
        Latitude of first point (degrees or radians).
    dLongitude2_in : float
        Longitude of middle point (degrees or radians).
    dLatitude2_in : float
        Latitude of middle point (degrees or radians).
    dLongitude3_in : float
        Longitude of third point (degrees or radians).
    dLatitude3_in : float
        Latitude of third point (degrees or radians).
    iFlag_radian : bool, optional
        If True, input coordinates are in radians. Default is False.

    Returns
    -------
    float
        Angle at point 2 in degrees, ranging from 0 to 360.
    """
    if iFlag_radian:
        lon1_rad, lat1_rad = dLongitude1_in, dLatitude1_in
        lon2_rad, lat2_rad = dLongitude2_in, dLatitude2_in
        lon3_rad, lat3_rad = dLongitude3_in, dLatitude3_in
    else:
        lon1_rad, lat1_rad = np.radians([dLongitude1_in, dLatitude1_in])
        lon2_rad, lat2_rad = np.radians([dLongitude2_in, dLatitude2_in])
        lon3_rad, lat3_rad = np.radians([dLongitude3_in, dLatitude3_in])

    point1_3d = convert_longitude_latitude_to_sphere_3d(
        lon1_rad, lat1_rad, iFlag_radian=True
    )
    point2_3d = convert_longitude_latitude_to_sphere_3d(
        lon2_rad, lat2_rad, iFlag_radian=True
    )
    point3_3d = convert_longitude_latitude_to_sphere_3d(
        lon3_rad, lat3_rad, iFlag_radian=True
    )

    vector_to_point1 = point1_3d - point2_3d
    vector_to_point3 = point3_3d - point2_3d

    dot_product = np.dot(vector_to_point1, vector_to_point3)
    cross_product = np.cross(vector_to_point1, vector_to_point3)
    determinant = np.dot(point2_3d, cross_product)

    angle_radians = np.arctan2(determinant, dot_product)
    angle_degrees = float(np.degrees(angle_radians))

    if angle_degrees < 0.0:
        angle_degrees += 360.0

    return angle_degrees
