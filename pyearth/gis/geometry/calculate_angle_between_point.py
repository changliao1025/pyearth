"""Calculate angle between points on a sphere."""

import numpy as np

from pyearth.gis.geometry.calculate_angle_between_vectors_degrees import (
    calculate_angle_between_vectors_degrees,
)
from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
)

__all__ = [
    "calculate_angle_between_point",
]


def calculate_angle_between_point(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> float:
    """Calculate the angle between three points on a sphere.

    Computes the angle at point 2 (middle vertex) formed by the vectors
    to point 1 and point 3.

    Parameters
    ----------
    dLongitude1_in : float
        Longitude of the first point.
    dLatitude1_in : float
        Latitude of the first point.
    dLongitude2_in : float
        Longitude of the second point (the vertex/middle point).
    dLatitude2_in : float
        Latitude of the second point (the vertex/middle point).
    dLongitude3_in : float
        Longitude of the third point.
    dLatitude3_in : float
        Latitude of the third point.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians. Defaults to False (degrees).

    Returns
    -------
    float
        The angle in degrees [0, 180] between the vectors from the middle point.
    """
    a3 = convert_longitude_latitude_to_sphere_3d(
        dLongitude1_in, dLatitude1_in, iFlag_radian=iFlag_radian
    )
    b3 = convert_longitude_latitude_to_sphere_3d(
        dLongitude2_in, dLatitude2_in, iFlag_radian=iFlag_radian
    )
    c3 = convert_longitude_latitude_to_sphere_3d(
        dLongitude3_in, dLatitude3_in, iFlag_radian=iFlag_radian
    )

    a3vec = a3 - b3
    c3vec = c3 - b3

    return calculate_angle_between_vectors_degrees(a3vec, c3vec)
