import numpy as np
from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
)
from pyearth.gis.geometry.calculate_angle_between_vectors_degrees import (
    calculate_angle_between_vectors_degrees,
)


def calculate_angle_between_point(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> float:
    """

    Calculates the angle between three points on a sphere.

    Args:
        dLongitude1_in (float): Longitude of the first point.
        dLatitude1_in (float): Latitude of the first point.
        dLongitude2_in (float): Longitude of the second point (the middle one).
        dLatitude2_in (float): Latitude of the second point (the middle one).
        dLongitude3_in (float): Longitude of the third point.
        dLatitude3_in (float): Latitude of the third point.
        iFlag_radian (bool, optional): If True, input coordinates are in radians. Defaults to False (degrees).

    Returns:
        float: The angle in degrees between the vectors from the middle point to the other two.
    """

    # Determine whether the downstream conversion routine should treat inputs as radians.
    iFlag_convert_radian = True if iFlag_radian else None

    # The points in 3D space
    a3 = convert_longitude_latitude_to_sphere_3d(
        dLongitude1_in, dLatitude1_in, iFlag_convert_radian
    )
    b3 = convert_longitude_latitude_to_sphere_3d(
        dLongitude2_in, dLatitude2_in, iFlag_convert_radian
    )
    c3 = convert_longitude_latitude_to_sphere_3d(
        dLongitude3_in, dLatitude3_in, iFlag_convert_radian
    )

    # Vectors in 3D space
    a3vec = a3 - b3
    c3vec = c3 - b3

    angle3deg = calculate_angle_between_vectors_degrees(a3vec, c3vec)
    return angle3deg
