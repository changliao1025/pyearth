"""
Calculate angle between point normal for geographic coordinates.

This module provides functionality to calculate the angle at a point (middle point)
formed by three points on a sphere, useful for geographic calculations involving
lines and polygons on Earth's surface.
"""

import numpy as np
from typing import Union

from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
)


def calculate_angle_between_point_normal(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> float:
    """
    Calculate the angle at point (point 2) formed by three points on a sphere.

    This function computes the angle at the middle point of three geographic points
    by converting them to 3D Cartesian coordinates on a unit sphere, then calculating
    the angle using vector operations.

    Args:
        dLongitude1_in: Longitude of first point (degrees or radians)
        dLatitude1_in: Latitude of first point (degrees or radians)
        dLongitude2_in: Longitude of middle point (degrees or radians)
        dLatitude2_in: Latitude of middle point (degrees or radians)
        dLongitude3_in: Longitude of third point (degrees or radians)
        dLatitude3_in: Latitude of third point (degrees or radians)
        iFlag_radian: If True, input coordinates are in radians; if False, in degrees (default: False)

    Returns:
    float: Angle at point (point 2) in degrees, ranging from 0 to 360

    Note:
    - The angle is measured from point 1 to point 3 around point 2 (the point)
    - Returns values in range [0, 360) degrees
    - Point 2 is the point where the angle is measured
    - Uses 3D vector cross product and dot product for accurate spherical calculations

    Example:
    >>> # Calculate angle at point for three geographic points
    >>> angle = calculate_angle_between_point_normal(
    ...     -122.0, 37.0,  # Point 1: San Francisco area
    ...     -121.0, 37.5,  # Point 2: Middle point
    ...     -120.0, 38.0   # Point 3: Sacramento area
    ... )
    >>> print(f"Angle at point: {angle:.2f} degrees")

    Algorithm:
    1. Convert lat/lon to radians if needed
    2. Project points onto 3D unit sphere
    3. Create vectors from point to other two points
    4. Calculate angle using arctan2(det, dot) where:
           - dot = dot product of the two vectors
           - det = determinant using cross product
        5. Convert to degrees and normalize to [0, 360)
    """
    # Convert to radians if input is in degrees
    if iFlag_radian:
        dLongitude1_radian = dLongitude1_in
        dLatitude1_radian = dLatitude1_in
        dLongitude2_radian = dLongitude2_in
        dLatitude2_radian = dLatitude2_in
        dLongitude3_radian = dLongitude3_in
        dLatitude3_radian = dLatitude3_in
    else:
        dLongitude1_radian, dLatitude1_radian = np.radians(
            [dLongitude1_in, dLatitude1_in]
        )
        dLongitude2_radian, dLatitude2_radian = np.radians(
            [dLongitude2_in, dLatitude2_in]
        )
        dLongitude3_radian, dLatitude3_radian = np.radians(
            [dLongitude3_in, dLatitude3_in]
        )

    # Convert geographic coordinates to 3D Cartesian coordinates on unit sphere
    point1_3d = convert_longitude_latitude_to_sphere_3d(
        dLongitude1_radian, dLatitude1_radian
    )
    point2_3d = convert_longitude_latitude_to_sphere_3d(
        dLongitude2_radian, dLatitude2_radian
    )  # Middle point
    point3_3d = convert_longitude_latitude_to_sphere_3d(
        dLongitude3_radian, dLatitude3_radian
    )

    # Create vectors from point (point 2) to the other two points
    vector_to_point1 = point1_3d - point2_3d
    vector_to_point3 = point3_3d - point2_3d

    # Calculate angle using dot product and cross product
    dot_product = np.dot(vector_to_point1, vector_to_point3)
    cross_product = np.cross(vector_to_point1, vector_to_point3)
    determinant = np.dot(point2_3d, cross_product)

    # Calculate angle in radians using atan2 for proper quadrant handling
    angle_radians = np.arctan2(determinant, dot_product)

    # Convert to degrees
    angle_degrees = np.degrees(angle_radians)

    # Normalize to [0, 360) range
    if angle_degrees < 0:
        angle_degrees += 360.0

    return angle_degrees
