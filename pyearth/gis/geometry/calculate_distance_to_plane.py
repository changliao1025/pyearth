"""
Calculate distance from a point to a plane on Earth's surface.

This module provides functionality to calculate the perpendicular distance from a
geographic point to a plane defined by three other points, using 3D spherical geometry.
The plane passes through the Earth's center and the three defining points.
"""

import math
import numpy as np
from typing import Union

from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
)


def calculate_distance_to_plane(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> float:
    """
    Calculate the distance from a point to a plane defined by three points in 3D space.

    This function computes the perpendicular distance from a geographic point (point 2)
    to a plane that passes through the Earth's center and three defining points
    (points 1, 2, and 3). The calculation uses 3D Cartesian coordinates on a unit sphere.

    Args:
        dLongitude1_in: Longitude of first plane-defining point (degrees or radians)
        dLatitude1_in: Latitude of first plane-defining point (degrees or radians)
        dLongitude2_in: Longitude of query point (degrees or radians)
        dLatitude2_in: Latitude of query point (degrees or radians)
        dLongitude3_in: Longitude of third plane-defining point (degrees or radians)
        dLatitude3_in: Latitude of third plane-defining point (degrees or radians)
        iFlag_radian: If True, input coordinates are in radians; if False, in degrees (default: False)

    Returns:
        float: Perpendicular distance from point 2 to the plane defined by points 1, 2, and 3
               in unit sphere coordinates (dimensionless, typically < 2)

    Note:
        - Point 2 is the query point whose distance to the plane is calculated
        - Points 1 and 3 (along with point 2) define the plane
        - The plane passes through the Earth's center (origin)
        - Returns 0.0 if the three points are collinear (cannot define a unique plane)
        - Distance is in unit sphere coordinates, not meters
        - For great circle distance checking, see related functions

    Example:
        >>> # Check if a point lies on a great circle defined by two other points
        >>> distance = calculate_distance_to_plane(
        ...     -120.0, 37.0,  # Point 1: defines plane
        ...     -122.0, 38.0,  # Point 2: query point
        ...     -118.0, 36.0   # Point 3: defines plane
        ... )
        >>> if distance < 0.001:
        ...     print("Point 2 is approximately on the great circle through points 1 and 3")

    Algorithm:
        1. Convert all three geographic coordinates to 3D Cartesian (x,y,z) on unit sphere
        2. Create two vectors in the plane:
           - v1: from point 1 to point 2
           - v2: from point 1 to point 3
        3. Calculate plane normal vector using cross product: n = v1 × v2
        4. Check if normal is zero (points are collinear) → return 0.0
        5. Compute plane equation: Ax + By + Cz + D = 0 where:
           - (A, B, C) = normal vector
           - D = -(A*x1 + B*y1 + C*z1)
        6. Calculate perpendicular distance using point-to-plane formula:
           distance = |Ax2 + By2 + Cz2 + D| / sqrt(A² + B² + C²)

    References:
        - Point-to-plane distance: https://mathworld.wolfram.com/Point-PlaneDistance.html
        - Great circle checking: https://stackoverflow.com/questions/8204998/
    """
    # Convert to radians if input is in degrees
    if iFlag_radian:
        lon1_rad, lat1_rad = dLongitude1_in, dLatitude1_in
        lon2_rad, lat2_rad = dLongitude2_in, dLatitude2_in
        lon3_rad, lat3_rad = dLongitude3_in, dLatitude3_in
    else:
        lon1_rad, lat1_rad = np.radians([dLongitude1_in, dLatitude1_in])
        lon2_rad, lat2_rad = np.radians([dLongitude2_in, dLatitude2_in])
        lon3_rad, lat3_rad = np.radians([dLongitude3_in, dLatitude3_in])

    # Convert geographic coordinates to 3D Cartesian coordinates on unit sphere
    x1, y1, z1 = convert_longitude_latitude_to_sphere_3d(
        lon1_rad, lat1_rad
    )  # Plane point 1
    x2, y2, z2 = convert_longitude_latitude_to_sphere_3d(
        lon2_rad, lat2_rad
    )  # Query point
    x3, y3, z3 = convert_longitude_latitude_to_sphere_3d(
        lon3_rad, lat3_rad
    )  # Plane point 2

    # Calculate two vectors in the plane
    # v1: from point 1 to point 2
    v1_x, v1_y, v1_z = x2 - x1, y2 - y1, z2 - z1
    # v2: from point 1 to point 3
    v2_x, v2_y, v2_z = x3 - x1, y3 - y1, z3 - z1

    # Compute the normal vector using cross product: n = v1 × v2
    normal_x = v1_y * v2_z - v1_z * v2_y
    normal_y = v1_z * v2_x - v1_x * v2_z
    normal_z = v1_x * v2_y - v1_y * v2_x

    # Check if the normal vector is zero (points are collinear)
    # If collinear, the three points don't define a unique plane
    if abs(normal_x) < 1e-10 and abs(normal_y) < 1e-10 and abs(normal_z) < 1e-10:
        return 0.0

    # Calculate the plane equation coefficients: Ax + By + Cz + D = 0
    A, B, C = normal_x, normal_y, normal_z
    D = -(A * x1 + B * y1 + C * z1)

    # Calculate the perpendicular distance from point 2 to the plane
    # Formula: distance = |Ax + By + Cz + D| / sqrt(A² + B² + C²)
    numerator = abs(A * x2 + B * y2 + C * z2 + D)
    denominator = math.sqrt(A**2 + B**2 + C**2)
    distance = numerator / denominator

    return distance
