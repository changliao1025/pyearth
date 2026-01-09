"""
Calculate distance from a point to a line segment on Earth's surface.

This module provides functionality to calculate the shortest perpendicular distance
from a geographic point to a line segment defined by two endpoints, using 3D
spherical geometry for accuracy.
"""

import numpy as np
from typing import Union

from pyearth.system.define_global_variables import earth_radius
from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
)


def calculate_distance_to_line(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> float:
    """
    Calculate the shortest distance from a point to a line segment on Earth's surface.

    This function computes the perpendicular distance from a geographic point to a
    line segment defined by two endpoints. It uses 3D Cartesian coordinates on a
    unit sphere for accurate calculations, then scales the result to meters using
    Earth's radius.

    Args:
        dLongitude1_in: Longitude of line segment start point (degrees or radians)
        dLatitude1_in: Latitude of line segment start point (degrees or radians)
        dLongitude2_in: Longitude of query point (degrees or radians)
        dLatitude2_in: Latitude of query point (degrees or radians)
        dLongitude3_in: Longitude of line segment end point (degrees or radians)
        dLatitude3_in: Latitude of line segment end point (degrees or radians)
        iFlag_radian: If True, input coordinates are in radians; if False, in degrees (default: False)

    Returns:
        float: Shortest distance from the query point (point 2) to the line segment
               (point 1 to point 3) in meters

    Note:
        - Point 2 is the query point whose distance is being measured
        - Points 1 and 3 define the line segment endpoints
        - The distance is perpendicular to the line segment
        - If the perpendicular projection falls outside the segment, the distance
          to the nearest endpoint is returned
        - Uses Earth's radius from global variables

    Example:
        >>> # Calculate distance from San Francisco to LA-Vegas line segment
        >>> distance = calculate_distance_to_line(
        ...     -118.24, 34.05,   # Point 1: Los Angeles
        ...     -122.42, 37.77,   # Point 2: San Francisco (query point)
        ...     -115.14, 36.17    # Point 3: Las Vegas
        ... )
        >>> print(f"Distance: {distance/1000:.2f} km")

    Algorithm:
        1. Convert geographic coordinates to 3D Cartesian (x,y,z) on unit sphere
        2. Create vectors:
           - line_vec: from point 1 to point 3 (defines the line segment)
           - point_vec: from point 1 to point 2 (query point)
        3. Project point_vec onto line_vec to find closest point on infinite line
        4. Clamp projection to stay within line segment [0, segment_length]
        5. Calculate Euclidean distance from query point to projection point
        6. Scale by Earth's radius to get distance in meters

    References:
        - Distance from point to line: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
        - Vector projection: https://en.wikipedia.org/wiki/Vector_projection
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
    )  # Line start
    x2, y2, z2 = convert_longitude_latitude_to_sphere_3d(
        lon2_rad, lat2_rad
    )  # Query point
    x3, y3, z3 = convert_longitude_latitude_to_sphere_3d(lon3_rad, lat3_rad)  # Line end

    # Create vectors for distance calculation
    line_vector = np.array([x3 - x1, y3 - y1, z3 - z1])  # Line segment vector
    point_vector = np.array([x2 - x1, y2 - y1, z2 - z1])  # Vector to query point

    # Calculate line segment length and unit vector
    line_length = np.linalg.norm(line_vector)

    # Handle edge case: line segment has zero length (start == end)
    if line_length < 1e-10:
        # Distance to the single point
        distance_on_sphere = np.linalg.norm(point_vector)
    else:
        line_unit_vector = line_vector / line_length

        # Project query point onto the line to find closest point on infinite line
        projection_length = np.dot(point_vector, line_unit_vector)

        # Clamp projection to line segment bounds [0, line_length]
        # This ensures we find the closest point ON THE SEGMENT, not the infinite line
        projection_length_clamped = np.clip(projection_length, 0.0, line_length)

        # Find the actual projection point on the line segment
        line_start_point = np.array([x1, y1, z1])
        projection_point = (
            line_start_point + projection_length_clamped * line_unit_vector
        )

        # Calculate Euclidean distance from query point to projection point
        query_point = np.array([x2, y2, z2])
        distance_on_sphere = np.linalg.norm(query_point - projection_point)

    # Convert from unit sphere distance to meters using Earth's radius
    distance_meters = distance_on_sphere * earth_radius

    return float(distance_meters)
