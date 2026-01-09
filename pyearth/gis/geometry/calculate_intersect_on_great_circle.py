"""
Great circle calculations and geometric operations on a sphere.

This module provides utilities for working with great circles (geodesic paths)
on a sphere. Great circles represent the shortest distance between two points
on a sphere and are fundamental in geospatial calculations.

Main Functions
--------------
calculate_intersect_on_great_circle : Find closest point on great circle arc to query point
find_great_circle_intersection_with_meridian : Find latitude where great circle crosses a meridian
project_point_onto_plane : Project 3D point onto plane


Use Cases
---------
- Mesh operations: Finding nearest point on a line to a point
- Polygon splitting: Determining where polygons cross the International Date Line
- Distance calculations: Computing shortest paths on Earth's surface
- Geometric operations: Projections and intersections on spherical surfaces

Notes
-----
All calculations assume a perfect sphere (Earth's ellipsoidal shape is ignored).
Input/output angles are in degrees unless specified otherwise.
The convert_sphere_3d_to_longitude_latitude function is imported from
convert_between_longitude_latitude_and_sphere_3d module.

See Also
--------
calculate_distance_to_line : Calculate perpendicular distance to line segment
calculate_distance_to_plane : Calculate distance from point to plane
convert_between_longitude_latitude_and_sphere_3d : Coordinate conversion utilities

References
----------
.. [1] Weisstein, Eric W. "Great Circle." From MathWorld--A Wolfram Web Resource.
       https://mathworld.wolfram.com/GreatCircle.html
.. [2] Williams, E. "Aviation Formulary V1.47"
       http://www.edwilliams.org/avform147.htm
"""

import numpy as np
from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
    convert_sphere_3d_to_longitude_latitude,
)


def project_point_onto_plane(p: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Project a 3D point onto a plane defined by its normal vector.

    Given a point and a plane (defined by a normal vector passing through
    the origin), this function calculates the orthogonal projection of the
    point onto the plane.

    Parameters
    ----------
    p : np.ndarray
        3D point to be projected, shape (3,).
        Example: [x, y, z]
    normal : np.ndarray
        Normal vector of the plane, shape (3,).
        Should be a unit vector (magnitude = 1) for correct results.
        The plane passes through the origin with this normal direction.

    Returns
    -------
    np.ndarray
        Projected point on the plane, shape (3,).
        This is the point on the plane closest to the input point.

    Raises
    ------
    ValueError
        If normal vector has zero magnitude.

    Notes
    -----
    - The plane is defined by all points r such that: normal · r = 0
      (plane through origin perpendicular to normal)
    - Projection formula: p_proj = p - (p · n) * n, where n is unit normal
    - If normal is not unit length, it will be normalized automatically
    - Distance from point to plane is: |p · n| (when n is unit vector)

    Examples
    --------
    >>> # Project point onto xy-plane (normal = z-axis)
    >>> p = np.array([1, 2, 3])
    >>> normal = np.array([0, 0, 1])
    >>> proj = project_point_onto_plane(p, normal)
    >>> np.allclose(proj, [1, 2, 0])
    True

    See Also
    --------
    calculate_distance_to_plane : Calculate distance from point to plane
    """
    # Validate and normalize the normal vector
    normal_magnitude = np.linalg.norm(normal)

    if normal_magnitude < 1e-10:
        raise ValueError(
            "Normal vector has zero or near-zero magnitude. "
            f"Got magnitude: {normal_magnitude}. Cannot define a plane."
        )

    # Normalize to unit vector if not already
    if not np.isclose(normal_magnitude, 1.0):
        normal_unit = normal / normal_magnitude
    else:
        normal_unit = normal

    # Calculate the projection
    # Distance from point to plane along normal direction
    distance_to_plane = np.dot(p, normal_unit)

    # Subtract the normal component to get projection
    projection = p - distance_to_plane * normal_unit

    return projection


def calculate_intersect_on_great_circle(
    dLongitude1_in: float,
    dLatitude1_in: float,
    dLongitude2_in: float,
    dLatitude2_in: float,
    dLongitude3_in: float,
    dLatitude3_in: float,
    iFlag_radian: bool = False,
) -> tuple:
    """Calculate the closest point on a great circle arc to a query point.

    This function finds the point on the great circle arc between point1 and
    point3 that is closest to point2 (the query point). This is useful for
    finding the nearest location on a line to a given point in mesh operations.

    Parameters
    ----------
    dLongitude1_in : float
        Longitude of first point on the great circle arc.
        In degrees by default, or radians if iFlag_radian=True.
    dLatitude1_in : float
        Latitude of first point on the great circle arc.
        In degrees by default, or radians if iFlag_radian=True.
    dLongitude2_in : float
        Longitude of query point (point to find nearest location to).
        In degrees by default, or radians if iFlag_radian=True.
    dLatitude2_in : float
        Latitude of query point (point to find nearest location to).
        In degrees by default, or radians if iFlag_radian=True.
    dLongitude3_in : float
        Longitude of second point on the great circle arc.
        In degrees by default, or radians if iFlag_radian=True.
    dLatitude3_in : float
        Latitude of second point on the great circle arc.
        In degrees by default, or radians if iFlag_radian=True.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians. If False (default),
        input coordinates are in degrees.

    Returns
    -------
    tuple
        (longitude, latitude) of the closest point on the great circle arc.
        Always returned in degrees.

    Notes
    -----
        - The returned point is the orthogonal projection onto the great circle,
            which may lie outside the arc segment between point1 and point3
        - All calculations are performed on a unit sphere
        - The great circle is defined by the plane through point1, point3,
            and Earth's center

    Examples
    --------
    >>> # Find closest point on equator arc to North Pole
    >>> lon, lat = calculate_intersect_on_great_circle(
    ...     0, 0,      # Point 1: (0°E, 0°N) on equator
    ...     0, 90,     # Query: North Pole
    ...     90, 0      # Point 3: (90°E, 0°N) on equator
    ... )
    >>> # Result should be close to (0°, 0°) - nearest point on equator

    Algorithm
    ---------
    1. Convert all points to 3D Cartesian coordinates on unit sphere
    2. Calculate great circle plane normal: n = normalize(cross(point1, point3))
    3. Project query point onto plane: projected = query - dot(query, n) * n
    4. Normalize projected point to sphere surface
    5. Convert back to longitude/latitude

    See Also
    --------
    calculate_distance_to_line : Calculate perpendicular distance to line segment
    calculate_distance_to_plane : Calculate distance from point to plane
    """
    # Convert all points to 3D Cartesian coordinates
    point1_3d = convert_longitude_latitude_to_sphere_3d(
        dLongitude1_in, dLatitude1_in, iFlag_radian
    )
    query_point_3d = convert_longitude_latitude_to_sphere_3d(
        dLongitude2_in, dLatitude2_in, iFlag_radian
    )
    point3_3d = convert_longitude_latitude_to_sphere_3d(
        dLongitude3_in, dLatitude3_in, iFlag_radian
    )

    # Calculate the normal vector to the great circle plane
    # The plane passes through point1, point3, and Earth's center
    plane_normal = np.cross(point1_3d, point3_3d)
    plane_normal /= np.linalg.norm(plane_normal)

    # Project the query point onto the great circle plane
    distance_to_plane = np.dot(query_point_3d, plane_normal)
    projected_point = query_point_3d - distance_to_plane * plane_normal

    # Normalize to sphere surface (the projection may not be on unit sphere)
    projected_point /= np.linalg.norm(projected_point)

    # Convert back to longitude/latitude (always in degrees)
    longitude_intersect, latitude_intersect = convert_sphere_3d_to_longitude_latitude(
        *projected_point
    )

    return float(longitude_intersect), float(latitude_intersect)


def find_great_circle_intersection_with_meridian(
    lon1: float, lat1: float, lon2: float, lat2: float, target_lon: float
) -> tuple:
    """Find the latitude on a great circle arc at a specified longitude.

    Given two points defining a great circle, find the latitude where
    the great circle crosses a specified longitude meridian. This is commonly
    used for polygon splitting operations at the International Date Line (±180°).

    Parameters
    ----------
    lon1 : float
        Longitude of the first point in degrees.
    lat1 : float
        Latitude of the first point in degrees.
    lon2 : float
        Longitude of the second point in degrees.
    lat2 : float
        Latitude of the second point in degrees.
    target_lon : float
        The target longitude in degrees where intersection is sought.

    Returns
    -------
    tuple
        (target_lon, target_lat) where target_lat is the latitude at the
        target longitude.

    Raises
    ------
    ValueError
        If the great circle does not cross the target longitude meridian,
        or if the two input points define the same location.

    Notes
    -----
    - A great circle may cross a meridian at 0, 1, or 2 points
    - This function returns one intersection point
    - The target longitude is returned unchanged in the result tuple
    - All angles are in degrees

    Examples
    --------
    >>> # Great circle from (0°E, 0°N) to (180°E, 0°N) crosses 90°E at equator
    >>> lon, lat = calculate_great_circle_latitude_at_longitude(0, 0, 180, 0, 90)
    >>> np.allclose([lon, lat], [90, 0])
    True

    >>> # Great circle crossing International Date Line
    >>> lon, lat = calculate_great_circle_latitude_at_longitude(170, 45, -170, 50, 180)
    >>> # Returns latitude where circle crosses ±180°

    Algorithm
    ---------
    1. Convert both points to 3D Cartesian coordinates
    2. Calculate great circle plane normal: n = cross(p1, p2)
    3. Solve plane equation for latitude at target longitude:
       n_x * cos(lat) * cos(lon) + n_y * cos(lat) * sin(lon) + n_z * sin(lat) = 0
    4. Simplify to: tan(lat) = -(n_x*cos(lon) + n_y*sin(lon)) / n_z

    See Also
    --------
    calculate_intersect_on_great_circle : Find closest point on great circle
    split_polygon_cross_idl : Uses this function for polygon splitting

    References
    ----------
    .. [1] Williams, E. "Aviation Formulary V1.47"
           http://www.edwilliams.org/avform147.htm
    """
    # Convert to radians
    lon1_rad = np.radians(lon1)
    lat1_rad = np.radians(lat1)
    lon2_rad = np.radians(lon2)
    lat2_rad = np.radians(lat2)
    target_lon_rad = np.radians(target_lon)

    # Convert both points to 3D Cartesian coordinates
    p1 = convert_longitude_latitude_to_sphere_3d(lon1, lat1, iFlag_radian=False)
    p2 = convert_longitude_latitude_to_sphere_3d(lon2, lat2, iFlag_radian=False)

    # Calculate the great circle plane normal
    normal = np.cross(p1, p2)
    normal_magnitude = np.linalg.norm(normal)

    # Check if points are the same or antipodal
    if normal_magnitude < 1e-10:
        raise ValueError(
            f"Points ({lon1}, {lat1}) and ({lon2}, {lat2}) are too close "
            "or antipodal to define a unique great circle."
        )

    normal = normal / normal_magnitude

    # The great circle plane equation is: normal · r = 0
    # For a point at longitude=target_lon: r = [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]
    # Substituting into plane equation and solving for lat:
    # n_x * cos(lat) * cos(target_lon) + n_y * cos(lat) * sin(target_lon) + n_z * sin(lat) = 0
    # cos(lat) * (n_x * cos(target_lon) + n_y * sin(target_lon)) + n_z * sin(lat) = 0

    nx, ny, nz = normal
    cos_target = np.cos(target_lon_rad)
    sin_target = np.sin(target_lon_rad)

    # Coefficient of cos(lat)
    A = nx * cos_target + ny * sin_target
    # Coefficient of sin(lat)
    B = nz

    # We have: A * cos(lat) + B * sin(lat) = 0
    # This can be written as: tan(lat) = -A / B (if B != 0)

    if abs(B) < 1e-10:
        # Great circle plane is perpendicular to z-axis (equatorial great circle)
        if abs(A) < 1e-10:
            # Degenerate case - return equator point
            return (target_lon, 0.0)
        else:
            raise ValueError(
                f"Great circle does not intersect meridian at {target_lon}° "
                "(great circle is equatorial)"
            )

    # Calculate latitude
    lat_rad = np.arctan2(-A, B)
    lat_deg = np.degrees(lat_rad)

    return (target_lon, float(lat_deg))
