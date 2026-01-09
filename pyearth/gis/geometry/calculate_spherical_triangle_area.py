"""
Calculate the area of a spherical triangle using L'Huilier's theorem.

This module provides utilities for computing the area of a triangle on a
spherical surface, such as Earth. The calculation uses L'Huilier's formula,
which is a robust method for calculating spherical excess.

Main Function
-------------
calculate_spherical_triangle_area : Compute area of spherical triangle

Algorithm
---------
L'Huilier's Theorem computes the spherical excess (E) of a triangle, which
is directly related to the triangle's area:

    Area = E * R²

where:
    - E is the spherical excess in radians
    - R is the sphere radius
    - E = 4 * arctan(sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)))
    - s = (a + b + c) / 2 (semi-perimeter)
    - a, b, c are the great circle arc lengths between vertices

Use Cases
---------
- Polygon area calculation (triangulation method)
- Mesh quality metrics on spherical surfaces
- Geographic area computation
- Spherical geometry calculations

See Also
--------
calculate_polygon_area : Calculate area of spherical polygon
calculate_distance_based_on_longitude_latitude : Great circle distance

References
----------
.. [1] L'Huilier, S.-A.-J. "Mémoire sur la polyèdrométrie", 1812.
.. [2] https://mathworld.wolfram.com/SphericalTriangle.html
.. [3] https://mathworld.wolfram.com/LHuiliersTheorem.html
"""

import numpy as np
from typing import Union, List, Optional
from pyearth.system.define_global_variables import earth_radius
from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
    calculate_distance_based_on_longitude_latitude,
)


def calculate_spherical_triangle_area(
    aLongitude_in: Union[List[float], np.ndarray],
    aLatitude_in: Union[List[float], np.ndarray],
    iFlag_radian: bool = False,
    dRadius_in: Optional[float] = None,
) -> float:
    """Calculate the area of a spherical triangle using L'Huilier's theorem.

    Computes the area of a triangle on a spherical surface given the
    longitude/latitude coordinates of its three points. Uses L'Huilier's
    formula for calculating the spherical excess, which is numerically
    stable for triangles of all sizes.

    Parameters
    ----------
    aLongitude_in : list or np.ndarray
        Longitudes of the three triangle points.
        Must contain exactly 3 values.
        In degrees by default, or radians if iFlag_radian=True.
    aLatitude_in : list or np.ndarray
        Latitudes of the three triangle points.
        Must contain exactly 3 values.
        In degrees by default, or radians if iFlag_radian=True.
    iFlag_radian : bool, optional
        If True, input coordinates are in radians and output is spherical
        excess in radians. If False (default), input is in degrees and
        output is area in square meters (or square units of dRadius_in).
    dRadius_in : float, optional
        Sphere radius in meters. If None (default), uses Earth's mean
        radius from global variables (approximately 6371229 m).
        Only used when iFlag_radian=False.

    Returns
    -------
    float
        If iFlag_radian=True: Spherical excess in radians.
        If iFlag_radian=False: Triangle area in square meters
        (or square units of dRadius_in).

    Raises
    ------
    ValueError
        If aLongitude_in or aLatitude_in don't contain exactly 3 values.
        If the three points are collinear (degenerate triangle).

    Notes
    -----
    - The algorithm uses L'Huilier's theorem for numerical stability
    - All calculations are performed on a perfect sphere
    - The three points must not be collinear
    - Triangle points should be ordered (clockwise or counterclockwise)

    Algorithm
    ---------
     1. Calculate great circle distances between each pair of points:
         a = distance(point0, point1)
         b = distance(point1, point2)
         c = distance(point2, point0)

    2. Compute semi-perimeter:
       s = (a + b + c) / 2

    3. Apply L'Huilier's formula for spherical excess:
       tan(E/4) = sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2))
       E = 4 * arctan(tan(E/4))

    4. Calculate area:
       Area = E * R²

    Examples
    --------
    >>> # Equilateral triangle on Earth's surface
    >>> # Vertices at (0°, 0°), (1°, 0°), (0.5°, 0.866°)
    >>> lons = [0, 1, 0.5]
    >>> lats = [0, 0, 0.866]
    >>> area = calculate_spherical_triangle_area(lons, lats)
    >>> # Returns area in square meters

    >>> # Using radians and getting spherical excess
    >>> import numpy as np
    >>> lons_rad = np.radians([0, 10, 5])
    >>> lats_rad = np.radians([0, 0, 8.66])
    >>> excess = calculate_spherical_triangle_area(
    ...     lons_rad, lats_rad, iFlag_radian=True
    ... )
    >>> # Returns spherical excess in radians

    >>> # Custom sphere radius
    >>> area = calculate_spherical_triangle_area(
    ...     [0, 1, 0.5], [0, 0, 1], dRadius_in=1000.0
    ... )
    >>> # Returns area in square units of radius 1000

    See Also
    --------
    calculate_polygon_area : Calculate area of multi-sided spherical polygon
    calculate_distance_based_on_longitude_latitude : Great circle distance

    References
    ----------
    .. [1] L'Huilier, S.-A.-J. "Mémoire sur la polyèdrométrie", 1812.
    .. [2] Girard's theorem: https://mathworld.wolfram.com/GirardsSphericalExcessFormula.html
    .. [3] L'Huilier's theorem: https://mathworld.wolfram.com/LHuiliersTheorem.html
    """
    # Validate inputs
    aLongitude_in = np.asarray(aLongitude_in)
    aLatitude_in = np.asarray(aLatitude_in)

    if len(aLongitude_in) != 3:
        raise ValueError(
            f"Triangle requires exactly 3 vertices. "
            f"Got {len(aLongitude_in)} longitude values."
        )

    if len(aLatitude_in) != 3:
        raise ValueError(
            f"Triangle requires exactly 3 vertices. "
            f"Got {len(aLatitude_in)} latitude values."
        )

    # Convert to radians if needed
    if iFlag_radian:
        aLongitude_radian = aLongitude_in
        aLatitude_radian = aLatitude_in
    else:
        aLongitude_radian = np.radians(aLongitude_in)
        aLatitude_radian = np.radians(aLatitude_in)

    # Calculate the three great circle arc lengths (edges of the triangle)
    # Edge a: from point 0 to point 1
    a = calculate_distance_based_on_longitude_latitude(
        aLongitude_radian[0],
        aLatitude_radian[0],
        aLongitude_radian[1],
        aLatitude_radian[1],
        iFlag_radian=True,
    )

    # Edge b: from point 1 to point 2
    b = calculate_distance_based_on_longitude_latitude(
        aLongitude_radian[1],
        aLatitude_radian[1],
        aLongitude_radian[2],
        aLatitude_radian[2],
        iFlag_radian=True,
    )

    # Edge c: from point 2 to point 0
    c = calculate_distance_based_on_longitude_latitude(
        aLongitude_radian[2],
        aLatitude_radian[2],
        aLongitude_radian[0],
        aLatitude_radian[0],
        iFlag_radian=True,
    )

    # Check for degenerate triangle (collinear points)
    if np.isclose(a + b, c) or np.isclose(b + c, a) or np.isclose(a + c, b):
        raise ValueError(
            "Triangle points are collinear (degenerate triangle). "
            f"Edge lengths: a={a}, b={b}, c={c}"
        )

    # Calculate semi-perimeter
    s = 0.5 * (a + b + c)

    # L'Huilier's formula for spherical excess
    # tan(E/4) = sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2))
    tan_quarter_excess = np.sqrt(
        np.tan(0.5 * s)
        * np.tan(0.5 * (s - a))
        * np.tan(0.5 * (s - b))
        * np.tan(0.5 * (s - c))
    )

    # Spherical excess in radians
    spherical_excess = 4.0 * np.arctan(tan_quarter_excess)

    # Return based on flag
    if iFlag_radian:
        # Return spherical excess in radians
        return float(spherical_excess)
    else:
        # Calculate area in square meters (or square units of radius)
        if dRadius_in is not None:
            area = spherical_excess * dRadius_in**2
        else:
            area = spherical_excess * earth_radius**2

        return float(area)
