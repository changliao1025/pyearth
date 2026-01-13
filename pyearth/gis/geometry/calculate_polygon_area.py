"""
Spherical polygon area calculations on Earth's surface.

This module provides functions for calculating the area of polygons on a
spherical Earth model. Three different algorithms are available, each with
different trade-offs between accuracy, performance, and implementation complexity.

Main Functions
--------------
calculate_polygon_area : Calculate area of spherical polygon with multiple algorithms
calculate_polygon_file_area : Calculate total area from polygon file (GeoJSON)
spherical_polygon_area : Karney's method for spherical polygon area (most accurate)
haversine : Haversine function for spherical trigonometry

Algorithms Available
--------------------
Algorithm 0: Green's Theorem Line Integral
    - Based on Green's theorem using line integrals
    - Fast and memory efficient
    - Good for regular polygons

Algorithm 1: L'Huilier's Theorem (Spherical Triangulation)
    - Decomposes polygon into triangles from first vertex
    - Uses L'Huilier's formula for spherical excess
    - Good for convex polygons

Algorithm 2: Karney's Method (Default, Recommended)
    - Most accurate for all polygon shapes
    - Handles edge cases and degenerate polygons
    - Based on JPL/NASA research
    - Reference: https://trs.jpl.nasa.gov/handle/2014/41271

Use Cases
---------
- Geographic area calculation for land parcels, watersheds, countries
- Mesh generation and quality metrics
- Geospatial analysis and statistics
- Validation of polygon simplification algorithms

Notes
-----
- All algorithms assume a perfect sphere (not WGS84 ellipsoid)
- For ellipsoidal calculations, consider using nvector API
- Input coordinates should be in degrees unless iFlag_radian=True
- Polygons are automatically closed if not already closed
- Degenerate polygons (line-like shapes) can be detected with dLine_threshold

See Also
--------
calculate_spherical_triangle_area : Calculate area of single spherical triangle
visvalingam_whyatt_geodetic : Polygon simplification using area metrics

References
----------
.. [1] Beyer, W.H. "CRC Standard Mathematical Tables", 28th ed.
       CRC Press, 1987, p. 132.
.. [2] Karney, C.F.F. "Algorithms for geodesics", 2013.
       https://trs.jpl.nasa.gov/handle/2014/41271
.. [3] L'Huilier, S.-A.-J. "Mémoire sur la polyèdrométrie", 1812.
"""

import math
import numpy as np
from typing import Union, Optional, Tuple
from osgeo import ogr
from pyearth.system.define_global_variables import earth_radius

from pyearth.gis.geometry.calculate_spherical_triangle_area import (
    calculate_spherical_triangle_area,
)
from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
    calculate_distance_based_on_longitude_latitude,
)


def haversine(x: float) -> float:
    """Calculate the haversine function: hav(x) = (1 - cos(x)) / 2.

    The haversine function is commonly used in spherical trigonometry
    calculations, particularly for great circle distances and areas.

    Parameters
    ----------
    x : float
        Angle in radians.

    Returns
    -------
    float
        The haversine of x, in range [0, 1].
        - hav(0) = 0
        - hav(π) = 1
        - hav(π/2) = 0.5

    Notes
    -----
    This is equivalent to: sin²(x/2)

    The haversine function was historically important for navigation
    before the advent of electronic calculators, as it avoids issues
    with numerical precision when calculating small angles.

    Examples
    --------
    >>> haversine(0)
    0.0
    >>> haversine(np.pi)
    1.0
    >>> np.isclose(haversine(np.pi/2), 0.5)
    True

    See Also
    --------
    spherical_polygon_area : Uses haversine in Karney's algorithm
    """
    return (1.0 - math.cos(x)) / 2.0


def calculate_polygon_area(
    aLongitude_in: Union[list, np.ndarray],
    aLatitude_in: Union[list, np.ndarray],
    iFlag_algorithm: int = 2,
    iFlag_radian: bool = False,
    dRadius_in: Optional[float] = None,
    dLine_threshold: Optional[float] = None,
) -> float:
    """Calculate area of a spherical polygon on Earth's surface.

    Computes the area of a polygon on a spherical Earth using one of three
    available algorithms. The polygon is automatically closed if the first
    and last points differ.

    Parameters
    ----------
    aLongitude_in : list or np.ndarray
        Longitude coordinates of polygon vertices.
        In degrees by default (or radians if iFlag_radian=True).
    aLatitude_in : list or np.ndarray
        Latitude coordinates of polygon vertices.
        In degrees by default (or radians if iFlag_radian=True).
        Must have same length as aLongitude_in.
    iFlag_algorithm : int, optional
        Algorithm selection (default: 2, recommended):
        - 0: Green's Theorem line integral (fast, good for regular polygons)
        - 1: L'Huilier's theorem triangulation (good for convex polygons)
        - 2: Karney's method (most accurate, handles edge cases)
    iFlag_radian : bool, optional
        If True, input coordinates are in radians. If False (default),
        coordinates are in degrees.
    dRadius_in : float, optional
        Sphere radius in meters. If None (default), uses Earth's mean
        radius from global variables (approximately 6371229 m).
    dLine_threshold : float, optional
        Threshold for detecting degenerate line-like polygons.
        If max_edge_length / sum_other_edges > (1 - threshold),
        the polygon is considered degenerate and area = 0.
        If None, no check is performed.

    Returns
    -------
    float
        Polygon area in square meters (or square units of dRadius_in).
        Returns 0.0 for degenerate polygons when dLine_threshold is set.

    Raises
    ------
    ValueError
        If less than 3 points are provided (not a valid polygon).
    ValueError
        If aLongitude_in and aLatitude_in have different lengths.

    Notes
    -----
    - The polygon is automatically closed if not already closed
    - All algorithms assume a perfect sphere, not the WGS84 ellipsoid
    - For more accurate ellipsoidal calculations, consider using nvector API
    - Algorithm 2 (Karney) is recommended for general use
    - The area is always positive (uses absolute value internally)

    Examples
    --------
    >>> # Square on equator (approx 111km × 111km)
    >>> lon = [0, 1, 1, 0]
    >>> lat = [0, 0, 1, 1]
    >>> area = calculate_polygon_area(lon, lat, iFlag_algorithm=2)
    >>> # Result: approximately 12,364 km²

    >>> # Triangle with custom radius
    >>> lon = [0, 1, 0.5]
    >>> lat = [0, 0, 1]
    >>> area = calculate_polygon_area(lon, lat, dRadius_in=1000.0)

    >>> # Detect degenerate polygon
    >>> lon = [0, 100, 0.01]  # Nearly a line
    >>> lat = [0, 0, 0]
    >>> area = calculate_polygon_area(lon, lat, dLine_threshold=0.01)
    >>> # Returns 0.0 if polygon is too narrow

    Algorithm Details
    -----------------
    Algorithm 0 (Green's Theorem):
        Based on line integral formulation. Computes colatitudes and
        azimuths, then integrates (1-cos(colat)) * daz. Fast but may
        have precision issues for very small or irregular polygons.

    Algorithm 1 (L'Huilier):
        Triangulates polygon from first vertex, computes spherical
        excess for each triangle using L'Huilier's formula, then sums.
        Good accuracy for convex polygons.

    Algorithm 2 (Karney):
        Uses Karney's spherical polygon area formula from JPL.
        Most robust and accurate, especially for edge cases like
        polygons crossing poles or anti-meridian.

    See Also
    --------
    calculate_polygon_file_area : Calculate area from GeoJSON file
    calculate_spherical_triangle_area : Area of single spherical triangle
    spherical_polygon_area : Karney's algorithm implementation

    References
    ----------
    .. [1] Karney, C.F.F. "Algorithms for geodesics", 2013.
           https://trs.jpl.nasa.gov/handle/2014/41271
    .. [2] https://mathworld.wolfram.com/SphericalPolygon.html
    """
    # Convert to numpy arrays for easier manipulation
    aLongitude_in = np.asarray(aLongitude_in)
    aLatitude_in = np.asarray(aLatitude_in)

    # Validate inputs
    npoint = len(aLongitude_in)
    if npoint < 3:
        raise ValueError(
            f"A polygon requires at least 3 points. Got {npoint} point(s)."
        )

    if len(aLatitude_in) != npoint:
        raise ValueError(
            f"Longitude and latitude arrays must have same length. "
            f"Got {npoint} longitudes and {len(aLatitude_in)} latitudes."
        )

    if len(aLatitude_in) != npoint:
        raise ValueError(
            f"Longitude and latitude arrays must have same length. "
            f"Got {npoint} longitudes and {len(aLatitude_in)} latitudes."
        )

    # Close polygon if not already closed
    if aLatitude_in[0] != aLatitude_in[-1] or aLongitude_in[0] != aLongitude_in[-1]:
        aLatitude_in = np.append(aLatitude_in, aLatitude_in[0])
        aLongitude_in = np.append(aLongitude_in, aLongitude_in[0])
        npoint = len(aLongitude_in)

    # Convert to radians if needed
    if iFlag_radian:
        aLongitude_radian_in = aLongitude_in
        aLatitude_radian_in = aLatitude_in
    else:
        aLongitude_radian_in = np.radians(aLongitude_in)
        aLatitude_radian_in = np.radians(aLatitude_in)

    # Optional: Check for degenerate line-like polygons
    if dLine_threshold is not None:
        aLength = np.zeros(npoint - 1)
        for i in range(npoint - 1):
            dLength = calculate_distance_based_on_longitude_latitude(
                aLongitude_in[i],
                aLatitude_in[i],
                aLongitude_in[i + 1],
                aLatitude_in[i + 1],
            )
            aLength[i] = dLength

        dLength_max = np.max(aLength)
        dlength_rest = np.sum(aLength) - dLength_max

        # Check if polygon is too narrow (close to a line)
        if dlength_rest > 0 and (dLength_max / dlength_rest) > (1 - dLine_threshold):
            # Degenerate polygon - essentially a line
            return 0.0

    # Algorithm selection
    if iFlag_algorithm == 0:
        # Algorithm 0: Green's Theorem Line Integral
        # Get colatitude (surface distance as angle from origin)
        a = (
            np.sin(aLatitude_radian_in / 2) ** 2
            + np.cos(aLatitude_radian_in) * np.sin(aLongitude_radian_in / 2) ** 2
        )
        colat = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

        # Azimuth of each point from arbitrary origin
        az = np.arctan2(
            np.cos(aLatitude_radian_in) * np.sin(aLongitude_radian_in),
            np.sin(aLatitude_radian_in),
        ) % (2 * np.pi)

        # Calculate azimuth differences
        daz = np.diff(az)
        daz = (daz + np.pi) % (2 * np.pi) - np.pi

        # Average surface distance for each step
        deltas = np.diff(colat) / 2
        colat = colat[0:-1] + deltas

        # Integral: (1 - cos(colat)) * daz
        integrands = (1 - np.cos(colat)) * daz

        # Sum and take absolute value
        area = abs(np.sum(integrands))

        # Could be area of inside or outside - take minimum
        area = min(area, 1 - area)

    elif iFlag_algorithm == 1:
        # Algorithm 1: L'Huilier's Theorem (Spherical Triangulation)
        dLongitude_root = aLongitude_radian_in[0]
        dLatitude_root = aLatitude_radian_in[0]
        dLongitude_b = aLongitude_radian_in[1]
        dLatitude_b = aLatitude_radian_in[1]

        nTriangle = npoint - 2
        aArea = np.zeros(nTriangle)

        for i in np.arange(1, nTriangle + 1, 1):
            # Define triangle vertices
            dLongitude_a = dLongitude_b
            dLatitude_a = dLatitude_b
            dLongitude_b = aLongitude_radian_in[i + 1]
            dLatitude_b = aLatitude_radian_in[i + 1]

            # Calculate spherical triangle area
            aLongitude_temp = [dLongitude_root, dLongitude_a, dLongitude_b]
            aLatitude_temp = [dLatitude_root, dLatitude_a, dLatitude_b]
            dArea_triangle = calculate_spherical_triangle_area(
                aLongitude_temp, aLatitude_temp, iFlag_radian=True
            )
            aArea[i - 1] = dArea_triangle

        area = np.sum(aArea)

    elif iFlag_algorithm == 2:
        # Algorithm 2: Karney's Method (JPL/NASA - Most Accurate)
        if dRadius_in is not None:
            dRadius = dRadius_in
        else:
            dRadius = earth_radius

        dArea_m = spherical_polygon_area(
            aLatitude_radian_in, aLongitude_radian_in, dRadius
        )
        return float(dArea_m)

    else:
        raise ValueError(
            f"Invalid algorithm selection: {iFlag_algorithm}. "
            "Valid options are 0 (Green), 1 (L'Huilier), or 2 (Karney)."
        )

    # Convert from fraction of sphere to square meters
    if iFlag_radian:
        # Return as fraction of sphere
        return float(area)
    else:
        # Convert to square meters
        if dRadius_in is not None:
            dArea_m = area * dRadius_in**2
        else:
            dArea_m = area * earth_radius**2

        return float(dArea_m)


def calculate_polygon_file_area(sFilename_polygon_in: str) -> float:
    """Calculate total area of all polygons in a GeoJSON file.

    Reads a GeoJSON file and computes the sum of areas for all polygon
    features. Supports both POLYGON and MULTIPOLYGON geometry types.

    Parameters
    ----------
    sFilename_polygon_in : str
        Path to the GeoJSON file containing polygon geometries.

    Returns
    -------
    float
        Total area of all polygons in square meters.
        Returns 0.0 if no valid polygons are found.

    Raises
    ------
    FileNotFoundError
        If the specified file cannot be opened.

    Notes
    -----
    - Only POLYGON and MULTIPOLYGON geometries are processed
    - Empty geometries are skipped
    - Uses algorithm 2 (Karney) for all area calculations
    - Requires GDAL/OGR to be installed

    Examples
    --------
    >>> # Calculate total area from a GeoJSON file
    >>> area = calculate_polygon_file_area('watersheds.geojson')
    >>> print(f"Total area: {area / 1e6:.2f} km²")

    See Also
    --------
    calculate_polygon_area : Calculate area of single polygon
    get_geometry_coordinates : Extract coordinates from OGR geometry
    """
    from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates

    pDriver = ogr.GetDriverByName("GeoJSON")
    pDataSource = pDriver.Open(sFilename_polygon_in, 0)

    if pDataSource is None:
        raise FileNotFoundError(
            f"Could not open file: {sFilename_polygon_in}. "
            "Please check the file path and ensure it's a valid GeoJSON file."
        )

    pLayer = pDataSource.GetLayer()
    dArea = 0.0

    pFeature = pLayer.GetNextFeature()
    while pFeature:
        pGeometry = pFeature.GetGeometryRef()

        if pGeometry is None:
            pFeature = pLayer.GetNextFeature()
            continue

        if pGeometry.IsEmpty():
            pFeature = pLayer.GetNextFeature()
            continue

        sGeometryType = pGeometry.GetGeometryName()

        if sGeometryType == "POLYGON":
            aCoords_gcs = get_geometry_coordinates(pGeometry)
            dArea += calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1])
        elif sGeometryType == "MULTIPOLYGON":
            for i in range(pGeometry.GetGeometryCount()):
                pGeometry_temp = pGeometry.GetGeometryRef(i)
                aCoords_gcs = get_geometry_coordinates(pGeometry_temp)
                dArea += calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1])
        # Skip unsupported geometry types

        pFeature = pLayer.GetNextFeature()

    pDataSource = None
    return dArea


def spherical_polygon_area(
    lat: Union[list, np.ndarray], lon: Union[list, np.ndarray], r: float
) -> float:
    """Calculate area of spherical polygon using Karney's method.

    This is the most accurate algorithm for spherical polygon area calculation,
    based on research from JPL/NASA. It handles edge cases well including
    polygons that cross the anti-meridian or contain poles.

    The algorithm computes the spherical excess using a robust formulation
    based on L'Huilier's theorem applied to each edge's spherical triangle.

    Parameters
    ----------
    lat : list or np.ndarray
        Latitudes of all vertices in radians.
        Must be in range [-π/2, π/2].
    lon : list or np.ndarray
        Longitudes of all vertices in radians.
        Must be in range [-π, π] or [0, 2π].
    r : float
        Spherical radius in meters (typically Earth's mean radius ≈ 6371229 m).

    Returns
    -------
    float
        Area of the spherical polygon in square units of r.
        Always returns a positive value.

    Notes
    -----
    - Input coordinates must be in radians
    - Polygon should be closed (first point = last point)
    - The algorithm is numerically stable even for very small or large polygons
    - Handles degenerate cases (e.g., consecutive identical vertices)

    Algorithm
    ---------
    For each edge of the polygon:
    1. Compute the spherical triangle formed by:
       - The two edge vertices
       - The North pole (or origin point)
    2. Calculate the spherical excess using L'Huilier's formula:
       E = 4 * arctan(sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)))
       where s = (a+b+c)/2 is the semi-perimeter
    3. Sum the signed excesses (sign depends on edge direction)
    4. Multiply by r² to get area

    Examples
    --------
    >>> # Spherical cap at North pole (small circle at 45°N)
    >>> lats = np.radians([45, 45, 45, 45])
    >>> lons = np.radians([0, 90, 180, 270])
    >>> r = 6371229.0  # Earth's mean radius
    >>> area = spherical_polygon_area(lats, lons, r)

    See Also
    --------
    calculate_polygon_area : High-level interface with algorithm selection
    haversine : Haversine function used in this algorithm

    References
    ----------
    .. [1] Karney, C.F.F. "Algorithms for geodesics", Journal of Geodesy, 2013.
           https://trs.jpl.nasa.gov/handle/2014/41271
    .. [2] Beyer, W.H. "CRC Standard Mathematical Tables", 28th ed.
    """
    lam1 = lam2 = beta1 = beta2 = cosB1 = cosB2 = 0.0
    hav = 0.0
    sum_excess = 0.0

    for j in range(len(lat)):
        k = (j + 1) % len(lat)

        if j == 0:
            lam1 = lon[j]
            beta1 = lat[j]
            lam2 = lon[j + 1]
            beta2 = lat[j + 1]
            cosB1 = math.cos(beta1)
            cosB2 = math.cos(beta2)
        else:
            # Reuse previous endpoint as new starting point
            lam1 = lam2
            beta1 = beta2
            lam2 = lon[k]
            beta2 = lat[k]
            cosB1 = cosB2
            cosB2 = math.cos(beta2)

        # Skip if edge has zero length (same longitude)
        if lam1 != lam2:
            # Calculate haversine distance between the two points
            hav = haversine(beta2 - beta1) + cosB1 * cosB2 * haversine(lam2 - lam1)
            a = 2 * math.asin(math.sqrt(hav))

            # Compute complementary latitudes (colatitudes)
            b = math.pi / 2 - beta2
            c = math.pi / 2 - beta1

            # Semi-perimeter of spherical triangle
            s = 0.5 * (a + b + c)

            # L'Huilier's formula for spherical excess
            # tan(E/4) = sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2))
            t = (
                math.tan(s / 2)
                * math.tan((s - a) / 2)
                * math.tan((s - b) / 2)
                * math.tan((s - c) / 2)
            )

            excess = abs(4 * math.atan(math.sqrt(abs(t))))

            # Sign of excess depends on edge direction
            if lam2 < lam1:
                excess = -excess

            sum_excess += excess

    # Area = |spherical excess| * r²
    return abs(sum_excess) * r * r
