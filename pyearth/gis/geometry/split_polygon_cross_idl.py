"""
International Date Line Polygon Splitting
==========================================

This module provides functionality to split polygons that cross the International Date Line
(IDL) at 180°/-180° longitude into separate left (Eastern) and right (Western) hemisphere
polygons.

When polygons span across the International Date Line, they need to be split into separate
geometries for proper rendering and spatial operations in many GIS systems. This module
identifies where a polygon crosses the IDL and creates two valid closed polygons, one for
each hemisphere.

Key Features
------------
- Automatic detection of IDL crossing edges
- Great circle intersection calculation for accurate split points
- Separate left (Eastern, lon > 0) and right (Western, lon < 0) polygons
- Automatic polygon closure for resulting geometries
- Preserves topology and vertex order (CCW)

Main Functions
--------------
split_polygon_cross_idl : Split an IDL-crossing polygon into two hemisphere polygons

Algorithm Overview
------------------
The splitting algorithm works by:

1. **Find crossing edges**: Identify the two edges where longitude changes sign across 180°
2. **Calculate intersections**: Compute exact intersection points with the IDL meridian using
   great circle calculations
3. **Determine top/bottom**: Identify which intersection is at higher latitude
4. **Partition vertices**: Assign each vertex to left (Eastern) or right (Western) polygon
5. **Add boundary points**: Insert IDL boundary points (±180°) at intersections
6. **Close polygons**: Ensure both resulting polygons are properly closed rings

The algorithm assumes the polygon vertices are in counter-clockwise (CCW) order.

Use Cases
---------
- Splitting country boundaries that cross the IDL (e.g., Russia, Fiji, Kiribati)
- Processing Pacific Ocean regions for mapping
- Preparing polygons for renderers that don't handle IDL automatically
- Converting global datasets to hemisphere-specific projections
- Fixing topology issues in IDL-crossing geometries

Coordinate Systems
------------------
**Input/Output**: Geographic coordinates (longitude, latitude) in decimal degrees
    - Longitude: [-180, 180] range (negative = West, positive = East)
    - Latitude: [-90, 90] range (negative = South, positive = North)

**IDL Location**: 180° longitude (equivalent to -180°)
    - Eastern Hemisphere: (0, 180) → Left polygon uses ~179.9999999999°
    - Western Hemisphere: (-180, 0) → Right polygon uses ~-179.9999999999°

Limitations
-----------
- Assumes exactly TWO crossing edges (one polygon crossing the IDL)
- Requires counter-clockwise (CCW) vertex order
- Does not handle multi-part polygons or polygons with holes
- May fail for complex self-intersecting polygons

See Also
--------
pyearth.gis.geometry.reorder_idl_polygon : Alternative IDL handling approach
pyearth.gis.geometry.convert_idl_polygon_to_valid_polygon : Longitude shift approach
pyearth.gis.geometry.calculate_intersect_on_great_circle : Great circle intersection

References
----------
.. [1] International Date Line: https://en.wikipedia.org/wiki/International_Date_Line
.. [2] Great Circle Navigation: Intersection with meridian calculations

Examples
--------
>>> # Example 1: Simple rectangle crossing IDL
>>> coords = np.array([
...     [170.0, -10.0],
...     [170.0, 10.0],
...     [-170.0, 10.0],
...     [-170.0, -10.0],
...     [170.0, -10.0]  # Closed
... ])
>>> left, right = split_polygon_cross_idl(coords)
>>> # left polygon: Eastern hemisphere (170° to ~180°)
>>> # right polygon: Western hemisphere (~-180° to -170°)

>>> # Example 2: Complex polygon crossing IDL
>>> coords = np.array([
...     [175.0, 0.0],
...     [178.0, 5.0],
...     [-178.0, 5.0],
...     [-175.0, 0.0],
...     [-178.0, -5.0],
...     [178.0, -5.0],
...     [175.0, 0.0]
... ])
>>> left, right = split_polygon_cross_idl(coords)
"""

import os
import numpy as np
from typing import List, Tuple, Optional
from osgeo import ogr, gdal, osr
from pyearth.gis.geometry.calculate_intersect_on_great_circle import find_great_circle_intersection_with_meridian


def split_polygon_cross_idl(aCoord_gcs):
    """
    Split a polygon crossing the International Date Line into left and right hemisphere polygons.

    This function takes a polygon that crosses the IDL (180°/-180° longitude) and splits it
    into two separate polygons: one for the Eastern hemisphere (positive longitudes) and one
    for the Western hemisphere (negative longitudes). The split is performed along the IDL
    meridian using great circle intersection calculations.

    The algorithm identifies edges that cross the IDL, calculates exact intersection points,
    and partitions vertices into left (Eastern) and right (Western) polygons, adding boundary
    points at ±180° as needed.

    Parameters
    ----------
    aCoord_gcs : array-like
        Input polygon coordinates as Nx2 array of (longitude, latitude) pairs in decimal degrees.

        Requirements:
        - Should be a closed polygon (first point equals last point)
        - Vertices should be in counter-clockwise (CCW) order
        - Must cross the IDL exactly twice (one continuous crossing)
        - Minimum 4 points (3 unique + 1 closure)

        Format: [[lon1, lat1], [lon2, lat2], ..., [lonN, latN]]

        Longitude range: [-180, 180] degrees
        Latitude range: [-90, 90] degrees

    Returns
    -------
    list of numpy.ndarray
        List containing two polygons: [left_polygon, right_polygon]

        - **left_polygon** (numpy.ndarray): Eastern hemisphere polygon (lon > 0)
          - Coordinates with positive longitudes
          - Includes boundary points near +180°
          - Closed polygon (first == last)

        - **right_polygon** (numpy.ndarray): Western hemisphere polygon (lon < 0)
          - Coordinates with negative longitudes
          - Includes boundary points near -180°
          - Closed polygon (first == last)

        Returns None if exactly 2 IDL crossings are not found.

    Raises
    ------
    None
        Prints warning and returns None if polygon doesn't cross IDL exactly twice

    Notes
    -----
    1. **Crossing Detection**: The function identifies edges where longitude changes sign
       across the IDL. It checks both the polygon perimeter and the closing edge.

    2. **Intersection Calculation**: Uses great circle calculations via
       `find_great_circle_intersection_with_meridian` to find exact latitude where
       each edge crosses the 180° meridian.

    3. **Boundary Points**: Adds points slightly inside ±180° to avoid numerical issues:
       - Left polygon: 180° - 1e-10 (Eastern boundary)
       - Right polygon: -180° + 1e-10 (Western boundary)

    4. **Vertex Partitioning**: Vertices are assigned based on longitude sign:
       - Positive longitude → left polygon (Eastern hemisphere)
       - Negative longitude → right polygon (Western hemisphere)
       - Vertices at crossing points are included in both polygons

    5. **Polygon Closure**: Both resulting polygons are automatically closed by
       appending the first vertex to the end if not already closed.

    6. **CCW Order Requirement**: The algorithm assumes counter-clockwise vertex order.
       Clockwise polygons may produce incorrect results.

    7. **Two Crossings Only**: If the polygon crosses the IDL more than twice or
       doesn't cross at all, the function returns None with a warning.

    8. **Numerical Precision**: Uses 1e-10 offset from ±180° for boundary points and
       1e-10 tolerance for closure checks to handle floating-point precision.

    Warnings
    --------
    - Function prints warning and returns None if exactly 2 crossings not found
    - Does not validate CCW order - incorrect order will produce wrong results
    - Does not handle polygons with holes or multi-part geometries
    - May fail for self-intersecting or very complex polygons
    - Assumes continuous crossing (not multiple separate crossings)

    Examples
    --------
    >>> # Example 1: Simple rectangular polygon crossing IDL
    >>> coords = np.array([
    ...     [170.0, -10.0],
    ...     [170.0, 10.0],
    ...     [-170.0, 10.0],
    ...     [-170.0, -10.0],
    ...     [170.0, -10.0]
    ... ])
    >>> result = split_polygon_cross_idl(coords)
    >>> left, right = result
    >>> # Left polygon spans from 170° to ~180°
    >>> # Right polygon spans from ~-180° to -170°
    >>> np.all(left[:, 0] > 0)  # All Eastern
    True
    >>> np.all(right[:, 0] < 0)  # All Western
    True

    >>> # Example 2: Check polygon closure
    >>> left, right = split_polygon_cross_idl(coords)
    >>> np.allclose(left[0], left[-1])  # First == last
    True
    >>> np.allclose(right[0], right[-1])  # First == last
    True

    >>> # Example 3: Polygon that doesn't cross IDL
    >>> coords_no_cross = np.array([
    ...     [10.0, 0.0],
    ...     [20.0, 0.0],
    ...     [20.0, 10.0],
    ...     [10.0, 10.0],
    ...     [10.0, 0.0]
    ... ])
    >>> result = split_polygon_cross_idl(coords_no_cross)
    Warning: no intersection found
    >>> result is None
    True

    >>> # Example 4: Complex polygon with multiple vertices
    >>> coords_complex = np.array([
    ...     [175.0, 0.0],
    ...     [178.0, 2.0],
    ...     [179.0, 5.0],
    ...     [-179.0, 5.0],
    ...     [-178.0, 2.0],
    ...     [-175.0, 0.0],
    ...     [-178.0, -2.0],
    ...     [-179.0, -5.0],
    ...     [179.0, -5.0],
    ...     [178.0, -2.0],
    ...     [175.0, 0.0]
    ... ])
    >>> left, right = split_polygon_cross_idl(coords_complex)
    >>> # Both polygons contain vertices from original plus boundary points

    See Also
    --------
    reorder_idl_polygon : Alternative approach to handling IDL polygons
    convert_idl_polygon_to_valid_polygon : Longitude shifting approach
    find_great_circle_intersection_with_meridian : Great circle intersection calculation
    """
    # Convert input to numpy array for efficient operations
    aCoord_gcs = np.array(aCoord_gcs)

    # Get number of points in the polygon
    nPoint = len(aCoord_gcs)

    # Find indices where edges cross the International Date Line
    # The algorithm requires coordinates to be in counter-clockwise (CCW) order
    aIndex = []

    # Iterate through all edges (excluding the closing edge)
    # An edge crosses the IDL when longitude changes sign (+ to - or - to +)
    for i in range(nPoint - 1):
        dLongitude = aCoord_gcs[i,0]
        dLongitude_next = aCoord_gcs[i + 1,0]

        # Check for eastward crossing: positive to negative longitude (0° → 180° → -180°)
        if dLongitude > 0 and dLongitude < 180.0 and dLongitude_next < 0:
            aIndex.append(i)
            continue

        # Check for westward crossing: negative to positive longitude (-180° → 180° → 0°)
        if dLongitude < 0 and dLongitude_next > 0:
            aIndex.append(i)
            continue

    # Check the closing edge (from last point back to first point)
    # This ensures closed polygons are handled correctly
    dLongitude = aCoord_gcs[nPoint - 1,0]
    dLongitude_next = aCoord_gcs[0,0]
    if dLongitude > 0 and  dLongitude < 180.0 and dLongitude_next < 0:
        aIndex.append(nPoint - 1)
    if dLongitude < 0 and dLongitude_next > 0:
        aIndex.append(nPoint - 1)

    # Verify we found exactly 2 crossings
    # A polygon crossing the IDL should have exactly 2 edges that cross it
    if len(aIndex) != 2:
        print('Warning: no intersection found')
        print(aCoord_gcs)
        return

    # Calculate the first great circle intersection point with the IDL (180° meridian)
    # Get the edge endpoints for the first crossing
    lon1= aCoord_gcs[aIndex[0],0]
    lat1= aCoord_gcs[aIndex[0],1]
    lon2= aCoord_gcs[aIndex[0] + 1,0]
    lat2= aCoord_gcs[aIndex[0] + 1,1]
    target_lon = 180.0

    # Find where this great circle arc intersects the 180° meridian
    d, dLat0 = find_great_circle_intersection_with_meridian(lon1, lat1, lon2, lat2, target_lon)

    # Calculate the second great circle intersection point with the IDL
    # Get the edge endpoints for the second crossing
    lon1= aCoord_gcs[aIndex[1],0]
    lat1= aCoord_gcs[aIndex[1],1]

    # Handle the closing edge case (last point connects to first point)
    if aIndex[1] == nPoint - 1:
        lon2= aCoord_gcs[0,0]
        lat2= aCoord_gcs[0,1]
    else:
        lon2= aCoord_gcs[aIndex[1] + 1,0]
        lat2= aCoord_gcs[aIndex[1] + 1,1]
    target_lon = 180.0

    # Find where this great circle arc intersects the 180° meridian
    d, dLat1 = find_great_circle_intersection_with_meridian(lon1, lat1, lon2, lat2, target_lon)

    # Determine which intersection is north (top) and which is south (bottom)
    # This ordering is crucial for correctly partitioning the polygon vertices
    if dLat0 > dLat1:
        dLat_top = dLat0
        dLat_bottom = dLat1
    else:
        dLat_top = dLat1
        dLat_bottom = dLat0

    # Initialize lists to hold the two resulting polygons
    # Left polygon: Contains vertices with positive longitude (Eastern hemisphere side)
    # Right polygon: Contains vertices with negative longitude (Western hemisphere side)
    aCoord_gcs_left = list()
    aCoord_gcs_right = list()

    # Partition vertices based on longitude sign
    # Strategy: Traverse all vertices and assign to left or right based on longitude
    # At IDL crossings, add boundary points at exactly ±180° to close the polygons
    aIndex_dummy = np.array(aIndex)
    iFlag_added_right = 0  # Track if boundary points added to right polygon
    iFlag_added_left = 0   # Track if boundary points added to left polygon

    for i in range(nPoint):
        dLongitude = aCoord_gcs[i,0]

        # Process vertices with positive longitude (Eastern hemisphere)
        if dLongitude > 0:
            # Before the first crossing: add to left polygon
            if i <= np.min(aIndex_dummy):
                aCoord_gcs_left.append(aCoord_gcs[i])
                # At the crossing edge, add boundary points to close the polygon
                if i in aIndex:
                    if iFlag_added_left == 0:
                        iFlag_added_left = 1
                        # Add two points along the 180° meridian (slightly offset to avoid exactly 180°)
                        # These connect the split polygon edges: from bottom crossing to top crossing
                        aCoord_gcs_left.append([180-1.0E-10, dLat_bottom])
                        aCoord_gcs_left.append([180-1.0E-10, dLat_top])
            # After the first crossing: check if this is a crossing edge
            else:
                if i in aIndex:
                    aCoord_gcs_left.append(aCoord_gcs[i])
                    # At the second crossing, add boundary points if not already added
                    if iFlag_added_left == 0:
                        iFlag_added_left = 1
                        aCoord_gcs_left.append([180-1.0E-10, dLat_bottom])
                        aCoord_gcs_left.append([180-1.0E-10, dLat_top])
                else:
                    # Regular vertex with positive longitude: add to left polygon
                    aCoord_gcs_left.append(aCoord_gcs[i])

        # Process vertices with negative longitude (Western hemisphere)
        if dLongitude < 0:
            # Before the first crossing: add to right polygon
            if i <= np.min(aIndex_dummy):
                aCoord_gcs_right.append(aCoord_gcs[i])
                # At the crossing edge, add boundary points to close the polygon
                if i in aIndex:
                    if iFlag_added_right == 0:
                        iFlag_added_right = 1
                        # Add two points along the -180° meridian (slightly offset)
                        # Note: order is reversed (top to bottom) compared to left polygon
                        aCoord_gcs_right.append([-180+1.0E-10, dLat_top])
                        aCoord_gcs_right.append([-180+1.0E-10, dLat_bottom])
                    else:
                        pass
            # After the first crossing: check if this is a crossing edge
            else:
                if i in aIndex:
                    aCoord_gcs_right.append(aCoord_gcs[i])
                    # At the second crossing, add boundary points if not already added
                    if iFlag_added_right == 0:
                        aCoord_gcs_right.append([-180+1.0E-10, dLat_top])
                        aCoord_gcs_right.append([-180+1.0E-10, dLat_bottom])
                else:
                    # Regular vertex with negative longitude: add to right polygon
                    aCoord_gcs_right.append(aCoord_gcs[i])

    # Convert lists to numpy arrays for consistent output format
    aCoord_gcs_left = np.array(aCoord_gcs_left)
    aCoord_gcs_right = np.array(aCoord_gcs_right)

    # Ensure both polygons are properly closed (first point == last point)
    # This is critical for valid polygon geometry in GIS operations
    if len(aCoord_gcs_left) > 0:
        # Check if first and last points are the same (within tolerance)
        if not np.allclose(aCoord_gcs_left[0], aCoord_gcs_left[-1], atol=1e-10):
            # Close the polygon by appending the first point
            aCoord_gcs_left = np.vstack([aCoord_gcs_left, aCoord_gcs_left[0:1]])

    if len(aCoord_gcs_right) > 0:
        # Check if first and last points are the same (within tolerance)
        if not np.allclose(aCoord_gcs_right[0], aCoord_gcs_right[-1], atol=1e-10):
            # Close the polygon by appending the first point
            aCoord_gcs_right = np.vstack([aCoord_gcs_right, aCoord_gcs_right[0:1]])

    # Return both split polygons as a list
    # First element: Left polygon (positive longitudes, Eastern hemisphere)
    # Second element: Right polygon (negative longitudes, Western hemisphere)
    return [aCoord_gcs_left, aCoord_gcs_right]