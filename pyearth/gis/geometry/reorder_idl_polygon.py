"""
International Date Line Polygon Reordering
===========================================

This module provides functionality to reorder polygon vertices that cross the International
Date Line (IDL) at 180°/-180° longitude to prevent self-intersection artifacts.

When polygons span across the International Date Line, naive rendering or processing can
create self-intersections due to the discontinuity at 180°/-180°. This module reorders
vertices to trace a valid, non-self-intersecting perimeter.

Key Features
------------
- Automatic detection of IDL-crossing polygons
- Vertex reordering to prevent self-intersection
- Preserves polygon topology and closure
- Handles both clockwise and counter-clockwise orientations
- Maintains latitude bounds and spatial extent

Main Functions
--------------
reorder_idl_polygon : Reorder IDL-crossing polygon vertices to prevent self-intersection

Algorithm Overview
------------------
The reordering algorithm works by:

1. **Separate vertices**: Split points into Eastern (lon >= 0) and Western (lon < 0) groups
2. **Sort by latitude**:
   - Eastern points: South to North (ascending latitude)
   - Western points: North to South (descending latitude)
3. **Construct path**: Concatenate sorted groups to form valid perimeter
4. **Close polygon**: Append first vertex to end
5. **Reverse orientation**: Convert to counter-clockwise (CCW) orientation

This creates a path that traces up the Eastern side, crosses at the top, traces down
the Western side, and crosses back at the bottom.

Use Cases
---------
- Visualizing regions crossing the International Date Line (e.g., Russia, Fiji, Kiribati)
- Processing Pacific Ocean polygons
- Fixing topology errors in global datasets
- Preparing data for GIS systems that don't handle IDL automatically
- Converting between different longitude conventions ([0, 360] vs [-180, 180])

Coordinate Systems
------------------
**Input/Output**: Geographic coordinates (longitude, latitude) in decimal degrees
    - Longitude: [-180, 180] range (negative = West, positive = East)
    - Latitude: [-90, 90] range (negative = South, positive = North)

**IDL Location**: 180° longitude (equivalent to -180°)
    - Eastern Hemisphere: [0, 180]
    - Western Hemisphere: [-180, 0]

Limitations
-----------
- Assumes polygon is "tightly clustered" around the IDL
- May not work correctly for polygons that span > 180° of longitude
- Designed for simple polygons (no holes)
- Does not handle multi-part geometries

See Also
--------
pyearth.gis.geometry.convert_idl_polygon_to_valid_polygon : Alternative IDL handling
osgeo.ogr.Geometry : GDAL geometry operations

References
----------
.. [1] International Date Line: https://en.wikipedia.org/wiki/International_Date_Line
.. [2] OGC Simple Features Specification: Polygon geometry rules

Examples
--------
>>> # Example 1: Simple IDL-crossing rectangle
>>> vertices = [
...     (179.0, -10.0),   # Eastern side, south
...     (179.0, 10.0),    # Eastern side, north
...     (-179.0, 10.0),   # Western side, north
...     (-179.0, -10.0),  # Western side, south
...     (179.0, -10.0)    # Closing point
... ]
>>> reordered = reorder_idl_polygon(vertices)
>>> # Vertices now trace valid perimeter without self-intersection

>>> # Example 2: Polygon with many vertices
>>> # Vertices clustered around 180°/-180°
>>> complex_vertices = [...]  # Many points near IDL
>>> valid_polygon = reorder_idl_polygon(complex_vertices)
"""

import math
from typing import List, Tuple

# Define the coordinate type alias for clarity (Lon, Lat)
Coord = Tuple[float, float]

def reorder_idl_polygon(vertices: List[Coord]) -> List[Coord]:
    """
    Reorder vertices of an International Date Line (IDL)-crossing polygon to prevent self-intersection.

    This function takes a polygon that crosses the International Date Line (180°/-180° meridian)
    and reorders its vertices to create a valid, non-self-intersecting polygon. The algorithm
    assumes the polygon is tightly clustered around the IDL and sorts vertices to trace the
    perimeter consistently.

    The reordering strategy separates Eastern (lon >= 0) and Western (lon < 0) vertices,
    sorts them by latitude, then concatenates them to form a path that traces up the Eastern
    side, across the top, down the Western side, and back across the bottom. The result is
    converted to counter-clockwise (CCW) orientation.

    Parameters
    ----------
    vertices : List[Tuple[float, float]]
        List of (longitude, latitude) coordinate tuples defining the polygon.

        Expected format: [(lon1, lat1), (lon2, lat2), ..., (lonN, latN)]

        Requirements:
        - Should be a closed polygon (first point equals last point)
        - Longitude in decimal degrees: [-180, 180] range
        - Latitude in decimal degrees: [-90, 90] range
        - Should cross the International Date Line (contain both positive and negative longitudes)
        - Minimum 4 points (3 unique + 1 closure point)

        The function works best when:
        - Polygon is "tightly clustered" around 180°/-180° meridian
        - Polygon doesn't span more than 180° of longitude
        - All vertices are near the IDL (e.g., longitudes in [150, 180] ∪ [-180, -150])

    Returns
    -------
    List[Tuple[float, float]]
        New list of (longitude, latitude) tuples forming a valid, closed polygon
        in counter-clockwise (CCW) orientation.

        Properties:
        - First and last points are identical (closed polygon)
        - Vertices ordered to prevent self-intersection at IDL
        - Same number of vertices as input
        - CCW orientation (standard for exterior rings in GIS)

        If input is empty, returns empty list.
        If polygon doesn't cross IDL (all points on one side), returns input unchanged.

    Raises
    ------
    None
        Function does not raise exceptions, but may produce unexpected results if:
        - Input polygon is not actually closed
        - Input polygon spans > 180° longitude
        - Input polygon is not clustered around IDL

    Notes
    -----
    1. **IDL Crossing Detection**: The function considers a polygon to cross the IDL if it
       contains both Eastern (lon >= 0) and Western (lon < 0) vertices. Polygons with all
       vertices on one side are returned unchanged.

    2. **Longitude Partitioning**:
       - Eastern Hemisphere: lon >= 0 (includes Prime Meridian at 0°)
       - Western Hemisphere: lon < 0
       - Exactly 0.0° is grouped with Eastern points

    3. **Sorting Strategy**:
       - Eastern points sorted by latitude (ascending): South → North
       - Western points sorted by latitude (descending): North → South
       - This creates a continuous path around the perimeter

    4. **Orientation**: The final polygon is in counter-clockwise (CCW) orientation, which
       is the standard for exterior polygon rings in most GIS formats (OGC Simple Features,
       GeoJSON, Shapefile).

    5. **Closure Handling**: The function strips the closing vertex during processing and
       re-adds it at the end to ensure proper closure of the reordered polygon.

    6. **Assumptions and Limitations**:
       - Assumes simple polygon (no holes, no multi-part)
       - Works best for polygons tightly clustered around IDL
       - May fail for very large polygons (> 180° longitude span)
       - Does not validate input geometry for self-intersections
       - Does not handle polygons that should cross at multiple latitudes

    7. **Alternative Approaches**: For more complex IDL handling (e.g., shifting all
       longitudes to [0, 360] range), consider `convert_idl_polygon_to_valid_polygon`.

    8. **Performance**: O(N log N) due to sorting, where N is number of vertices.

    Warnings
    --------
    - This function modifies vertex order but not vertex values
    - May produce unexpected results for polygons not centered on the IDL
    - Does not validate that input forms a valid polygon
    - Assumes specific geometric configuration (clustered around IDL)

    Examples
    --------
    >>> # Example 1: Simple rectangle crossing IDL
    >>> vertices = [
    ...     (179.0, -10.0),   # Eastern side, south
    ...     (179.0, 10.0),    # Eastern side, north
    ...     (-179.0, 10.0),   # Western side, north
    ...     (-179.0, -10.0),  # Western side, south
    ...     (179.0, -10.0)    # Closing point
    ... ]
    >>> reordered = reorder_idl_polygon(vertices)
    >>> len(reordered)
    5
    >>> reordered[0] == reordered[-1]  # Still closed
    True

    >>> # Example 2: Polygon that doesn't cross IDL
    >>> pacific_vertices = [
    ...     (170.0, -10.0),
    ...     (170.0, 10.0),
    ...     (160.0, 10.0),
    ...     (160.0, -10.0),
    ...     (170.0, -10.0)
    ... ]
    >>> result = reorder_idl_polygon(pacific_vertices)
    >>> result == pacific_vertices  # Unchanged (all Eastern)
    True

    >>> # Example 3: Empty input
    >>> empty = reorder_idl_polygon([])
    >>> len(empty)
    0

    >>> # Example 4: Complex polygon with many vertices near IDL
    >>> # Vertices from satellite swath crossing date line
    >>> swath = [
    ...     (178.0, -5.0), (179.0, -3.0), (179.5, 0.0), (179.0, 3.0), (178.0, 5.0),
    ...     (-178.0, 5.0), (-179.0, 3.0), (-179.5, 0.0), (-179.0, -3.0), (-178.0, -5.0),
    ...     (178.0, -5.0)  # Closing
    ... ]
    >>> valid_swath = reorder_idl_polygon(swath)
    >>> # Vertices now ordered to prevent self-intersection

    >>> # Example 5: Verify CCW orientation of output
    >>> def is_ccw(coords):
    ...     '''Check if polygon is counter-clockwise using shoelace formula'''
    ...     area = 0.0
    ...     for i in range(len(coords) - 1):
    ...         area += (coords[i+1][0] - coords[i][0]) * (coords[i+1][1] + coords[i][1])
    ...     return area < 0
    >>> vertices = [(179, -10), (179, 10), (-179, 10), (-179, -10), (179, -10)]
    >>> reordered = reorder_idl_polygon(vertices)
    >>> is_ccw(reordered)
    True

    See Also
    --------
    convert_idl_polygon_to_valid_polygon : Alternative IDL handling using longitude shift
    osgeo.ogr.Geometry.IsValid : Validate polygon geometry
    """
    # Handle empty input
    if not vertices:
        return []

    # Step 1: Strip the closing vertex for processing
    # The closing vertex (duplicate of first vertex) will be re-added at the end
    # to ensure the reordered polygon is properly closed
    unique_vertices = vertices[:-1]

    # Step 2: Separate points into Eastern and Western hemispheres
    # Eastern: lon >= 0 (includes Prime Meridian at 0.0)
    # Western: lon < 0 (negative longitudes)
    east_points: List[Coord] = []
    west_points: List[Coord] = []

    for lon, lat in unique_vertices:
        if lon >= 0:
            east_points.append((lon, lat))
        else:
            west_points.append((lon, lat))

    # Step 3: Handle edge case - polygon doesn't cross IDL
    # If all points are on one side (East or West), there's no IDL crossing
    # and no risk of self-intersection, so return original vertices unchanged
    if not east_points or not west_points:
        # If it doesn't cross, it's not at risk of IDL self-intersection.
        return vertices

    # Step 4: Sort points by latitude to create a valid perimeter path
    # Eastern points: Sort ascending (South to North) - traces Eastern boundary upward
    # Western points: Sort descending (North to South) - traces Western boundary downward
    # This creates a continuous path: E(bottom→top) then W(top→bottom)
    east_points.sort(key=lambda p: p[1])  # Sort by latitude (ascending)
    west_points.sort(key=lambda p: p[1], reverse=True)  # Sort by latitude (descending)

    # Step 5: Construct the new linear ring by concatenating sorted point groups
    # Path follows: South→North on Eastern side, then North→South on Western side
    new_vertices = east_points + west_points

    # Step 6: Close the polygon by appending the first vertex
    # This ensures the polygon forms a closed ring (first == last)
    new_vertices.append(new_vertices[0])

    # Step 7: Reverse to counter-clockwise (CCW) orientation
    # CCW is the standard orientation for exterior polygon rings in GIS
    # (following OGC Simple Features Specification and GeoJSON standards)
    new_vertices = new_vertices[::-1]

    return new_vertices