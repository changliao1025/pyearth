"""
Douglas-Peucker Polyline Simplification for Geographic Coordinates
====================================================================

This module provides geodesic-aware implementation of the Douglas-Peucker algorithm
for simplifying polylines and polygons in geographic coordinate systems (longitude/latitude).

The Douglas-Peucker algorithm is a classical line simplification algorithm that reduces
the number of points in a curve while preserving its essential shape. This implementation
uses geodesic distance calculations to properly handle coordinates on a spherical Earth.

Key Features
------------
- Geodesic distance calculations for accurate simplification on sphere
- Automatic polygon detection and handling (closed rings)
- Recursive divide-and-conquer algorithm
- Preserves topology while reducing vertex count

Main Functions
--------------
douglas_peucker_geodetic : Simplify polyline/polygon using Douglas-Peucker algorithm

Algorithm Overview
------------------
The Douglas-Peucker algorithm works by:

1. Drawing a line between the first and last point of the polyline
2. Finding the point with maximum perpendicular distance from this line
3. If this distance exceeds the tolerance:
   - Recursively simplify the line segment before this point
   - Recursively simplify the line segment after this point
4. If distance is within tolerance:
   - Keep only the first and last point, discard all intermediate points

For geographic coordinates, perpendicular distance is calculated using great circle
(geodesic) distance on the sphere, not Euclidean distance.

Use Cases
---------
- Reducing file size of GPS tracks and routes
- Simplifying coastlines and political boundaries for map display
- Generalizing geographic features at different scales
- Reducing computational complexity for spatial operations

See Also
--------
pyearth.gis.geometry.calculate_distance_to_line : Calculate geodesic distance from point to line

References
----------
.. [1] Douglas, D. H., & Peucker, T. K. (1973). "Algorithms for the reduction of the
       number of points required to represent a digitized line or its caricature".
       Cartographica: The International Journal for Geographic Information and
       Geovisualization, 10(2), 112-122.

Examples
--------
>>> # Simplify a polyline
>>> polyline = [(0.0, 0.0), (0.1, 0.05), (0.2, 0.0), (0.3, -0.05), (0.4, 0.0)]
>>> simplified = douglas_peucker_geodetic(polyline, tolerance=1000.0)  # 1 km tolerance
>>> len(simplified) < len(polyline)
True

>>> # Simplify a polygon (automatically detected from closed ring)
>>> polygon = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 1.05), (0.0, 1.0), (0.0, 0.0)]
>>> simplified_polygon = douglas_peucker_geodetic(polygon, tolerance=5000.0)
>>> simplified_polygon[0] == simplified_polygon[-1]  # Still closed
True
"""

from typing import List, Tuple
from pyearth.gis.geometry.calculate_distance_to_line import calculate_distance_to_line


def douglas_peucker_geodetic(aCoords_gcs, tolerance):
    """
    Simplify a polyline or polygon using the Douglas-Peucker algorithm with geodesic distance.

    This function reduces the number of vertices in a polyline or polygon while preserving
    its essential shape. It uses geodesic (great circle) distance calculations to properly
    handle geographic coordinates on a spherical Earth.

    The algorithm recursively subdivides the line, keeping only vertices that are more than
    the specified tolerance distance away from the line connecting their neighbors.

    Parameters
    ----------
    aCoords_gcs : array-like of tuples or lists
        Input coordinates as a sequence of (longitude, latitude) pairs in decimal degrees.
        For polygons, the first and last coordinate should be identical (closed ring).
        Minimum 2 points required (3 for polygons including duplicate closure point).

        Format: [(lon1, lat1), (lon2, lat2), ..., (lonN, latN)]

        Longitude range: typically [-180, 180] or [0, 360] degrees
        Latitude range: [-90, 90] degrees

    tolerance : float
        Maximum perpendicular distance (in meters) allowed for a vertex to be removed.
        Vertices with perpendicular distance to the simplified line greater than this
        value will be kept. Smaller tolerance = more vertices retained = less simplification.

        Must be non-negative. Common values:
        - 10-100 m: High detail (minimal simplification)
        - 100-1000 m: Medium detail (moderate simplification)
        - 1000-10000 m: Low detail (aggressive simplification)

    Returns
    -------
    list of tuples
        Simplified coordinates as a list of (longitude, latitude) tuples.

        Properties:
        - Always includes first and last points of input
        - Number of points: 2 <= len(output) <= len(input)
        - For polygons: first and last point are identical (closed ring preserved)
        - Points are in same order as input
        - Coordinate format matches input format

    Raises
    ------
    ValueError
        If aCoords_gcs has fewer than 2 points
        If tolerance is negative

    Notes
    -----
    1. **Polygon Detection**: If the first and last coordinates are identical, the input
       is treated as a closed polygon. The closure point is temporarily removed during
       processing and restored in the output.

    2. **Distance Calculation**: Uses geodesic distance via `calculate_distance_to_line`,
       which accounts for Earth's spherical shape. This is more accurate than Euclidean
       distance, especially for:
       - Long line segments (> 100 km)
       - High latitudes (near poles)
       - Lines crossing the International Date Line

    3. **Algorithm Complexity**:
       - Best case: O(n) when all points are removed
       - Average case: O(n log n)
       - Worst case: O(n²) when all points are kept
       where n = number of input vertices

    4. **Tolerance Selection**: The optimal tolerance depends on:
       - Map scale: larger scale (zoomed in) → smaller tolerance
       - Data quality: noisy GPS tracks → larger tolerance
       - Performance: larger tolerance → fewer points → faster operations
       - Visual quality: balance between simplification and shape preservation

    5. **Topology Preservation**: The algorithm preserves:
       - Start and end points
       - Overall shape and direction
       - Closed ring property for polygons
       But may NOT preserve:
       - Self-intersections (can introduce new ones or remove existing ones)
       - Minimum bounding rectangle
       - Exact area or perimeter

    6. **Coordinate System**: Input coordinates should be in geographic coordinate system
       (GCS) with longitude/latitude in decimal degrees. Projected coordinates (meters)
       should be converted to GCS first.

    7. **Degenerate Cases**:
       - If all points are within tolerance: returns [first_point, last_point]
       - If input has 2 points: returns input unchanged
       - If input is a single point repeated: returns first and last point

    Warnings
    --------
    - Large tolerance values can drastically change polygon shape and area
    - For very long polylines (> 1000 points), consider iterative simplification
    - Geodesic distance calculation assumes WGS84 ellipsoid
    - Does not handle coordinates crossing International Date Line specially

    Examples
    --------
    >>> # Example 1: Simplify a simple polyline
    >>> coords = [(0.0, 0.0), (0.1, 0.05), (0.2, 0.0)]
    >>> simplified = douglas_peucker_geodetic(coords, tolerance=1000.0)
    >>> # Middle point removed if deviation < 1 km

    >>> # Example 2: Simplify a GPS track with noise
    >>> gps_track = [
    ...     (-122.4194, 37.7749),  # San Francisco
    ...     (-122.4195, 37.7750),  # Slight deviation
    ...     (-122.4196, 37.7751),  # Slight deviation
    ...     (-122.4197, 37.7752),  # Slight deviation
    ...     (-122.4000, 37.8000)   # Oakland
    ... ]
    >>> clean_track = douglas_peucker_geodetic(gps_track, tolerance=50.0)
    >>> # Removes GPS noise while keeping overall path

    >>> # Example 3: Simplify a polygon (closed ring)
    >>> polygon = [
    ...     (0.0, 0.0),
    ...     (0.5, 0.01),  # Nearly collinear
    ...     (1.0, 0.0),
    ...     (1.0, 1.0),
    ...     (0.0, 1.0),
    ...     (0.0, 0.0)   # Closure point
    ... ]
    >>> simplified_poly = douglas_peucker_geodetic(polygon, tolerance=1000.0)
    >>> simplified_poly[0] == simplified_poly[-1]  # Still closed
    True
    >>> # Nearly collinear point at (0.5, 0.01) likely removed

    >>> # Example 4: Coastline simplification at different scales
    >>> coastline = [...]  # Complex coastline with many vertices
    >>> # High detail for zoomed-in view
    >>> detail_high = douglas_peucker_geodetic(coastline, tolerance=100.0)
    >>> # Medium detail for regional view
    >>> detail_medium = douglas_peucker_geodetic(coastline, tolerance=1000.0)
    >>> # Low detail for continental view
    >>> detail_low = douglas_peucker_geodetic(coastline, tolerance=10000.0)
    >>> len(detail_high) > len(detail_medium) > len(detail_low)
    True

    >>> # Example 5: Check simplification ratio
    >>> original_count = len(coastline)
    >>> simplified = douglas_peucker_geodetic(coastline, tolerance=5000.0)
    >>> simplified_count = len(simplified)
    >>> ratio = simplified_count / original_count
    >>> print(f"Reduced to {ratio:.1%} of original vertices")

    See Also
    --------
    calculate_distance_to_line : Calculate geodesic perpendicular distance
    osgeo.ogr.Geometry.Simplify : GDAL's simplification (uses different algorithm)
    """
    # Input validation
    if not hasattr(aCoords_gcs, "__len__"):
        raise ValueError("aCoords_gcs must be a sequence of coordinate pairs")

    if len(aCoords_gcs) < 2:
        raise ValueError(
            f"aCoords_gcs must contain at least 2 points, got {len(aCoords_gcs)}"
        )

    if tolerance < 0:
        raise ValueError(f"tolerance must be non-negative, got {tolerance}")

    # Handle empty or single-point edge cases
    if len(aCoords_gcs) == 2:
        return [(coord[0], coord[1]) for coord in aCoords_gcs]
    # Handle empty or single-point edge cases
    if len(aCoords_gcs) == 2:
        return [(coord[0], coord[1]) for coord in aCoords_gcs]

    # Convert input coordinates to list of tuples
    aPoint = list()
    for aCoord in aCoords_gcs:
        aPoint.append((aCoord[0], aCoord[1]))

    # Detect if input is a closed polygon
    is_polygon = aPoint[0] == aPoint[-1]
    if is_polygon:
        # Remove closure point temporarily; will be restored at end
        aPoint = aPoint[:-1]

    # Nested helper function: Find point with maximum distance from line segment
    def find_furthest_point(points, start, end):
        """
        Find the point in points[start+1:end] with maximum distance from line points[start]-points[end].

        Parameters
        ----------
        points : list of tuples
            All points in the polyline
        start : int
            Index of line segment start point
        end : int
            Index of line segment end point

        Returns
        -------
        index : int
            Index of the furthest point between start and end
        max_distance : float
            Maximum perpendicular distance in meters
        """
        max_distance = 0
        index = 0
        for i in range(start + 1, end):
            distance = point_line_distance(points[i], points[start], points[end])
            if distance > max_distance:
                max_distance = distance
                index = i
        return index, max_distance

    # Nested helper function: Calculate geodesic perpendicular distance
    def point_line_distance(point, start, end):
        """
        Calculate perpendicular distance from point to line segment using geodesic distance.

        Parameters
        ----------
        point : tuple
            (longitude, latitude) of the point to measure
        start : tuple
            (longitude, latitude) of line segment start
        end : tuple
            (longitude, latitude) of line segment end

        Returns
        -------
        float
            Perpendicular distance in meters from point to line
        """
        # Extract coordinates
        dLongitude0 = point[0]  # Point to measure (middle point)
        dLatitude0 = point[1]
        dLongitude1 = start[0]  # Line segment start
        dLatitude1 = start[1]
        dLongitude2 = end[0]  # Line segment end
        dLatitude2 = end[1]

        # Calculate geodesic perpendicular distance
        # Note: calculate_distance_to_line handles the geodesic calculation
        distance = calculate_distance_to_line(
            dLongitude1, dLatitude1, dLongitude0, dLatitude0, dLongitude2, dLatitude2
        )

        return distance

    # Nested recursive function: Simplify polyline segment
    def simplify(points, start, end, tolerance, simplified):
        """
        Recursively simplify a segment of the polyline using Douglas-Peucker algorithm.

        Parameters
        ----------
        points : list of tuples
            All points in the polyline
        start : int
            Index of segment start point
        end : int
            Index of segment end point
        tolerance : float
            Distance tolerance in meters
        simplified : list of int
            Accumulated list of indices to keep (modified in place)
        """
        # Find point with maximum perpendicular distance
        index, max_distance = find_furthest_point(points, start, end)

        if max_distance > tolerance:
            # Point exceeds tolerance - keep it and recursively process sub-segments
            simplify(points, start, index, tolerance, simplified)
            simplify(points, index, end, tolerance, simplified)
        else:
            # All points within tolerance - keep only start and end
            if start not in simplified:
                simplified.append(start)
            if end not in simplified:
                simplified.append(end)

    # Initialize and run simplification
    simplified = []
    simplify(aPoint, 0, len(aPoint) - 1, tolerance, simplified)

    # Sort indices to maintain original point order
    simplified.sort()

    # Extract simplified points using kept indices
    simplified_points = [aPoint[i] for i in simplified]

    # Restore polygon closure if input was a closed ring
    if is_polygon:
        simplified_points.append(simplified_points[0])

    return simplified_points
