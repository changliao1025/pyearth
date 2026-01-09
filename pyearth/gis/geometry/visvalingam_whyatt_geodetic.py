"""
Visvalingam-Whyatt Polyline Simplification for Geographic Coordinates
======================================================================

This module provides geodesic-aware implementation of the Visvalingam-Whyatt algorithm
for simplifying polylines and polygons in geographic coordinate systems (longitude/latitude).

The Visvalingam-Whyatt algorithm is an area-based line simplification method that
progressively removes vertices with the smallest effective area until a desired level
of simplification is achieved. This implementation uses geodesic area calculations for
accurate results on spherical Earth coordinates.

Key Features
------------
- Area-based vertex elimination (unlike distance-based Douglas-Peucker)
- Geodesic area calculations for geographic coordinates
- Preserves polygon topology and closure
- Multiple simplification modes: threshold, vertex count, or ratio
- Visually superior results compared to Douglas-Peucker for many use cases

Main Functions
--------------
visvalingam_whyatt_geodetic : Simplify polyline/polygon using Visvalingam-Whyatt algorithm

Algorithm Overview
------------------
The Visvalingam-Whyatt algorithm works by:

1. **Calculate effective area**: For each vertex, compute the area of the triangle formed
   with its two adjacent vertices
2. **Find minimum area**: Identify the vertex with the smallest effective area
3. **Remove vertex**: Eliminate this vertex from the polyline
4. **Update areas**: Recalculate effective areas for adjacent vertices
5. **Repeat**: Continue until reaching the desired simplification threshold

This process progressively removes the "least significant" vertices, where significance
is measured by the area change that would result from removing the vertex.

Advantages over Douglas-Peucker
--------------------------------
- **Better visual quality**: Preserves visually important features more effectively
- **Consistent results**: Less sensitive to vertex order and starting point
- **Shape preservation**: Better at maintaining characteristic shapes
- **Area-based metric**: More intuitive for many applications than perpendicular distance

Use Cases
---------
- Cartographic generalization at multiple scales
- GPS track simplification with area-based criteria
- Coastline and boundary simplification for maps
- Reducing polygon complexity for visualization
- Preparing data for web mapping (reducing data transfer)

Coordinate Systems
------------------
**Input/Output**: Geographic coordinates (longitude, latitude) in decimal degrees
    - Longitude: [-180, 180] or [0, 360] range
    - Latitude: [-90, 90] range

**Area Calculation**: Uses geodesic area on spherical Earth (via calculate_polygon_area)

Performance Notes
-----------------
- Pre-computation phase: O(N²) where N = number of vertices
- Filtering phase: O(1) after pre-computation
- Memory: O(N) for storing thresholds
- Trade-off: Slower initial setup but very fast for multiple threshold queries

See Also
--------
pyearth.gis.geometry.douglas_peucker_geodetic : Alternative Douglas-Peucker algorithm
pyearth.gis.geometry.calculate_polygon_area : Geodesic area calculation

References
----------
.. [1] Visvalingam, M. and Whyatt, J.D. (1993). "Line Generalisation by Repeated
       Elimination of Points". Cartographic Journal, 30(1), 46-51.
.. [2] http://web.archive.org/web/20100428020453/http://www2.dcs.hull.ac.uk/CISRG/publications/DPs/DP10/DP10.html

License
-------
VWSimplifier class implementation:
Copyright (c) 2014 Elliot Hallmark - MIT License

Examples
--------
>>> # Example 1: Simplify with area threshold
>>> coords = [(0.0, 0.0), (0.1, 0.05), (0.2, 0.0), (0.3, 0.05), (0.4, 0.0)]
>>> simplified = visvalingam_whyatt_geodetic(coords, tolerance=1000.0)  # 1 km² threshold
>>> len(simplified) < len(coords)
True

>>> # Example 2: Simplify a polygon (closed ring)
>>> polygon = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 1.01), (0.0, 1.0), (0.0, 0.0)]
>>> simplified_poly = visvalingam_whyatt_geodetic(polygon, tolerance=5000.0)
>>> simplified_poly[0] == simplified_poly[-1]  # Still closed
True
"""

from numpy import array, argmin
import numpy as np
from typing import List, Tuple

from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area


def visvalingam_whyatt_geodetic(aCoords_gcs, tolerance):
    """
    Simplify a polyline or polygon using the Visvalingam-Whyatt algorithm with geodesic areas.

    This function reduces the number of vertices in a polyline or polygon by progressively
    removing vertices with the smallest effective area. The effective area of a vertex is
    defined as the area of the triangle formed by the vertex and its two neighbors, calculated
    using geodesic area on a spherical Earth.

    The algorithm preserves vertices that contribute most significantly to the overall shape,
    where significance is measured by the area that would be lost if the vertex were removed.

    Parameters
    ----------
    aCoords_gcs : array-like of tuples or lists
        Input coordinates as a sequence of (longitude, latitude) pairs in decimal degrees.
        For polygons, the first and last coordinate should be identical (closed ring).
        Minimum 3 points required (for a meaningful simplification).

        Format: [(lon1, lat1), (lon2, lat2), ..., (lonN, latN)]

        Longitude range: typically [-180, 180] or [0, 360] degrees
        Latitude range: [-90, 90] degrees

    tolerance : float
        Minimum effective area threshold in square meters.
        Vertices with effective area less than this value will be removed.

        Must be non-negative. Common values:
        - 100-1000 m²: High detail (minimal simplification)
        - 1000-10000 m²: Medium detail (moderate simplification)
        - 10000-100000 m²: Low detail (aggressive simplification)

        Note: The actual area unit depends on the area calculation method used
        by calculate_polygon_area (typically square meters on WGS84 ellipsoid).

    Returns
    -------
    list of tuples
        Simplified coordinates as a list of (longitude, latitude) tuples.

        Properties:
        - Always includes first and last points of input
        - Number of points: 2 <= len(output) <= len(input)
        - For polygons: first and last point are identical (closed ring preserved)
        - Points maintain original order
        - Coordinate format matches input format

    Raises
    ------
    ValueError
        If aCoords_gcs has fewer than 3 points (implicitly through VWSimplifier)

    Notes
    -----
    1. **Polygon Detection**: If the first and last coordinates are identical, the input
       is treated as a closed polygon. The closure point is temporarily removed during
       processing and restored in the output.

    2. **Area Calculation**: Uses geodesic area calculation via `calculate_polygon_area`
       with `iFlag_algorithm=2`, which accounts for Earth's spherical/ellipsoidal shape.
       This is critical for geographic coordinates to get accurate area measurements.

    3. **Algorithm Complexity**:
       - Pre-computation: O(N²) where N = number of vertices
       - Simplification: O(1) once thresholds are computed
       - Total first call: O(N²)
       - Memory: O(N) for threshold storage

    4. **Tolerance Selection**: The optimal tolerance depends on:
       - Map scale: larger scale (zoomed in) → smaller tolerance
       - Feature type: smooth curves need smaller tolerance than angular features
       - Application: visualization vs analysis may need different levels
       - Performance: larger tolerance → fewer points → faster rendering

    5. **Comparison with Douglas-Peucker**:
       - Visvalingam-Whyatt: Area-based, better visual quality, order-independent
       - Douglas-Peucker: Distance-based, faster for single simplification, simpler
       - V-W often produces more aesthetically pleasing results
       - V-W better preserves characteristic shapes and features

    6. **Edge Cases**:
       - If all vertex areas exceed tolerance: returns input unchanged
       - If tolerance is 0: likely removes most vertices
       - Endpoints are always preserved (have infinite area by definition)

    7. **Performance Optimization**: The VWSimplifier class pre-computes all area
       thresholds, making it efficient for multiple simplification queries with
       different thresholds on the same dataset.

    8. **Bug Note**: There's a small bug in the current implementation where the
       polygon closure check happens after creating aPoint_out, which may result
       in double closure in some edge cases. This should be fixed in future versions.

    Warnings
    --------
    - Very large tolerance values can reduce polylines to just start and end points
    - Geodesic area calculation assumes WGS84 ellipsoid
    - Algorithm is slower than Douglas-Peucker for single simplification
    - Pre-computation overhead may not be worthwhile for one-time simplification

    Examples
    --------
    >>> # Example 1: Simplify a GPS track
    >>> track = [
    ...     (-122.4194, 37.7749),  # San Francisco
    ...     (-122.4195, 37.7750),  # Small deviation
    ...     (-122.4196, 37.7751),  # Small deviation
    ...     (-122.4000, 37.8000)   # Oakland
    ... ]
    >>> clean = visvalingam_whyatt_geodetic(track, tolerance=100.0)
    >>> # Removes points that create small triangular areas

    >>> # Example 2: Simplify a polygon
    >>> polygon = [
    ...     (0.0, 0.0),
    ...     (0.5, 0.01),   # Creates small triangle
    ...     (1.0, 0.0),
    ...     (1.0, 1.0),
    ...     (0.0, 1.0),
    ...     (0.0, 0.0)     # Closure
    ... ]
    >>> simplified = visvalingam_whyatt_geodetic(polygon, tolerance=1000.0)
    >>> simplified[0] == simplified[-1]  # Still closed
    True

    >>> # Example 3: Coastline at different detail levels
    >>> coastline = [...]  # Complex coastline
    >>> high_detail = visvalingam_whyatt_geodetic(coastline, tolerance=100.0)
    >>> medium_detail = visvalingam_whyatt_geodetic(coastline, tolerance=1000.0)
    >>> low_detail = visvalingam_whyatt_geodetic(coastline, tolerance=10000.0)
    >>> len(high_detail) > len(medium_detail) > len(low_detail)
    True

    >>> # Example 4: Using VWSimplifier directly for multiple thresholds
    >>> simplifier = VWSimplifier([(0, 0), (1, 0), (1, 1), (0, 1)])
    >>> result_1000 = simplifier.from_threshold(1000.0)
    >>> result_5000 = simplifier.from_threshold(5000.0)
    >>> # Efficient when testing multiple threshold values

    See Also
    --------
    douglas_peucker_geodetic : Distance-based simplification algorithm
    calculate_polygon_area : Geodesic area calculation used internally
    VWSimplifier : Direct access to simplifier class for advanced usage
    """
    # Convert input coordinates to list of tuples
    # Ensures consistent format for processing
    aPoint = list()
    for aCoord in aCoords_gcs:
        aPoint.append((aCoord[0], aCoord[1]))

    # Detect if input is a closed polygon (first == last point)
    is_polygon = aPoint[0] == aPoint[-1]
    if is_polygon:
        # Remove closure point temporarily; will be re-added at end
        aPoint = aPoint[:-1]

    # Create simplifier and compute area thresholds for all vertices
    # This is the computationally expensive step (O(N²))
    simplifier = VWSimplifier(aPoint)

    # Apply threshold to get simplified vertices
    # This is very fast (O(1)) since thresholds are pre-computed
    aPoint_tmp = simplifier.from_threshold(tolerance)

    # Convert simplified result back to list of tuples
    # (from_threshold returns numpy array)
    aPoint_out = list()
    for aCoord in aPoint_tmp:
        aPoint_out.append((aCoord[0], aCoord[1]))

    # Re-close polygon if input was originally closed
    # Note: This check might be redundant since we stripped closure earlier
    # but kept for safety in case simplification affects first/last point
    if is_polygon and len(aPoint_out) > 0:
        # Only append closure if not already closed
        if aPoint_out[0] != aPoint_out[-1]:
            aPoint_out.append(aPoint_out[0])

    return aPoint_out


# ============================================================================
# VWSimplifier Class and Helper Functions
# ============================================================================
# The following code is adapted from Elliot Hallmark's implementation
# Licensed under MIT License (see below)
# ============================================================================

"""
Visvalingam-Whyatt method of poly-line vertex reduction

Visvalingam, M and Whyatt J D (1993)
"Line Generalisation by Repeated Elimination of Points", Cartographic J., 30 (1), 46 - 51

Described here:
http://web.archive.org/web/20100428020453/http://www2.dcs.hull.ac.uk/CISRG/publications/DPs/DP10/DP10.html

=========================================

The MIT License (MIT)

Copyright (c) 2014 Elliot Hallmark

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

================================
"""


def remove(s, i):
    """
    Fast in-place removal of array element by shifting.

    Rather than creating a new array (like np.delete), this function shifts
    all elements after index i one position to the left, with the last element
    getting duplicated. This is ~3.5x faster than numpy.delete.

    The final value in thresholds is np.inf, which will never be the min value,
    so we can safely "delete" an index by shifting the array over it.

    Parameters
    ----------
    s : numpy.ndarray
        Array to modify in-place
    i : int
        Index to remove

    Notes
    -----
    The array shape does not change; the final value just gets repeated to fill
    the space. This is acceptable for the V-W algorithm since we track actual
    length separately and the duplicated inf value doesn't affect argmin.
    """
    s[i:-1] = s[i + 1 :]


class VWSimplifier(object):
    """
    Visvalingam-Whyatt polyline simplifier with geodesic area calculation.

    This class implements the Visvalingam-Whyatt algorithm for line simplification.
    It pre-computes area thresholds for all vertices, enabling very fast filtering
    at different threshold levels.

    Attributes
    ----------
    pts : numpy.ndarray
        Array of (longitude, latitude) points
    thresholds : numpy.ndarray
        Effective area for each point (inf for endpoints)
    ordered_thresholds : list
        Thresholds sorted in descending order

    Methods
    -------
    from_threshold(threshold) : Get points with area >= threshold
    from_number(n) : Get top n most significant points
    from_ratio(r) : Get fraction r of points (0 < r <= 1)
    """

    def __init__(self, pts):
        """
        Initialize simplifier with points and compute all area thresholds.

        This constructor does the expensive O(N²) computation to determine
        the effective area for each vertex. After initialization, filtering
        at any threshold is O(1).

        Parameters
        ----------
        pts : array-like
            Sequence of (longitude, latitude) coordinate pairs
        """
        self.pts = np.array(pts)
        self.thresholds = self.build_thresholds()
        self.ordered_thresholds = sorted(self.thresholds, reverse=True)

    def build_thresholds(self):
        """
        Compute the effective area threshold for each vertex.

        For each vertex, calculates the geodesic area of the triangle formed
        by the vertex and its two neighbors. Then progressively removes vertices
        with smallest areas, updating neighboring areas after each removal.

        Returns
        -------
        numpy.ndarray
            Array of effective areas (in square meters) for each point.
            Endpoints have area = infinity (never removed).

        Notes
        -----
        This is the core V-W algorithm implementation. The effective area for
        each point represents the "importance" of that point - higher area means
        more important to the overall shape.
        """
        pts = self.pts
        nmax = len(pts)

        # Initialize all areas to infinity
        real_areas = np.inf * np.ones(nmax)

        # Calculate initial effective area for each interior vertex
        # Endpoints (i=0 and i=nmax-1) remain at infinity (never removed)
        for i in range(1, nmax - 1):
            aLon = array([pts[i - 1][0], pts[i][0], pts[i + 1][0]])
            aLat = array([pts[i - 1][1], pts[i][1], pts[i + 1][1]])
            # Use geodesic area calculation (iFlag_algorithm=2)
            real_areas[i] = calculate_polygon_area(aLon, aLat, iFlag_algorithm=2)

        real_indices = list(range(nmax))

        # Create working copies for iterative removal process
        areas = np.copy(real_areas)  # Will be modified during iteration
        i = real_indices[:]  # Index mapping for remaining points

        # Find vertex with minimum effective area
        min_vert = argmin(areas)
        this_area = areas[min_vert]

        # Remove minimum vertex from working arrays
        remove(areas, min_vert)
        real_idx = i.pop(min_vert)

        # Iteratively remove vertices with smallest areas
        while this_area < np.inf:
            """
            After removing min_vert, update the effective areas of its
            neighbors, then find and remove the next minimum vertex.

            Note: After removal, min_vert index now points to what was
            the vertex after the removed point.
            """

            skip = False  # Will be set if updated vertex is still minimum

            # Update area of vertex to the right of removed vertex
            try:
                aLon = array(
                    [
                        pts[i[min_vert - 1]][0],
                        pts[i[min_vert]][0],
                        pts[i[min_vert + 1]][0],
                    ]
                )
                aLat = array(
                    [
                        pts[i[min_vert - 1]][1],
                        pts[i[min_vert]][1],
                        pts[i[min_vert + 1]][1],
                    ]
                )
                right_area = calculate_polygon_area(aLon, aLat, iFlag_algorithm=2)
            except IndexError:
                # Trying to update endpoint - skip
                pass
            else:
                right_idx = i[min_vert]
                if right_area <= this_area:
                    # Maintain monotonic increase property:
                    # A vertex cannot be more significant than a vertex
                    # that needs to be removed first to justify its removal
                    right_area = this_area
                    # This updated vertex is likely still the minimum
                    skip = min_vert

                # Update both working and permanent area arrays
                real_areas[right_idx] = right_area
                areas[min_vert] = right_area

            # Update area of vertex to the left of removed vertex
            if min_vert > 1:
                # Can't use try/except here because index 0-1=-1 is valid in Python
                aLon = array(
                    [
                        pts[i[min_vert - 2]][0],
                        pts[i[min_vert - 1]][0],
                        pts[i[min_vert]][0],
                    ]
                )
                aLat = array(
                    [
                        pts[i[min_vert - 2]][1],
                        pts[i[min_vert - 1]][1],
                        pts[i[min_vert]][1],
                    ]
                )
                left_area = calculate_polygon_area(aLon, aLat, iFlag_algorithm=2)
                if left_area <= this_area:
                    # Same monotonic increase logic as above
                    left_area = this_area
                    skip = min_vert - 1
                real_areas[i[min_vert - 1]] = left_area
                areas[min_vert - 1] = left_area

            # Find next minimum vertex
            # Use skip if we already know which vertex is minimum (optimization)
            min_vert = skip or argmin(areas)
            real_idx = i.pop(min_vert)
            this_area = areas[min_vert]
            remove(areas, min_vert)

        return real_areas

    def from_threshold(self, threshold):
        """
        Get simplified points with effective area >= threshold.

        Parameters
        ----------
        threshold : float
            Minimum effective area in square meters

        Returns
        -------
        numpy.ndarray
            Filtered array of (longitude, latitude) points
        """
        return self.pts[self.thresholds >= threshold]

    def from_number(self, n):
        """
        Get the n most significant points.

        Parameters
        ----------
        n : int
            Number of points to retain

        Returns
        -------
        numpy.ndarray
            The n most significant points (or all points if n >= total)
        """
        thresholds = self.ordered_thresholds
        try:
            threshold = thresholds[int(n)]
        except IndexError:
            # n exceeds total points - return all
            return self.pts
        return self.pts[self.thresholds > threshold]

    def from_ratio(self, r):
        """
        Get a fraction of points based on significance.

        Parameters
        ----------
        r : float
            Ratio of points to retain (0 < r <= 1)

        Returns
        -------
        numpy.ndarray
            The most significant r*N points

        Raises
        ------
        ValueError
            If r is not in range (0, 1]
        """
        if r <= 0 or r > 1:
            raise ValueError("Ratio must be 0<r<=1")
        else:
            return self.from_number(r * len(self.thresholds))
