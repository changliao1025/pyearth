import numpy as np
from osgeo import ogr
from typing import List, Tuple, Optional
from pyearth.gis.geometry.calculate_intersect_on_great_circle import (
    find_great_circle_intersection_with_meridian,
)
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates

Coord = Tuple[float, float]


def unwrap_longitudes(coords: np.ndarray) -> np.ndarray:
    """Unwrap longitude coordinates to handle International Date Line crossings.

    Adjusts longitude values to ensure they are all within 180 degrees of the
    first coordinate, preventing artificial jumps when calculating polygon areas.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates in (longitude, latitude) format.

    Returns
    -------
    np.ndarray
        Array with unwrapped longitude coordinates, same shape as input.

    Notes
    -----
    This function uses the first longitude point as a reference and adjusts all
    subsequent longitudes to be within ±180 degrees of this reference point by
    adding or subtracting 360 degrees as needed.

    Examples
    --------
    >>> coords = np.array([[-170, 10], [170, 20], [-160, 30]])
    >>> unwrapped = unwrap_longitudes(coords)
    >>> # Longitudes adjusted to avoid large jumps across IDL
    """
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    coords_unwrapped = coords.copy()
    lons = coords_unwrapped[:, 0]
    ref_lon = lons[0]

    # Vectorized approach for better performance
    diff = lons - ref_lon
    # Adjust longitudes that are more than 180 degrees away
    lons[diff > 180] -= 360
    lons[diff < -180] += 360

    return coords_unwrapped


def convert_international_date_line_polygon_to_unwrapped_polygon(
    geometry_in: ogr.Geometry,
) -> Optional[ogr.Geometry]:
    """
    Convert a polygon crossing the International Date Line to a valid polygon.

    This function detects if a polygon crosses the International Date Line (IDL)
    by checking if its longitude range exceeds 180 degrees. If so, it adjusts
    the longitudes of vertices in the western hemisphere (< -150°) by adding
    360°, effectively shifting them to the eastern hemisphere to create a
    continuous polygon.

    Parameters
    ----------
    geometry_in : ogr.Geometry
        The input OGR polygon geometry to check and potentially convert.
        Must be a valid OGR Geometry object (typically a Polygon).

    Returns
    -------
    ogr.Geometry or None
        - If the polygon doesn't cross the IDL (longitude span < 180°):
          Returns the original geometry if valid, None if invalid.
        - If the polygon crosses the IDL (longitude span >= 180°):
          Returns a new polygon with adjusted longitudes that create a
          valid continuous polygon. Returns None if the converted geometry
          is invalid.

    Raises
    ------
    ValueError
        If geometry_in is None.
    AttributeError
        If geometry_in is not a valid OGR Geometry object.

    Notes
    -----
    **Algorithm**:

    1. Extract all coordinates from the input geometry
    2. Calculate longitude range (max - min)
    3. If range < 180°:
       - Polygon doesn't cross IDL
       - Return original if valid, None if invalid
    4. If range >= 180°:
       - Polygon likely crosses IDL
       - Shift longitudes < -150° by adding 360°
       - Reconstruct polygon with adjusted coordinates
       - Validate and return

    **IDL Crossing Detection**:

    A longitude range >= 180° indicates the polygon likely crosses the IDL.
    For example, a polygon with vertices at 170°E and 170°W has a range of
    340° in the [-180, 180] system, but represents a narrow strip across
    the IDL.

    **Longitude Adjustment**:

    The threshold of -150° is used to identify western hemisphere coordinates
    that should be shifted. This works for typical IDL-crossing polygons in
    the Pacific region:

    - Coordinates < -150° are shifted to [210°, 360°] range
    - Coordinates >= -150° remain unchanged
    - Result is a continuous polygon in approximately [0°, 360°] range

    **Spatial Reference**:

    The output geometry preserves the spatial reference system from the input
    geometry. Note that standard geographic CRS typically use [-180, 180]
    longitude range, so the output with [0, 360] range may need special
    handling in some applications.

    Warnings
    --------
    - This function assumes the polygon crosses the IDL if longitude range
      exceeds 180°. This heuristic may not work for all cases.
    - The -150° threshold is designed for Pacific region polygons and may
      not be appropriate for all IDL-crossing geometries.
    - The output geometry may have longitudes outside the standard [-180, 180]
      range, which could cause issues with some GIS software.
    - Only the outer ring is preserved; interior rings (holes) are not handled.

    Examples
    --------
    Convert a simple IDL-crossing polygon:

    >>> from osgeo import ogr
    >>> # Polygon spanning from 170°E to 170°W (crosses IDL)
    >>> wkt = "POLYGON((170 -10, 180 -10, -170 -10, -160 -10, 170 -10))"
    >>> geom = ogr.CreateGeometryFromWkt(wkt)
    >>> result = convert_idl_polygon_to_valid_polygon(geom)
    >>> # Result has continuous longitudes: 170, 180, 190, 200

    Handle a normal polygon that doesn't cross the IDL:

    >>> # Polygon entirely in western hemisphere
    >>> wkt = "POLYGON((-120 30, -100 30, -100 40, -120 40, -120 30))"
    >>> geom = ogr.CreateGeometryFromWkt(wkt)
    >>> result = convert_idl_polygon_to_valid_polygon(geom)
    >>> # Returns original geometry (range < 180)

    Handle a polygon near but not crossing the IDL:

    >>> # Polygon in eastern hemisphere near IDL
    >>> wkt = "POLYGON((160 20, 175 20, 175 25, 160 25, 160 20))"
    >>> geom = ogr.CreateGeometryFromWkt(wkt)
    >>> result = convert_idl_polygon_to_valid_polygon(geom)
    >>> # Returns original geometry (range = 15° < 180°)

    Convert a polygon with wide IDL crossing:

    >>> # Polygon spanning Pacific Ocean
    >>> wkt = "POLYGON((150 0, 180 0, -150 0, 150 0))"
    >>> geom = ogr.CreateGeometryFromWkt(wkt)
    >>> result = convert_idl_polygon_to_valid_polygon(geom)
    >>> # Shifts -150 to 210, creating continuous polygon

    See Also
    --------
    pyearth.gis.geometry.split_polygon_cross_idl : Split IDL polygons differently
    osgeo.ogr.Geometry.IsValid : Check if geometry is topologically valid
    osgeo.ogr.Geometry.GetSpatialReference : Get geometry's spatial reference
    """
    # Validate input
    if geometry_in is None:
        raise ValueError("geometry_in cannot be None")

    # Get the coordinates first
    coordinates_gcs = get_geometry_coordinates(geometry_in)

    # Check if we got valid coordinates
    if coordinates_gcs is None or len(coordinates_gcs) == 0:
        raise ValueError("Failed to extract coordinates from geometry")

    # Calculate longitude range
    longitude_min = np.min(coordinates_gcs[:, 0])
    longitude_max = np.max(coordinates_gcs[:, 0])
    longitude_range = longitude_max - longitude_min

    # If longitude range < 180, polygon doesn't cross IDL
    if longitude_range < 180:
        # Return original geometry if valid, None otherwise
        if geometry_in.IsValid():
            return geometry_in
        else:
            return None

    # Polygon crosses IDL - adjust longitudes
    # Shift western hemisphere longitudes (< -150°) to eastern hemisphere
    adjusted_coords = coordinates_gcs.copy()

    # Vectorized longitude adjustment for better performance
    western_mask = adjusted_coords[:, 0] < -150
    adjusted_coords[western_mask, 0] += 360.0

    # Reconstruct polygon with adjusted coordinates
    geometry_out = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)

    for coord in adjusted_coords:
        ring.AddPoint(float(coord[0]), float(coord[1]))

    # Close the ring (ensures first and last points are the same)
    ring.CloseRings()

    # Add ring to polygon
    geometry_out.AddGeometry(ring)

    # Preserve spatial reference from input
    spatial_ref = geometry_in.GetSpatialReference()
    if spatial_ref is not None:
        geometry_out.AssignSpatialReference(spatial_ref)

    # Validate the output geometry?
    if not geometry_out.IsValid():
        return None

    return geometry_out


def split_international_date_line_polygon_coordinates(
    aCoord_gcs: np.ndarray,
) -> Optional[List[np.ndarray]]:
    """
    Split a polygon crossing the International Date Line into left and right hemisphere polygons.

    This function takes a polygon that actually crosses the IDL (180°/-180° longitude) and splits it
    into two separate polygons: one for the Eastern hemisphere (positive longitudes) and one
    for the Western hemisphere (negative longitudes). The split is performed along the IDL
    meridian using great circle intersection calculations.

    **Important**: This function should only be called for polygons that have actual edge crossings
    of the IDL. Polygons with vertices exactly on the IDL (±180°) but no edge crossings should be
    handled by `check_cross_international_date_line_polygon` which returns adjusted coordinates.

    Parameters
    ----------
    aCoord_gcs : np.ndarray
        Input polygon coordinates as Nx2 array of (longitude, latitude) pairs in decimal degrees.

        Requirements:
        - Should be a closed polygon (first point equals last point)
        - Vertices should be in counter-clockwise (CCW) order
        - Must have actual edges crossing the IDL (exactly 2 crossings)
        - Minimum 4 points (3 unique + 1 closure)
        - No vertices exactly on ±180° meridian (handled separately)

        Format: [[lon1, lat1], [lon2, lat2], ..., [lonN, latN]]

        Longitude range: [-180, 180] degrees
        Latitude range: [-90, 90] degrees

    Returns
    -------
    Optional[List[np.ndarray]]
        List containing two polygons: [left_polygon, right_polygon]

        - **left_polygon** (numpy.ndarray): Eastern hemisphere polygon (lon > 0)
          - Coordinates with positive longitudes
          - Includes boundary points near +180°
          - Closed polygon (first == last)

        - **right_polygon** (numpy.ndarray): Western hemisphere polygon (lon < 0)
          - Coordinates with negative longitudes
          - Includes boundary points near -180°
          - Closed polygon (first == last)

    Raises
    ------
    ValueError
        If exactly 2 IDL crossings are not found, or if input coordinates are invalid.
    TypeError
        If input is not a valid array-like structure.

    Notes
    -----
    **Algorithm Overview:**

    1. **Preprocessing**: Excludes vertices exactly on IDL from crossing detection
    2. **Edge Crossing Detection**: Finds edges that transition between hemispheres
    3. **Intersection Calculation**: Uses great circle math to find exact crossing latitudes
    4. **Boundary Addition**: Adds meridian boundary points for polygon closure
    5. **Vertex Partitioning**: Assigns vertices to eastern/western polygons
    6. **Polygon Closure**: Ensures both result polygons are properly closed

    **IDL Vertex Handling:**

    - Vertices exactly on ±180° are excluded from crossing detection
    - Such cases should be pre-processed by `check_cross_international_date_line_polygon`
    - This function focuses only on actual edge crossings

    **Numerical Precision:**

    - Uses 1e-8 offset from ±180° for boundary points
    - 1e-10 tolerance for closure checks
    - Great circle calculations for accurate intersection latitudes

    Warnings
    --------
    - Only handles polygons with exactly 2 IDL edge crossings
    - Does not validate CCW order - incorrect order will produce wrong results
    - Does not handle polygons with holes or multi-part geometries
    - May fail for self-intersecting or very complex polygons
    - Vertices exactly on IDL should be handled before calling this function

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

    convert_idl_polygon_to_valid_polygon : Longitude shifting approach
    find_great_circle_intersection_with_meridian : Great circle intersection calculation
    """
    # Input validation
    try:
        aCoord_gcs = np.asarray(aCoord_gcs, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(f"Input coordinates must be convertible to numeric array: {e}")

    if aCoord_gcs.ndim != 2 or aCoord_gcs.shape[1] != 2:
        raise ValueError("Input must be a 2D array with shape (n, 2)")

    if len(aCoord_gcs) < 4:
        raise ValueError("Polygon must have at least 4 points (3 unique + closure)")

    # Validate longitude and latitude ranges
    lons, lats = aCoord_gcs[:, 0], aCoord_gcs[:, 1]
    if not np.all((-180 <= lons) & (lons <= 180)):
        raise ValueError("Longitudes must be in range [-180, 180]")
    if not np.all((-90 <= lats) & (lats <= 90)):
        raise ValueError("Latitudes must be in range [-90, 90]")

    # Ensure counter-clockwise order
    if not check_counter_clockwise_local(aCoord_gcs):
        aCoord_gcs = aCoord_gcs[::-1]

    # Get number of points in the polygon
    nPoint = len(aCoord_gcs)

    # Find indices where edges cross the International Date Line
    # The algorithm requires coordinates to be in counter-clockwise (CCW) order
    aIndex = []

    # Vectorized approach to find IDL crossings for better performance
    longitudes = aCoord_gcs[:, 0]
    longitudes_next = np.roll(longitudes, -1)  # Next longitude with wrap-around

    # Check if any point lies exactly on the IDL
    # Exclude the last point to avoid duplication in closed polygons (first == last)
    idl_points = np.zeros_like(longitudes, dtype=bool)
    idl_points[:-1] = (longitudes[:-1] == 180.0) | (longitudes[:-1] == -180.0)

    # More explicit IDL point detection for crossing logic
    current_on_idl = (np.abs(longitudes - 180.0) < 1e-10) | (
        np.abs(longitudes + 180.0) < 1e-10
    )
    next_on_idl = (np.abs(longitudes_next - 180.0) < 1e-10) | (
        np.abs(longitudes_next + 180.0) < 1e-10
    )

    # Find eastward crossings: positive to negative longitude (0° → 180° → -180°)
    # Exclude edges where either vertex is exactly on the IDL
    eastward_crossings = (
        (longitudes > 0)
        & (longitudes < 180.0)
        & (longitudes_next < 0)
        & ~current_on_idl
        & ~next_on_idl
    )

    # Find westward crossings: negative to positive longitude (-180° → 180° → 0°)
    # Exclude edges where either vertex is exactly on the IDL
    westward_crossings = (
        (longitudes < 0)
        & (longitudes > -180.0)
        & (longitudes_next > 0)
        & ~current_on_idl
        & ~next_on_idl
    )

    # Combine crossing conditions
    crossing_mask = idl_points | eastward_crossings | westward_crossings
    crossing_indices = np.where(crossing_mask)[0]

    # Convert to list for compatibility with existing code
    aIndex = crossing_indices.tolist()

    # Find vertices that are exactly on the IDL (±180°)
    idl_vertex_indices = np.where(idl_points)[0].tolist()

    # Verify we found exactly 2 crossings for normal splitting
    # A polygon crossing the IDL should have exactly 2 edges that cross it
    if len(aIndex) != 2:
        raise ValueError(
            f"Expected exactly 2 IDL crossings, found {len(aIndex)}. "
            f"This polygon may not cross the IDL properly or may have a complex crossing pattern."
        )

    # Helper function to calculate intersection latitude for an edge
    def calculate_intersection_latitude(edge_index: int) -> float:
        """Calculate intersection latitude for given edge crossing IDL."""
        start_lon = aCoord_gcs[edge_index, 0]
        start_lat = aCoord_gcs[edge_index, 1]

        # Handle wrap-around for closing edge
        if edge_index == nPoint - 1:
            end_lon = aCoord_gcs[0, 0]
            end_lat = aCoord_gcs[0, 1]
        else:
            end_lon = aCoord_gcs[edge_index + 1, 0]
            end_lat = aCoord_gcs[edge_index + 1, 1]

        target_longitude = 180.0

        # Determine correct order for intersection calculation
        if start_lon > 0 and end_lon < 0:
            _, intersection_lat = find_great_circle_intersection_with_meridian(
                start_lon, start_lat, end_lon, end_lat, target_longitude
            )
        else:
            _, intersection_lat = find_great_circle_intersection_with_meridian(
                end_lon, end_lat, start_lon, start_lat, target_longitude
            )

        return intersection_lat

    # Calculate intersection latitudes for both crossings
    intersection_lat_0 = calculate_intersection_latitude(aIndex[0])
    intersection_lat_1 = calculate_intersection_latitude(aIndex[1])

    # Determine northern and southern intersection points
    northern_lat = max(intersection_lat_0, intersection_lat_1)
    southern_lat = min(intersection_lat_0, intersection_lat_1)

    # Initialize containers for the split polygons
    eastern_polygon_coords = []  # Positive longitudes (Eastern hemisphere)
    western_polygon_coords = []  # Negative longitudes (Western hemisphere)

    # Constants for numerical stability
    LONGITUDE_OFFSET = 1.0e-8
    EASTERN_BOUNDARY = 180.0 - LONGITUDE_OFFSET
    WESTERN_BOUNDARY = -180.0 + LONGITUDE_OFFSET

    # Track if boundary points have been added to avoid duplication
    eastern_boundary_added = False
    western_boundary_added = False

    # Minimum crossing index for vertex partitioning logic
    min_crossing_index = min(aIndex)

    # Partition vertices based on longitude and crossing logic
    for vertex_index in range(nPoint):
        longitude = aCoord_gcs[vertex_index, 0]
        is_crossing_vertex = vertex_index in aIndex
        is_before_first_crossing = vertex_index <= min_crossing_index

        # Process Eastern hemisphere vertices (positive longitude)
        if longitude > 0:
            # Apply offset if vertex is exactly on +180° IDL
            if abs(longitude - 180.0) < 1e-10:
                coord_to_add = [EASTERN_BOUNDARY, aCoord_gcs[vertex_index, 1]]
            else:
                coord_to_add = aCoord_gcs[vertex_index]
            eastern_polygon_coords.append(coord_to_add)

            # Add boundary points at crossing vertices
            if is_crossing_vertex and not eastern_boundary_added:
                eastern_boundary_added = True
                # Add meridian boundary points (south to north)
                eastern_polygon_coords.extend(
                    [[EASTERN_BOUNDARY, southern_lat], [EASTERN_BOUNDARY, northern_lat]]
                )

        # Process Western hemisphere vertices (negative longitude)
        elif longitude < 0:
            # Apply offset if vertex is exactly on -180° IDL
            if abs(longitude + 180.0) < 1e-10:
                coord_to_add = [WESTERN_BOUNDARY, aCoord_gcs[vertex_index, 1]]
            else:
                coord_to_add = aCoord_gcs[vertex_index]
            western_polygon_coords.append(coord_to_add)

            # Add boundary points at crossing vertices
            if is_crossing_vertex and not western_boundary_added:
                western_boundary_added = True
                # Add meridian boundary points (north to south for correct winding)
                western_polygon_coords.extend(
                    [[WESTERN_BOUNDARY, northern_lat], [WESTERN_BOUNDARY, southern_lat]]
                )

    # Convert to numpy arrays and ensure proper closure
    def ensure_polygon_closure(coords_list: list) -> np.ndarray:
        """Convert coordinate list to numpy array and ensure polygon is closed."""
        if not coords_list:
            return np.array([])

        coords_array = np.array(coords_list)

        # Close polygon if not already closed
        if not np.allclose(coords_array[0], coords_array[-1], atol=1e-10):
            coords_array = np.vstack([coords_array, coords_array[0:1]])

        return coords_array

    eastern_polygon = ensure_polygon_closure(eastern_polygon_coords)
    western_polygon = ensure_polygon_closure(western_polygon_coords)

    # Ensure counter-clockwise orientation for both polygons
    if len(eastern_polygon) > 0 and not check_counter_clockwise_local(eastern_polygon):
        eastern_polygon = eastern_polygon[::-1]
    if len(western_polygon) > 0 and not check_counter_clockwise_local(western_polygon):
        western_polygon = western_polygon[::-1]

    return [eastern_polygon, western_polygon]


def check_cross_international_date_line_polygon(
    coords: np.ndarray,
) -> Tuple[bool, np.ndarray]:
    """Check if polygon coordinates cross the International Date Line.

    This function distinguishes between actual IDL edge crossings and vertices
    that happen to lie exactly on the IDL meridian (±180°). Vertices on the IDL
    without edge crossings are adjusted slightly to avoid numerical issues.

    Parameters
    ----------
    coords : np.ndarray
        Array of polygon coordinates in (longitude, latitude) format.

    Returns
    -------
    Tuple[bool, np.ndarray]
        A tuple containing:
        - bool: True if the polygon has actual edges crossing the IDL, False otherwise
        - np.ndarray or None:
          - If False: Adjusted coordinates with IDL vertices moved slightly off the meridian
          - If True: None (original coordinates should be used for splitting)

    Notes
    -----
    **IDL Detection Logic:**

    1. **Vertex Detection**: Identifies vertices exactly on ±180° meridian
    2. **Edge Crossing Detection**: Finds edges that actually cross the IDL
       (excludes edges involving IDL vertices)
    3. **Coordinate Adjustment**: If no edge crossings but IDL vertices exist,
       moves these vertices slightly (±1e-6°) based on neighboring vertices
    4. **Hemisphere Assignment**: IDL vertices are moved to the hemisphere
       containing the majority of their neighbors

    **Return Behavior:**

    - Returns (False, adjusted_coords) for polygons with IDL vertices but no crossings
    - Returns (True, None) for polygons with actual IDL edge crossings
    - Returns (False, coords) for polygons that don't involve the IDL

    This approach ensures that single IDL vertices don't cause false crossing
    detection while preserving the ability to split truly crossing polygons.
    """
    try:
        if (
            not isinstance(coords, np.ndarray)
            or coords.ndim != 2
            or coords.shape[1] != 2
        ):
            raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

        if len(coords) < 3:
            return False, coords.copy()

        # Validate coordinate ranges
        lons, lats = coords[:, 0], coords[:, 1]
        if not (
            np.all((-180 <= lons) & (lons <= 180))
            and np.all((-90 <= lats) & (lats <= 90))
        ):
            return False, coords.copy()

    except (ValueError, IndexError, TypeError):
        return False, coords.copy()  # Invalid coordinates can't cross IDL

    # Make a copy to avoid modifying the original
    coords_updated = coords.copy()
    lons = coords_updated[:, 0]

    # Check if any vertices are exactly on the IDL
    # Exclude the last point to avoid duplication in closed polygons (first == last)
    idl_vertices = np.zeros_like(lons, dtype=bool)
    idl_vertices[:-1] = np.abs(np.abs(lons[:-1]) - 180.0) < 1e-10

    # If there are vertices on IDL, check if there are actual edge crossings
    if np.any(idl_vertices):
        lons_next = np.roll(lons, -1)

        # Find eastward crossings: positive to negative longitude, excluding IDL vertices
        eastward_crossings = (
            (lons > 0)
            & (lons < 180.0)
            & (lons_next < 0)
            & ~idl_vertices
            & ~np.roll(idl_vertices, -1)
        )

        # Find westward crossings: negative to positive longitude, excluding IDL vertices
        westward_crossings = (
            (lons < 0) & (lons_next > 0) & ~idl_vertices & ~np.roll(idl_vertices, -1)
        )

        # If no actual edge crossings, adjust IDL vertices and return False
        if not (np.any(eastward_crossings) or np.any(westward_crossings)):
            # Adjust vertices that are exactly on the IDL by slightly moving them
            IDL_OFFSET = 1e-6  # Small offset to move points just off the IDL
            idl_indices = np.where(idl_vertices)[0]

            for idx in idl_indices:
                # Determine which hemisphere to move to based on neighboring vertices
                prev_idx = (idx - 1) % len(coords_updated)
                next_idx = (idx + 1) % len(coords_updated)

                prev_lon = coords_updated[prev_idx, 0]
                next_lon = coords_updated[next_idx, 0]

                # Count non-IDL neighbors to determine preferred hemisphere
                neighbor_lons = []
                if abs(abs(prev_lon) - 180.0) > 1e-10:  # Previous vertex not on IDL
                    neighbor_lons.append(prev_lon)
                if abs(abs(next_lon) - 180.0) > 1e-10:  # Next vertex not on IDL
                    neighbor_lons.append(next_lon)

                if neighbor_lons:
                    # Move to hemisphere containing majority of neighbors
                    positive_neighbors = sum(1 for lon in neighbor_lons if lon > 0)
                    if positive_neighbors >= len(neighbor_lons) / 2:
                        # Move to eastern hemisphere
                        coords_updated[idx, 0] = 180.0 - IDL_OFFSET
                    else:
                        # Move to western hemisphere
                        coords_updated[idx, 0] = -180.0 + IDL_OFFSET
                else:
                    # If neighbors are also on IDL, check overall polygon distribution
                    positive_lons = np.sum(lons > 0)
                    negative_lons = np.sum(lons < 0)

                    if positive_lons >= negative_lons:
                        coords_updated[idx, 0] = 180.0 - IDL_OFFSET
                    else:
                        coords_updated[idx, 0] = -180.0 + IDL_OFFSET

            # After adjustment, no crossing
            return False, coords_updated
        else:
            # There are actual edge crossings
            return True, None

    # Original logic for cases without IDL vertices
    # Check for large jumps between consecutive points (> 180 degrees)
    lon_diffs = np.abs(np.diff(lons))
    max_jump = np.max(lon_diffs)

    # Also check the wrap-around from last to first point
    wrap_jump = abs(lons[-1] - lons[0])

    crossing = max_jump > 180 or wrap_jump > 180
    return crossing, None


# Alias for backward compatibility
check_cross_idl = check_cross_international_date_line_polygon


def check_cross_international_date_line_geometry(geometry_in: ogr.Geometry) -> bool:
    """Check if an OGR geometry crosses the International Date Line.

    This function provides comprehensive IDL crossing detection for various OGR geometry types
    including Polygon, MultiPolygon, LineString, MultiLineString, Point, MultiPoint, and
    GeometryCollection. It uses optimized algorithms and proper error handling.

    Parameters
    ----------
    geometry_in : ogr.Geometry
        Input OGR geometry to check. Must be a valid OGR Geometry object.

    Returns
    -------
    bool
        True if the geometry crosses the IDL, False otherwise.

        - For POLYGON/MULTIPOLYGON: Uses robust coordinate analysis
        - For LINESTRING/MULTILINESTRING: Checks line segments for IDL crossings
        - For POINT/MULTIPOINT: Always returns False (points can't cross IDL)
        - For GEOMETRYCOLLECTION: Recursively checks all sub-geometries

    Raises
    ------
    ValueError
        If geometry_in is None, invalid, or if coordinate extraction fails.
    TypeError
        If geometry_in is not an OGR Geometry object.

    Notes
    -----
    **Algorithm Improvements**:

    1. **Comprehensive Type Support**: Handles all major OGR geometry types
    2. **Optimized Processing**: Early returns and vectorized operations where possible
    3. **Robust Error Handling**: Validates geometry before processing
    4. **Recursive Support**: Properly handles nested geometries in collections
    5. **Performance**: Caches geometry type checks and minimizes coordinate extraction

    **IDL Crossing Logic**:

    - **Polygons**: Uses coordinate span analysis and consecutive point jumps
    - **LineStrings**: Examines individual segments for longitude sign changes
    - **Multi-geometries**: Returns True if ANY sub-geometry crosses IDL
    - **Collections**: Recursively processes all contained geometries

    **Performance Optimizations**:

    - Early validation to avoid unnecessary processing
    - Geometry type caching to minimize OGR calls
    - Vectorized coordinate analysis where applicable
    - Short-circuit evaluation for multi-part geometries

    Examples
    --------
    Check a simple polygon:

    >>> from osgeo import ogr
    >>> wkt = "POLYGON((170 -10, 180 -10, -170 -10, -160 -10, 170 -10))"
    >>> geom = ogr.CreateGeometryFromWkt(wkt)
    >>> check_cross_international_date_line_geometry(geom)
    True

    Check a non-crossing polygon:

    >>> wkt = "POLYGON((10 10, 20 10, 20 20, 10 20, 10 10))"
    >>> geom = ogr.CreateGeometryFromWkt(wkt)
    >>> check_cross_international_date_line_geometry(geom)
    False

    Check a multipolygon:

    >>> wkt = "MULTIPOLYGON(((170 0, 180 0, -170 0, 170 0)), ((10 0, 20 0, 20 10, 10 10, 10 0)))"
    >>> geom = ogr.CreateGeometryFromWkt(wkt)
    >>> check_cross_international_date_line_geometry(geom)
    True

    Check a linestring:

    >>> wkt = "LINESTRING(170 0, -170 0)"
    >>> geom = ogr.CreateGeometryFromWkt(wkt)
    >>> check_cross_international_date_line_geometry(geom)
    True

    See Also
    --------
    check_cross_international_date_line_polygon : Polygon-specific IDL checking
    split_international_date_line_polygon_coordinates : Split IDL-crossing polygons
    """
    # Input validation
    if geometry_in is None:
        raise ValueError("geometry_in cannot be None")

    # Verify it's an OGR Geometry object
    if not hasattr(geometry_in, "GetGeometryName"):
        raise TypeError("geometry_in must be an OGR Geometry object")

    # Validate geometry before processing
    if geometry_in.IsEmpty():
        return False

    # Cache geometry type to minimize OGR calls
    sGeometry_type = geometry_in.GetGeometryName()

    # Handle different geometry types with optimized logic
    if sGeometry_type == "POLYGON":
        return _check_polygon_idl_crossing(geometry_in)

    elif sGeometry_type == "MULTIPOLYGON":
        return _check_multipolygon_idl_crossing(geometry_in)

    elif sGeometry_type == "LINESTRING":
        return _check_linestring_idl_crossing(geometry_in)

    elif sGeometry_type == "MULTILINESTRING":
        return _check_multilinestring_idl_crossing(geometry_in)

    elif sGeometry_type in ["POINT", "MULTIPOINT"]:
        # Points cannot cross the IDL by definition
        return False

    elif sGeometry_type == "GEOMETRYCOLLECTION":
        return _check_geometry_collection_idl_crossing(geometry_in)

    else:
        # Handle unknown or unsupported geometry types
        raise ValueError(
            f"Unsupported geometry type: {sGeometry_type}. "
            f"Supported types: POLYGON, MULTIPOLYGON, LINESTRING, "
            f"MULTILINESTRING, POINT, MULTIPOINT, GEOMETRYCOLLECTION"
        )


def _check_polygon_idl_crossing(polygon: ogr.Geometry) -> bool:
    """Check if a single polygon crosses the IDL."""
    try:
        coordinates_gcs = get_geometry_coordinates(polygon)
        if coordinates_gcs is None or len(coordinates_gcs) == 0:
            raise ValueError("Failed to extract coordinates from polygon")

        crosses_idl, _ = check_cross_international_date_line_polygon(coordinates_gcs)
        return crosses_idl
    except Exception as e:
        raise ValueError(f"Error processing polygon: {e}")


def _check_multipolygon_idl_crossing(multipolygon: ogr.Geometry) -> bool:
    """Check if any polygon in a multipolygon crosses the IDL."""
    try:
        nGeometries = multipolygon.GetGeometryCount()
        if nGeometries == 0:
            return False

        # Use short-circuit evaluation - return True as soon as one crossing is found
        for i in range(nGeometries):
            sub_geometry = multipolygon.GetGeometryRef(i)
            if sub_geometry is None:
                continue

            # Recursively check each polygon
            if check_cross_international_date_line_geometry(sub_geometry):
                return True
        return False
    except Exception as e:
        raise ValueError(f"Error processing multipolygon: {e}")


def _check_linestring_idl_crossing(linestring: ogr.Geometry) -> bool:
    """Check if a linestring crosses the IDL."""
    try:
        coordinates_gcs = get_geometry_coordinates(linestring)
        if coordinates_gcs is None or len(coordinates_gcs) < 2:
            return False

        # Check for large longitude jumps between consecutive points
        lons = coordinates_gcs[:, 0]

        # Vectorized approach for better performance
        lon_diffs = np.abs(np.diff(lons))
        return np.any(lon_diffs > 180.0)

    except Exception as e:
        raise ValueError(f"Error processing linestring: {e}")


def _check_multilinestring_idl_crossing(multilinestring: ogr.Geometry) -> bool:
    """Check if any linestring in a multilinestring crosses the IDL."""
    try:
        nGeometries = multilinestring.GetGeometryCount()
        if nGeometries == 0:
            return False

        # Use short-circuit evaluation
        for i in range(nGeometries):
            sub_geometry = multilinestring.GetGeometryRef(i)
            if sub_geometry is None:
                continue

            # Recursively check each linestring
            if check_cross_international_date_line_geometry(sub_geometry):
                return True
        return False
    except Exception as e:
        raise ValueError(f"Error processing multilinestring: {e}")


def _check_geometry_collection_idl_crossing(collection: ogr.Geometry) -> bool:
    """Check if any geometry in a collection crosses the IDL."""
    try:
        nGeometries = collection.GetGeometryCount()
        if nGeometries == 0:
            return False

        # Use short-circuit evaluation
        for i in range(nGeometries):
            sub_geometry = collection.GetGeometryRef(i)
            if sub_geometry is None:
                continue

            # Recursively check each geometry in the collection
            if check_cross_international_date_line_geometry(sub_geometry):
                return True
        return False
    except Exception as e:
        raise ValueError(f"Error processing geometry collection: {e}")


def _validate_coordinate_array(coords: np.ndarray, min_points: int = 3) -> None:
    """Validate coordinate array format and content.

    Parameters
    ----------
    coords : np.ndarray
        Array to validate
    min_points : int, default=3
        Minimum number of points required

    Raises
    ------
    ValueError
        If coordinates are invalid
    """
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    if len(coords) < min_points:
        raise ValueError(f"Coordinates must have at least {min_points} points")

    # Validate longitude and latitude ranges
    lons, lats = coords[:, 0], coords[:, 1]
    if not np.all((-180 <= lons) & (lons <= 180)):
        raise ValueError("Longitudes must be in range [-180, 180]")
    if not np.all((-90 <= lats) & (lats <= 90)):
        raise ValueError("Latitudes must be in range [-90, 90]")


def check_counter_clockwise_local(coords: np.ndarray) -> bool:
    """Local implementation of counter-clockwise check to avoid circular imports.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates.

    Returns
    -------
    bool
        True if vertices are in counter-clockwise order, False otherwise.
    """
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    if len(coords) < 3:
        return True  # Degenerate case

    # Check if polygon crosses the International Date Line
    crosses_idl, coords_adjusted = check_cross_international_date_line_polygon(coords)
    if crosses_idl:
        coords_unwrapped = unwrap_longitudes(coords)
        # Calculate signed area using optimized shoelace formula
        signed_area = calculate_signed_area_shoelace(coords_unwrapped)
    else:
        # Standard case: calculate signed area directly
        # Use adjusted coordinates if available (for cases where IDL vertices were moved)
        coords_to_use = coords_adjusted if coords_adjusted is not None else coords
        signed_area = calculate_signed_area_shoelace(coords_to_use)
    return signed_area > 0


def calculate_signed_area_shoelace(coords: np.ndarray) -> float:
    """Calculate the signed area of a polygon using the shoelace formula.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates.

    Returns
    -------
    float
        Signed area of the polygon. Positive for counter-clockwise, negative for clockwise.
    """
    x, y = coords[:, 0], coords[:, 1]

    # Vectorized shoelace formula - more efficient than loops
    # Handle the wrap-around (last point to first point) implicitly
    x_rolled = np.roll(x, -1)  # x[i+1] for all i, with wrap-around
    y_rolled = np.roll(y, -1)  # y[i+1] for all i, with wrap-around

    signed_area = 0.5 * np.sum(x * y_rolled - x_rolled * y)
    return signed_area


def reorder_international_date_line_polygon_vertices(
    vertices: List[Coord],
) -> List[Coord]:
    """Reorder polygon vertices to create a valid geometry by rotation.

    Parameters
    ----------
    vertices : List[Coord]
        List of coordinate tuples in (longitude, latitude) format.

    Returns
    -------
    List[Coord]
        Reordered vertices that form a valid polygon.

    Raises
    ------
    ValueError
        If vertices cannot be converted to valid coordinates or if no valid
        polygon can be formed by rotation.
    """
    # Handle empty input
    if not vertices:
        return []

    # Validate input format
    try:
        # Check that each vertex is a valid coordinate pair
        for i, vertex in enumerate(vertices):
            if not hasattr(vertex, "__len__") or len(vertex) != 2:
                raise ValueError(f"Vertex {i} must be a coordinate pair (lon, lat)")
            lon, lat = float(vertex[0]), float(vertex[1])
            if not (-180 <= lon <= 180):
                raise ValueError(
                    f"Longitude {lon} at vertex {i} out of range [-180, 180]"
                )
            if not (-90 <= lat <= 90):
                raise ValueError(f"Latitude {lat} at vertex {i} out of range [-90, 90]")
    except (TypeError, ValueError, IndexError) as e:
        raise ValueError(f"Invalid vertex format: {e}")

    if len(vertices) < 3:
        raise ValueError("Polygon must have at least 3 vertices")

    # Ensure polygon is closed; append first vertex if not closed.
    if vertices and (vertices[0] != vertices[-1]):
        vertices = vertices + [vertices[0]]

    unique_vertices = vertices[:-1]
    n_points = len(unique_vertices)
    current_vertices = unique_vertices.copy()

    # Use a rotating method to find the polygon that is valid
    is_valid = False
    rotation_count = 0

    while not is_valid and rotation_count <= n_points:
        rotation_count += 1

        try:
            polygon = ogr.Geometry(ogr.wkbPolygon)
            linear_ring = ogr.Geometry(ogr.wkbLinearRing)
            for lon, lat in current_vertices:
                linear_ring.AddPoint(float(lon), float(lat))
            linear_ring.CloseRings()
            polygon.AddGeometry(linear_ring)

            if polygon.IsValid():
                is_valid = True
                break
        except Exception:
            # OGR geometry creation failed, try next rotation
            pass

        # Rotate the first point to the end
        if rotation_count <= n_points:
            current_vertices = current_vertices[1:] + [current_vertices[0]]

    if not is_valid:
        raise ValueError(
            "Cannot form a valid polygon by rotating vertices. "
            "The input coordinates may be self-intersecting or invalid."
        )

    new_vertices = current_vertices + [current_vertices[0]]
    return new_vertices
