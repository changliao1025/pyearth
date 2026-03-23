import numpy as np
from osgeo import ogr
from typing import List, Tuple, Optional
from pyearth.gis.geometry.calculate_intersect_on_great_circle import (
    find_great_circle_intersection_with_meridian,
)
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.pole_check import polygon_includes_pole
Coord = Tuple[float, float]

# Module-level constants for IDL detection and adjustment
IDL_TOLERANCE = 1e-6  # Tolerance for detecting points on the IDL (±180°)
IDL_OFFSET = 1e-7     # Offset to move points off IDL (< IDL_TOLERANCE so they remain detectable)


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
    Split a polygon crossing the International Date Line into eastern and western polygons.

    Takes a polygon that actually crosses the IDL (180°/−180° longitude) and splits it into
    two closed polygons: one for the Eastern hemisphere (positive longitudes, bounded by
    +180°) and one for the Western hemisphere (negative longitudes, bounded by −180°).
    The split boundary is computed using great-circle intersection with the ±180° meridian.

    **Important**: Call this function only for polygons with actual *edge* crossings of the
    IDL.  Polygons whose vertices merely touch ±180° without crossing should first be
    normalised by :func:`check_cross_international_date_line_polygon`.

    Parameters
    ----------
    aCoord_gcs : np.ndarray
        Polygon coordinates as an Nx2 array of ``(longitude, latitude)`` pairs in decimal
        degrees.

        Requirements:

        - Closed polygon: first point equals last point.
        - Minimum 4 points (3 unique vertices + closure).
        - Longitudes in ``[−180, 180]``, latitudes in ``[−90, 90]``.
        - Exactly 2 edges that cross the IDL (no vertices exactly on ±180°).

    Returns
    -------
    list of np.ndarray
        ``[eastern_polygon, western_polygon]``

        *eastern_polygon* — closed Mx2 array with all longitudes near +180°.
        *western_polygon* — closed Kx2 array with all longitudes near −180°.

        Either element may be an empty array (shape ``(0,)``) when the polygon only
        touches the IDL without crossing it.

    Raises
    ------
    TypeError
        If *aCoord_gcs* cannot be converted to a numeric array.
    ValueError
        If the array shape is wrong, coordinate ranges are violated, or the number of
        detected IDL edge crossings is not exactly 2.

    Notes
    -----
    **Algorithm**

    1. *Validate & orient* — coerce to float64, check shape/ranges, reverse to CCW if needed.
    2. *IDL-vertex normalisation* — if vertices lie exactly on ±180° but no edge crosses,
       nudge them into the dominant hemisphere and return a single-polygon result.
    3. *Edge-crossing detection* — find the two edges whose endpoints straddle the IDL
       (positive→negative or negative→positive longitude), excluding edges that start or
       end exactly on ±180°.
    4. *Intersection calculation* — compute the great-circle latitude at which each
       crossing edge pierces the ±180° meridian.
    5. *Edge-by-edge traversal* — walk every edge of the polygon.  At each IDL crossing,
       append the intersection point to the *current* sub-polygon, then switch the active
       sub-polygon and open it with the same intersection point (opposite boundary
       longitude) before continuing with the next vertex.
    6. *Closure & orientation* — close both polygons and ensure CCW winding.

    **Why edge-by-edge traversal?**

    The previous vertex-only loop used a single ``eastern_boundary_added`` flag that
    prevented the second crossing from inserting its own boundary point, and hardcoded a
    global ``[southern_lat, northern_lat]`` pair regardless of traversal direction.  The
    edge-by-edge approach ties each intersection point to its specific crossing edge,
    producing geometrically correct results for all CCW polygons.

    **Numerical precision**

    - Boundary longitude offset: ``IDL_OFFSET = 1e-7`` (module constant).
    - Polygon closure tolerance: ``1e-10``.
    - Great-circle intersection via :func:`find_great_circle_intersection_with_meridian`.

    Warnings
    --------
    - Only handles polygons with exactly 2 IDL edge crossings.
    - Does not handle polygons with holes or multi-part geometries.
    - Vertices exactly on ±180° should be normalised before calling this function.

    Examples
    --------
    >>> import numpy as np
    >>> coords = np.array([
    ...     [170.0, -10.0],
    ...     [170.0,  10.0],
    ...     [-170.0, 10.0],
    ...     [-170.0,-10.0],
    ...     [170.0, -10.0],
    ... ])
    >>> eastern, western = split_international_date_line_polygon_coordinates(coords)
    >>> np.all(eastern[:, 0] > 0)   # all eastern longitudes positive
    True
    >>> np.all(western[:, 0] < 0)   # all western longitudes negative
    True
    >>> np.allclose(eastern[0], eastern[-1])  # closed
    True
    >>> np.allclose(western[0], western[-1])  # closed
    True

    See Also
    --------
    check_cross_international_date_line_polygon : Detect and normalise IDL-touching polygons.
    convert_international_date_line_polygon_to_unwrapped_polygon : Longitude-shift approach.
    find_great_circle_intersection_with_meridian : Great-circle meridian intersection.
    """
    # ------------------------------------------------------------------
    # 1. Input validation
    # ------------------------------------------------------------------
    try:
        aCoord_gcs = np.asarray(aCoord_gcs, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(f"Input coordinates must be convertible to numeric array: {e}")

    if aCoord_gcs.ndim != 2 or aCoord_gcs.shape[1] != 2:
        raise ValueError("Input must be a 2D array with shape (n, 2)")

    if len(aCoord_gcs) < 4:
        raise ValueError("Polygon must have at least 4 points (3 unique + closure)")

    lons = aCoord_gcs[:, 0]
    lats = aCoord_gcs[:, 1]
    if not np.all((-180 <= lons) & (lons <= 180)):
        raise ValueError("Longitudes must be in range [-180, 180]")
    if not np.all((-90 <= lats) & (lats <= 90)):
        raise ValueError("Latitudes must be in range [-90, 90]")

    # Ensure counter-clockwise winding so the traversal direction is deterministic.
    if not check_counter_clockwise_local(aCoord_gcs):
        aCoord_gcs = aCoord_gcs[::-1]

    nPoint = len(aCoord_gcs)
    lons = aCoord_gcs[:, 0]  # refresh after possible reversal

    # ------------------------------------------------------------------
    # 2. Classify vertices and edges
    # ------------------------------------------------------------------
    lons_next = np.roll(lons, -1)  # longitude of the *next* vertex for each edge i→i+1

    # Vertices that lie exactly on the IDL (±180°).
    # Exclude the closing duplicate (last == first) to avoid double-counting.
    on_idl = np.zeros(nPoint, dtype=bool)
    on_idl[:-1] = np.abs(np.abs(lons[:-1]) - 180.0) < IDL_TOLERANCE

    # Edges that actually cross the IDL: one endpoint strictly positive, the other
    # strictly negative, and neither endpoint sits exactly on the IDL.
    not_on_idl = ~on_idl
    not_next_on_idl = ~np.roll(on_idl, -1)

    eastward = (  # east→west: positive lon → negative lon
        (lons > 0) & (lons < 180.0) & (lons_next < 0) & not_on_idl & not_next_on_idl
    )
    westward = (  # west→east: negative lon → positive lon
        (lons < 0) & (lons > -180.0) & (lons_next > 0) & not_on_idl & not_next_on_idl
    )
    crossing_edges = eastward | westward
    crossing_edge_indices = np.where(crossing_edges)[0].tolist()

    idl_vertex_indices = np.where(on_idl)[0].tolist()

    # ------------------------------------------------------------------
    # 3. IDL-vertex-only path: nudge vertices, return single polygon
    # ------------------------------------------------------------------
    if idl_vertex_indices and not crossing_edge_indices:
        non_idl_lons = lons[~on_idl]
        eastern_count = int(np.sum((non_idl_lons > 0) & (non_idl_lons < 180.0)))
        western_count = int(np.sum((non_idl_lons < 0) & (non_idl_lons > -180.0)))

        coords_norm = aCoord_gcs.copy()
        for idx in idl_vertex_indices:
            if eastern_count >= western_count:
                coords_norm[idx, 0] = 180.0 - IDL_OFFSET
            else:
                coords_norm[idx, 0] = -180.0 + IDL_OFFSET

        if eastern_count >= western_count:
            return [coords_norm, np.array([])]
        else:
            return [np.array([]), coords_norm]

    # ------------------------------------------------------------------
    # 4. Validate crossing count
    # ------------------------------------------------------------------
    if len(crossing_edge_indices) != 2:
        raise ValueError(
            f"Expected exactly 2 IDL edge crossings, found {len(crossing_edge_indices)}. "
            "The polygon may not cross the IDL properly or may have a complex crossing pattern."
        )

    # ------------------------------------------------------------------
    # 5. Compute great-circle intersection latitudes for each crossing edge
    # ------------------------------------------------------------------
    EASTERN_BOUNDARY = 180.0 - IDL_OFFSET
    WESTERN_BOUNDARY = -180.0 + IDL_OFFSET

    def _intersection_lat(edge_idx: int) -> float:
        """Return the latitude at which edge *edge_idx* → *edge_idx+1* crosses ±180°."""
        p0 = aCoord_gcs[edge_idx]
        p1 = aCoord_gcs[(edge_idx + 1) % nPoint]
        lon0, lat0 = float(p0[0]), float(p0[1])
        lon1, lat1 = float(p1[0]), float(p1[1])
        # find_great_circle_intersection_with_meridian expects the positive-lon point first.
        if lon0 > 0:
            _, lat_cross = find_great_circle_intersection_with_meridian(
                lon0, lat0, lon1, lat1, 180.0
            )
        else:
            _, lat_cross = find_great_circle_intersection_with_meridian(
                lon1, lat1, lon0, lat0, 180.0
            )
        return lat_cross

    idx_A, idx_B = crossing_edge_indices  # A comes before B in CCW traversal order
    lat_A = _intersection_lat(idx_A)      # latitude where edge A crosses the IDL
    lat_B = _intersection_lat(idx_B)      # latitude where edge B crosses the IDL

    # ------------------------------------------------------------------
    # 6. Edge-by-edge traversal to build the two sub-polygons
    #
    # Walk every edge i → i+1 (skipping the redundant closing vertex).
    # At each IDL-crossing edge:
    #   a. Append the intersection point to the *current* sub-polygon
    #      (closing it at the IDL boundary).
    #   b. Switch the active sub-polygon.
    #   c. Open the new sub-polygon with the same intersection point
    #      (using the opposite boundary longitude).
    # Then add the next vertex to the now-active sub-polygon.
    #
    # This ties each intersection point to its specific crossing edge and
    # direction, producing geometrically correct results regardless of the
    # relative latitudes of the two crossings.
    # ------------------------------------------------------------------
    eastern_coords: List = []
    western_coords: List = []

    # Determine which hemisphere the first vertex belongs to.
    first_lon = float(lons[0])
    if first_lon >= 0:
        active, inactive = eastern_coords, western_coords
    else:
        active, inactive = western_coords, eastern_coords

    def _boundary_lon_for(coords_list: list) -> float:
        return EASTERN_BOUNDARY if coords_list is eastern_coords else WESTERN_BOUNDARY

    def _snap_and_append(coords_list: list, lon: float, lat: float) -> None:
        """Append vertex, snapping ±180° exactly to the boundary constant."""
        if abs(abs(lon) - 180.0) < IDL_TOLERANCE:
            lon = _boundary_lon_for(coords_list)
        coords_list.append([lon, lat])

    # Add the first vertex.
    _snap_and_append(active, float(aCoord_gcs[0, 0]), float(aCoord_gcs[0, 1]))

    crossing_lat_map = {idx_A: lat_A, idx_B: lat_B}

    # Iterate over edges 0→1, 1→2, …, (nPoint-2)→(nPoint-1).
    # The last vertex is the closure duplicate of vertex 0; we still add it so
    # the polygon ring is explicitly closed.
    for i in range(nPoint - 1):
        if i in crossing_lat_map:
            lat_cross = crossing_lat_map[i]
            # Close the current sub-polygon at the IDL.
            active.append([_boundary_lon_for(active), lat_cross])
            # Switch sub-polygons.
            active, inactive = inactive, active
            # Open the new sub-polygon at the IDL.
            active.append([_boundary_lon_for(active), lat_cross])

        # Add the destination vertex of this edge.
        next_idx = i + 1
        _snap_and_append(active, float(aCoord_gcs[next_idx, 0]), float(aCoord_gcs[next_idx, 1]))

    # ------------------------------------------------------------------
    # 7. Close both polygons and enforce CCW winding
    # ------------------------------------------------------------------
    def _close_polygon(coords_list: list) -> np.ndarray:
        """Convert to ndarray and close the ring if necessary."""
        if not coords_list:
            return np.array([])
        arr = np.array(coords_list, dtype=float)
        if not np.allclose(arr[0], arr[-1], atol=1e-10):
            arr = np.vstack([arr, arr[0:1]])
        return arr

    eastern_polygon = _close_polygon(eastern_coords)
    western_polygon = _close_polygon(western_coords)

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
    # Use module-level IDL_TOLERANCE and IDL_OFFSET constants
    idl_vertices = np.zeros_like(lons, dtype=bool)
    idl_vertices[:-1] = np.abs(np.abs(lons[:-1]) - 180.0) < IDL_TOLERANCE

    # Track hemisphere support using only non-IDL vertices.
    non_idl_lons = lons[~idl_vertices]
    has_eastern = np.any((non_idl_lons > 0) & (non_idl_lons < 180.0))
    has_western = np.any((non_idl_lons < 0) & (non_idl_lons > -180.0))
    spans_both_hemispheres = has_eastern and has_western

    idl_touch_positive = np.any(np.abs(lons[idl_vertices] - 180.0) < IDL_TOLERANCE)
    idl_touch_negative = np.any(np.abs(lons[idl_vertices] + 180.0) < IDL_TOLERANCE)
    touches_both_idl_sides = idl_touch_positive and idl_touch_negative

    # A polygon enclosing either pole necessarily crosses the IDL.
    if polygon_includes_pole(coords_updated, pole="north") or polygon_includes_pole(
        coords_updated, pole="south"
    ):
        return True, None

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

        # If no actual edge crossings but polygon spans both hemispheres, it crosses IDL
        if touches_both_idl_sides or spans_both_hemispheres:
            # This is a true IDL crossing with vertices on the meridian
            return True, None
        elif not (np.any(eastward_crossings) or np.any(westward_crossings)):
            # Adjust vertices that are exactly on the IDL by slightly moving them
            idl_indices = np.where(idl_vertices)[0]
            for idx in idl_indices:
                # Determine which hemisphere to move to based on neighboring vertices
                prev_idx = (idx - 1) % len(coords_updated)
                next_idx = (idx + 1) % len(coords_updated)

                prev_lon = coords_updated[prev_idx, 0]
                next_lon = coords_updated[next_idx, 0]

                # Count non-IDL neighbors to determine preferred hemisphere
                neighbor_lons = []
                if abs(abs(prev_lon) - 180.0) > IDL_TOLERANCE:  # Previous vertex not on IDL
                    neighbor_lons.append(prev_lon)
                if abs(abs(next_lon) - 180.0) > IDL_TOLERANCE:  # Next vertex not on IDL
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
    else:
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
