"""
Convert polygons crossing the International Date Line to valid geometries.

This module provides utilities for handling polygons that cross the International
Date Line (IDL, longitude ±180°). Such polygons often appear invalid in standard
geographic coordinate systems because they have longitude values that span more
than 180 degrees, causing the polygon to wrap around the globe in the wrong
direction.

The conversion process shifts longitudes in the western hemisphere (< -150°) to
the eastern hemisphere (by adding 360°), creating a continuous polygon in the
range [0°, 360°] instead of [-180°, 180°].

Functions
---------
convert_idl_polygon_to_valid_polygon
    Convert a polygon that crosses the International Date Line into a valid
    polygon by adjusting longitude values.

Notes
-----
The International Date Line roughly follows the 180° meridian but zigzags to
avoid splitting countries. Polygons crossing the IDL have a longitude range
greater than 180° and need special handling to be represented correctly.

This function uses a simple heuristic: if the polygon's longitude span exceeds
180°, it assumes the polygon crosses the IDL and shifts western longitudes
(< -150°) to the eastern hemisphere by adding 360°.

The threshold of -150° is chosen to handle typical IDL-crossing polygons near
the Pacific Ocean. Polygons in this region often have coordinates like:
    [..., (170, lat), (180, lat), (-170, lat), ...]

After conversion, these become:
    [..., (170, lat), (180, lat), (190, lat), ...]

Examples
--------
Convert a polygon crossing the IDL:

>>> from osgeo import ogr
>>> # Polygon crossing IDL: spans from 170°E to 170°W
>>> wkt = "POLYGON((170 10, 180 10, -170 10, -160 10, 170 10))"
>>> geom = ogr.CreateGeometryFromWkt(wkt)
>>> valid_geom = convert_idl_polygon_to_valid_polygon(geom)
>>> # Result has longitudes: 170, 180, 190, 200 (continuous)

Handle a normal polygon that doesn't cross IDL:

>>> # Normal polygon in Eastern hemisphere
>>> wkt = "POLYGON((100 10, 120 10, 120 20, 100 20, 100 10))"
>>> geom = ogr.CreateGeometryFromWkt(wkt)
>>> result = convert_idl_polygon_to_valid_polygon(geom)
>>> # Returns original geometry if valid and doesn't cross IDL

See Also
--------
pyearth.gis.geometry.split_polygon_cross_idl : Split IDL-crossing polygons
pyearth.gis.location.get_geometry_coordinates : Extract coordinates from geometry
"""

from typing import Optional
import numpy as np
from osgeo import ogr
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates


def convert_idl_polygon_to_valid_polygon(
    geometry_in: ogr.Geometry
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
    num_points = len(coordinates_gcs)
    adjusted_coords = coordinates_gcs.copy()

    for i in range(num_points):
        longitude = adjusted_coords[i, 0]
        if longitude < -150:
            # Shift to eastern hemisphere by adding 360°
            adjusted_coords[i, 0] = longitude + 360.0

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

    # Validate the output geometry
    if not geometry_out.IsValid():
        return None

    return geometry_out



