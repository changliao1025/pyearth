"""
Clean and repair OGR/GDAL geometry objects.

This module provides utilities for cleaning and repairing geometry objects that may
contain artifacts, self-intersections, duplicate vertices, or other topological
issues. These issues commonly occur after spatial operations like union, intersection,
or buffer operations.

The cleaning process is multi-step:
    1. Fix self-intersections and invalid topology (MakeValid or Buffer(0))
    2. Simplify to remove collinear points and tiny segments
    3. Remove slivers using negative/positive buffer operations
    4. Validate and fix any remaining issues

Functions
---------
clean_geometry
    Clean and repair an OGR geometry object by removing artifacts and fixing
    topological issues.

Notes
-----
The cleaning operations may slightly modify the geometry coordinates due to
buffering and simplification. The tolerance parameter controls the magnitude
of these modifications.

For GDAL 3.0+, the MakeValid() method is preferred as it uses more sophisticated
algorithms (often based on the GEOS MakeValid function). For older versions,
the Buffer(0) trick is used as a fallback.

Examples
--------
Clean a geometry after a union operation:

>>> from osgeo import ogr
>>> # Assume geom is an OGR geometry that resulted from a union
>>> clean_geom = clean_geometry(geom, tolerance=1e-6, verbose=True)

Clean a self-intersecting polygon:

>>> # Create a self-intersecting polygon (bowtie shape)
>>> wkt = "POLYGON((0 0, 2 2, 2 0, 0 2, 0 0))"
>>> geom = ogr.CreateGeometryFromWkt(wkt)
>>> clean_geom = clean_geometry(geom, tolerance=1e-8)
>>> print(clean_geom.IsValid())
True

See Also
--------
osgeo.ogr.Geometry.MakeValid : GDAL 3.0+ method for fixing invalid geometries
osgeo.ogr.Geometry.Buffer : Buffer operation used for cleaning
osgeo.ogr.Geometry.Simplify : Simplification to remove unnecessary vertices
"""

from typing import Optional
from osgeo import ogr


def clean_geometry(
    geometry: Optional[ogr.Geometry], tolerance: float = 1e-6, verbose: bool = False
) -> Optional[ogr.Geometry]:
    """
    Clean geometry artifacts and issues that can occur after spatial operations.

    This function applies a series of cleaning operations to fix common geometry
    problems including self-intersections, duplicate vertices, collinear points,
    and sliver polygons. The cleaning process is designed to be robust and will
    attempt multiple strategies to produce a valid geometry.

    Parameters
    ----------
    geometry : ogr.Geometry or None
        The OGR geometry object to clean. Can be any geometry type (Point,
        LineString, Polygon, MultiPolygon, etc.). If None, returns None.
    tolerance : float, optional
        Tolerance value (in geometry units) for simplification and buffer
        operations. Larger values result in more aggressive simplification.
        Must be non-negative. Default is 1e-6.
        - Used for Simplify() to remove vertices within this distance
        - Used to calculate buffer tolerance (tolerance * 10) for sliver removal
    verbose : bool, optional
        If True, print warning messages about geometry type changes and
        cleaning failures. Default is False.

    Returns
    -------
    ogr.Geometry or None
        The cleaned geometry object, or None if the input was None. If cleaning
        completely fails, returns the original geometry unchanged.

    Raises
    ------
    ValueError
        If tolerance is negative.

    Notes
    -----
    The cleaning process follows these steps:

    1. **Fix Self-Intersections**: Uses MakeValid() (GDAL 3.0+) or Buffer(0) to
       fix self-intersections and invalid topology. The Buffer(0) trick works
       because buffering by zero distance forces the geometry to be reconstructed
       using valid topology rules.

    2. **Simplify**: Removes collinear points and tiny segments using the
       Douglas-Peucker algorithm (via Simplify()). This reduces the number of
       vertices while maintaining the overall shape within the tolerance.

    3. **Remove Slivers**: Uses a negative buffer followed by a positive buffer
       to remove very thin areas (slivers) that can appear as lines. The buffer
       distance is tolerance * 10 (or 1e-8 if tolerance is 0). This operation
       may remove small features smaller than twice the buffer distance.

    4. **Validate Geometry Type**: Checks if the geometry type changed during
       cleaning (e.g., a polygon became a multipolygon due to splitting). This
       is logged if verbose=True but the cleaned geometry is still returned.

    5. **Final Validation**: If the geometry is still invalid after all steps,
       applies a final Buffer(0) operation as a last resort.

    The function is designed to be robust and will not raise exceptions during
    normal operation. If any step fails, it will either skip that step or return
    the geometry in its current state.

    Warnings
    --------
    - The cleaning process may change the geometry type (e.g., Polygon to
      MultiPolygon if the cleaning splits the geometry).
    - Small features smaller than 2 * (tolerance * 10) may be removed during
      sliver removal.
    - Coordinates may be slightly modified due to buffering and simplification.
    - For GDAL < 3.0, only the Buffer(0) method is available for fixing
      self-intersections, which may not handle all invalid geometry cases.

    Examples
    --------
    Clean a geometry with default settings:

    >>> from osgeo import ogr
    >>> wkt = "POLYGON((0 0, 1 1, 1 0, 0 1, 0 0))"  # Self-intersecting
    >>> geom = ogr.CreateGeometryFromWkt(wkt)
    >>> clean_geom = clean_geometry(geom)
    >>> print(clean_geom.IsValid())
    True

    Clean with custom tolerance and verbose output:

    >>> geom = ogr.CreateGeometryFromWkt("POLYGON((0 0, 1 0, 1 1, 0.5 1, 0 1, 0 0))")
    >>> clean_geom = clean_geometry(geom, tolerance=0.1, verbose=True)

    Handle None input gracefully:

    >>> result = clean_geometry(None)
    >>> print(result)
    None

    Clean a geometry after union operation:

    >>> poly1 = ogr.CreateGeometryFromWkt("POLYGON((0 0, 2 0, 2 2, 0 2, 0 0))")
    >>> poly2 = ogr.CreateGeometryFromWkt("POLYGON((1 1, 3 1, 3 3, 1 3, 1 1))")
    >>> union = poly1.Union(poly2)
    >>> clean_union = clean_geometry(union, tolerance=1e-8)

    See Also
    --------
    osgeo.ogr.Geometry.MakeValid : GDAL 3.0+ method for fixing invalid geometries
    osgeo.ogr.Geometry.Buffer : Buffer operation
    osgeo.ogr.Geometry.Simplify : Douglas-Peucker simplification
    osgeo.ogr.Geometry.IsValid : Check if geometry is topologically valid
    """
    # Validate inputs
    if tolerance < 0:
        raise ValueError(f"tolerance must be non-negative, got {tolerance}")

    if geometry is None:
        return None

    try:
        original_type = geometry.GetGeometryName()

        # Step 1: Remove duplicate vertices and fix self-intersections
        if hasattr(geometry, "MakeValid"):
            # Use MakeValid if available (GDAL 3.0+)
            # This is the preferred method as it handles complex cases better
            geometry = geometry.MakeValid()
        else:
            # Fallback: use buffer trick to fix self-intersections
            # Buffer(0) forces reconstruction with valid topology
            geometry = geometry.Buffer(0)

        # Check if geometry became None after MakeValid/Buffer
        if geometry is None:
            if verbose:
                print("Warning: Geometry became None after MakeValid/Buffer operation")
            return None

        # Step 2: Simplify to remove collinear points and tiny segments
        if tolerance > 0:
            simplified = geometry.Simplify(tolerance)
            if simplified is not None:
                geometry = simplified
            elif verbose:
                print(
                    "Warning: Simplify operation returned None, keeping previous geometry"
                )

        # Step 3: Remove slivers using small negative then positive buffer
        # This removes very thin areas that might appear as lines
        buffer_tolerance = tolerance * 10 if tolerance > 0 else 1e-8
        try:
            # Negative buffer removes slivers
            temp_geom = geometry.Buffer(-buffer_tolerance)
            if temp_geom is not None and not temp_geom.IsEmpty():
                # Positive buffer restores to approximate original size
                buffered = temp_geom.Buffer(buffer_tolerance)
                if buffered is not None:
                    geometry = buffered
                elif verbose:
                    print(
                        "Warning: Positive buffer returned None, keeping previous geometry"
                    )
            elif verbose:
                print(
                    "Warning: Negative buffer resulted in empty geometry, skipping sliver removal"
                )
        except Exception as e:
            if verbose:
                print(
                    f"Warning: Buffer operation failed: {str(e)}, skipping sliver removal"
                )

        # Check if geometry became None after buffering
        if geometry is None:
            if verbose:
                print("Warning: Geometry became None after buffer operations")
            return None

        # Step 4: Ensure we still have the correct geometry type
        current_type = geometry.GetGeometryName()
        if current_type != original_type:
            if verbose:
                print(
                    f"Warning: Geometry type changed from {original_type} to {current_type} during cleaning"
                )

        # Step 5: Final validation
        if not geometry.IsValid():
            if verbose:
                print(
                    "Warning: Geometry still invalid after cleaning, applying final buffer fix"
                )
            final_fix = geometry.Buffer(0)
            if final_fix is not None:
                geometry = final_fix
            elif verbose:
                print(
                    "Warning: Final buffer fix returned None, returning previous geometry"
                )

        return geometry

    except Exception as e:
        if verbose:
            print(f"Warning: Failed to clean geometry: {str(e)}")
        return geometry  # Return original if cleaning fails
