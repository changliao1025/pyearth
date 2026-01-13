"""
Coordinate reprojection utilities using GDAL/OSR.

Functions for reprojecting coordinates between different spatial reference systems
using Well-Known Text (WKT) projection definitions.
"""

from typing import List, Tuple, Optional

import osgeo
from osgeo import osr


def reproject_coordinates(
    dx_in: float,
    dy_in: float,
    pProjection_source_in: str,
    pProjection_target_in: Optional[str] = None,
) -> Tuple[float, float]:
    """Reproject a single coordinate pair from source to target projection.

    Transforms coordinates between spatial reference systems using WKT
    (Well-Known Text) projection definitions. If no target projection is
    specified, coordinates are reprojected to WGS84 (EPSG:4326).

    Parameters
    ----------
    dx_in : float
        X coordinate (longitude or easting) in source projection.
    dy_in : float
        Y coordinate (latitude or northing) in source projection.
    pProjection_source_in : str
        WKT string defining the source spatial reference system.
    pProjection_target_in : str, optional
        WKT string defining the target spatial reference system.
        If None, defaults to WGS84 (EPSG:4326).

    Returns
    -------
    Tuple[float, float]
        Transformed coordinates (x_out, y_out) in target projection.

    Raises
    ------
    ValueError
        If WKT strings are invalid or coordinate transformation fails.
    RuntimeError
        If spatial reference import or transformation setup fails.

    Notes
    -----
    - Automatically handles GDAL 3+ axis order changes using traditional GIS order
    - Z coordinate is computed but not returned (use for 2D transformations)
    - For batch operations, use `reproject_coordinates_batch` for better performance

    Examples
    --------
    >>> # Reproject UTM to WGS84
    >>> utm_wkt = 'PROJCS["WGS 84 / UTM zone 10N",...]'
    >>> x_wgs, y_wgs = reproject_coordinates(500000, 4649776, utm_wkt)
    >>> # x_wgs ≈ -122.0, y_wgs ≈ 42.0

    See Also
    --------
    reproject_coordinates_batch : Reproject multiple coordinate pairs efficiently
    """
    try:
        # Validate inputs
        x = float(dx_in)
        y = float(dy_in)
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid coordinate values: x={dx_in}, y={dy_in}. Error: {e}")

    if not isinstance(pProjection_source_in, str) or not pProjection_source_in.strip():
        raise ValueError("Source projection must be a non-empty WKT string.")

    # Setup target spatial reference
    spatial_ref_target = osr.SpatialReference()
    try:
        if pProjection_target_in is not None:
            if (
                not isinstance(pProjection_target_in, str)
                or not pProjection_target_in.strip()
            ):
                raise ValueError("Target projection must be a non-empty WKT string.")
            result = spatial_ref_target.ImportFromWkt(pProjection_target_in)
            if result != 0:
                raise RuntimeError(
                    f"Failed to import target projection WKT. Error code: {result}"
                )
        else:
            # Default to WGS84
            result = spatial_ref_target.ImportFromEPSG(4326)
            if result != 0:
                raise RuntimeError(
                    f"Failed to import WGS84 (EPSG:4326). Error code: {result}"
                )
    except Exception as e:
        raise RuntimeError(f"Error setting up target spatial reference: {e}")

    # Setup source spatial reference
    spatial_ref_source = osr.SpatialReference()
    try:
        result = spatial_ref_source.ImportFromWkt(pProjection_source_in)
        if result != 0:
            raise RuntimeError(
                f"Failed to import source projection WKT. Error code: {result}"
            )
    except Exception as e:
        raise RuntimeError(f"Error importing source projection: {e}")

    # Handle GDAL 3+ axis order (https://github.com/OSGeo/gdal/issues/1546)
    if int(osgeo.__version__[0]) >= 3:
        spatial_ref_source.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        spatial_ref_target.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    # Create coordinate transformation
    try:
        transform = osr.CoordinateTransformation(spatial_ref_source, spatial_ref_target)
    except Exception as e:
        raise RuntimeError(f"Failed to create coordinate transformation: {e}")

    # Transform coordinates
    try:
        x_out, y_out, z = transform.TransformPoint(x, y)
    except Exception as e:
        raise ValueError(f"Coordinate transformation failed for point ({x}, {y}): {e}")

    return x_out, y_out


def reproject_coordinates_batch(
    aX_in: List[float],
    aY_in: List[float],
    pProjection_source_in: str,
    pProjection_target_in: Optional[str] = None,
) -> Tuple[List[float], List[float]]:
    """Reproject multiple coordinate pairs from source to target projection.

    Efficiently transforms a batch of coordinates between spatial reference systems
    using WKT (Well-Known Text) projection definitions. This function creates the
    transformation once and applies it to all coordinate pairs, making it more
    efficient than calling `reproject_coordinates` repeatedly.

    Parameters
    ----------
    aX_in : List[float]
        List of X coordinates (longitude or easting) in source projection.
    aY_in : List[float]
        List of Y coordinates (latitude or northing) in source projection.
        Must have the same length as aX_in.
    pProjection_source_in : str
        WKT string defining the source spatial reference system.
    pProjection_target_in : str, optional
        WKT string defining the target spatial reference system.
        If None, defaults to WGS84 (EPSG:4326).

    Returns
    -------
    Tuple[List[float], List[float]]
        Tuple of (x_out, y_out) lists containing transformed coordinates
        in target projection.

    Raises
    ------
    ValueError
        If input coordinate lists have different lengths, are empty,
        or if WKT strings are invalid.
    RuntimeError
        If spatial reference import or transformation setup fails.

    Notes
    -----
    - Automatically handles GDAL 3+ axis order changes using traditional GIS order
    - More efficient than calling `reproject_coordinates` in a loop
    - Returns lists in the same order as inputs
    - Z coordinates are computed but not returned

    Examples
    --------
    >>> # Reproject multiple UTM points to WGS84
    >>> utm_wkt = 'PROJCS["WGS 84 / UTM zone 10N",...]'
    >>> x_utm = [500000, 510000, 520000]
    >>> y_utm = [4649776, 4650000, 4650224]
    >>> x_wgs, y_wgs = reproject_coordinates_batch(x_utm, y_utm, utm_wkt)
    >>> # Returns lists of transformed coordinates

    See Also
    --------
    reproject_coordinates : Reproject a single coordinate pair
    """
    # Validate input lists
    if not isinstance(aX_in, (list, tuple)) or not isinstance(aY_in, (list, tuple)):
        raise ValueError(
            "Input coordinates must be lists or tuples. "
            f"Got types: x={type(aX_in).__name__}, y={type(aY_in).__name__}"
        )

    n_points = len(aX_in)
    if n_points == 0:
        raise ValueError("Input coordinate lists cannot be empty.")

    if len(aY_in) != n_points:
        raise ValueError(
            f"Coordinate lists must have the same length. "
            f"Got x={n_points}, y={len(aY_in)}"
        )

    if not isinstance(pProjection_source_in, str) or not pProjection_source_in.strip():
        raise ValueError("Source projection must be a non-empty WKT string.")

    # Setup target spatial reference
    spatial_ref_target = osr.SpatialReference()
    try:
        if pProjection_target_in is not None:
            if (
                not isinstance(pProjection_target_in, str)
                or not pProjection_target_in.strip()
            ):
                raise ValueError("Target projection must be a non-empty WKT string.")
            result = spatial_ref_target.ImportFromWkt(pProjection_target_in)
            if result != 0:
                raise RuntimeError(
                    f"Failed to import target projection WKT. Error code: {result}"
                )
        else:
            # Default to WGS84
            result = spatial_ref_target.ImportFromEPSG(4326)
            if result != 0:
                raise RuntimeError(
                    f"Failed to import WGS84 (EPSG:4326). Error code: {result}"
                )
    except Exception as e:
        raise RuntimeError(f"Error setting up target spatial reference: {e}")

    # Setup source spatial reference
    spatial_ref_source = osr.SpatialReference()
    try:
        result = spatial_ref_source.ImportFromWkt(pProjection_source_in)
        if result != 0:
            raise RuntimeError(
                f"Failed to import source projection WKT. Error code: {result}"
            )
    except Exception as e:
        raise RuntimeError(f"Error importing source projection: {e}")

    # Handle GDAL 3+ axis order (https://github.com/OSGeo/gdal/issues/1546)
    if int(osgeo.__version__[0]) >= 3:
        spatial_ref_source.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        spatial_ref_target.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    # Create coordinate transformation
    try:
        transform = osr.CoordinateTransformation(spatial_ref_source, spatial_ref_target)
    except Exception as e:
        raise RuntimeError(f"Failed to create coordinate transformation: {e}")

    # Transform all coordinates
    aX_out = []
    aY_out = []

    for i in range(n_points):
        try:
            x = float(aX_in[i])
            y = float(aY_in[i])
        except (TypeError, ValueError) as e:
            raise ValueError(
                f"Invalid coordinate at index {i}: x={aX_in[i]}, y={aY_in[i]}. Error: {e}"
            )

        try:
            x_out, y_out, z = transform.TransformPoint(x, y)
            aX_out.append(x_out)
            aY_out.append(y_out)
        except Exception as e:
            raise ValueError(
                f"Coordinate transformation failed at index {i} for point ({x}, {y}): {e}"
            )

    return aX_out, aY_out
