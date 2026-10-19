"""Coordinate reprojection utilities using GDAL/OSR.

Functions for reprojecting coordinates between different spatial reference systems
using Well-Known Text (WKT), EPSG codes, or SpatialReference objects.
"""

from typing import List, Optional, Sequence, Tuple, Union

from osgeo import osr

osr.UseExceptions()

__all__ = [
    "reproject_coordinates",
    "reproject_coordinates_batch",
]

SpatialRefLike = Union[str, osr.SpatialReference, int]


def _build_spatial_reference(
    projection: Optional[SpatialRefLike], is_target: bool = False
) -> osr.SpatialReference:
    """Construct an osr.SpatialReference from WKT, EPSG code, or existing SpatialReference."""
    if projection is None:
        if is_target:
            sr = osr.SpatialReference()
            sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            sr.ImportFromEPSG(4326)
            return sr
        raise ValueError("Source projection cannot be None.")

    sr = osr.SpatialReference()
    sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    if isinstance(projection, osr.SpatialReference):
        cloned = projection.Clone()
        cloned.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        return cloned

    if isinstance(projection, int):
        try:
            sr.ImportFromEPSG(projection)
            return sr
        except Exception as e:
            raise RuntimeError(f"Failed to import EPSG:{projection}: {e}") from e

    if isinstance(projection, str):
        proj_str = projection.strip()
        if not proj_str:
            raise ValueError("Projection string cannot be empty.")

        # Check if it's an EPSG string
        proj_upper = proj_str.upper()
        if proj_upper.startswith("EPSG:"):
            try:
                code = int(proj_upper.split(":")[1].strip())
                sr.ImportFromEPSG(code)
                return sr
            except (ValueError, IndexError):
                pass
        elif proj_str.isdigit():
            try:
                sr.ImportFromEPSG(int(proj_str))
                return sr
            except Exception:
                pass

        # Try WKT import
        try:
            if sr.ImportFromWkt(proj_str) == 0:
                return sr
        except Exception:
            pass

        # Try user input (handles PROJ, URN, etc.)
        try:
            if sr.SetFromUserInput(proj_str) == 0:
                return sr
        except Exception as e:
            raise RuntimeError(
                f"Failed to parse projection '{proj_str[:50]}...': {e}"
            ) from e

    raise ValueError(f"Unsupported projection type: {type(projection).__name__}")


def reproject_coordinates(
    x: Optional[float] = None,
    y: Optional[float] = None,
    source_projection: Optional[SpatialRefLike] = None,
    target_projection: Optional[SpatialRefLike] = None,
    **kwargs,
) -> Tuple[float, float]:
    """Reproject a single coordinate pair from source to target projection.

    Transforms coordinates between spatial reference systems. If no target
    projection is specified, coordinates are reprojected to WGS84 (EPSG:4326).

    Parameters
    ----------
    x : float
        X coordinate (longitude or easting) in source projection.
    y : float
        Y coordinate (latitude or northing) in source projection.
    source_projection : str, int, or osr.SpatialReference
        Definition of the source spatial reference system (WKT, EPSG, etc.).
    target_projection : str, int, or osr.SpatialReference, optional
        Definition of the target spatial reference system.
        If None, defaults to WGS84 (EPSG:4326).

    Returns
    -------
    Tuple[float, float]
        Transformed coordinates (x_out, y_out) in target projection.

    Raises
    ------
    ValueError
        If coordinates are invalid or transformation fails.
    RuntimeError
        If spatial reference setup fails.
    """
    if x is None:
        x = kwargs.get("dx_in")
    if y is None:
        y = kwargs.get("dy_in")
    if source_projection is None:
        source_projection = kwargs.get("pProjection_source_in")
    if target_projection is None and "pProjection_target_in" in kwargs:
        target_projection = kwargs.get("pProjection_target_in")

    if x is None or y is None:
        raise ValueError("Both x and y coordinates must be provided.")

    try:
        x_val = float(x)
        y_val = float(y)
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid coordinate values: x={x}, y={y}. Error: {e}") from e

    sr_source = _build_spatial_reference(source_projection, is_target=False)
    sr_target = _build_spatial_reference(target_projection, is_target=True)

    try:
        transform = osr.CoordinateTransformation(sr_source, sr_target)
    except Exception as e:
        raise RuntimeError(f"Failed to create coordinate transformation: {e}") from e

    try:
        x_out, y_out, _ = transform.TransformPoint(x_val, y_val)
    except Exception as e:
        raise ValueError(
            f"Coordinate transformation failed for point ({x_val}, {y_val}): {e}"
        ) from e

    return float(x_out), float(y_out)


def reproject_coordinates_batch(
    x_coordinates: Optional[Sequence[float]] = None,
    y_coordinates: Optional[Sequence[float]] = None,
    source_projection: Optional[SpatialRefLike] = None,
    target_projection: Optional[SpatialRefLike] = None,
    **kwargs,
) -> Tuple[List[float], List[float]]:
    """Reproject multiple coordinate pairs from source to target projection.

    Efficiently transforms a batch of coordinates between spatial reference systems
    using GDAL's vectorized C++ transformation pipeline.

    Parameters
    ----------
    x_coordinates : Sequence[float]
        X coordinates (longitudes or eastings) in source projection.
    y_coordinates : Sequence[float]
        Y coordinates (latitudes or northings) in source projection.
        Must have the same length as x_coordinates.
    source_projection : str, int, or osr.SpatialReference
        Definition of the source spatial reference system.
    target_projection : str, int, or osr.SpatialReference, optional
        Definition of the target spatial reference system.
        If None, defaults to WGS84 (EPSG:4326).

    Returns
    -------
    Tuple[List[float], List[float]]
        Tuple of (x_out, y_out) lists containing transformed coordinates.

    Raises
    ------
    ValueError
        If coordinate sequences have mismatched lengths or are empty.
    RuntimeError
        If spatial reference setup fails.
    """
    if x_coordinates is None:
        x_coordinates = kwargs.get("aX_in")
    if y_coordinates is None:
        y_coordinates = kwargs.get("aY_in")
    if source_projection is None:
        source_projection = kwargs.get("pProjection_source_in")
    if target_projection is None and "pProjection_target_in" in kwargs:
        target_projection = kwargs.get("pProjection_target_in")

    if x_coordinates is None or y_coordinates is None:
        raise ValueError("Both x_coordinates and y_coordinates must be provided.")

    n_points = len(x_coordinates)
    if n_points == 0:
        raise ValueError("Input coordinate sequences cannot be empty.")
    if len(y_coordinates) != n_points:
        raise ValueError(
            f"Coordinate sequences must have the same length. "
            f"Got x={n_points}, y={len(y_coordinates)}"
        )

    try:
        pts = [(float(x), float(y)) for x, y in zip(x_coordinates, y_coordinates)]
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid coordinate in batch: {e}") from e

    sr_source = _build_spatial_reference(source_projection, is_target=False)
    sr_target = _build_spatial_reference(target_projection, is_target=True)

    try:
        transform = osr.CoordinateTransformation(sr_source, sr_target)
    except Exception as e:
        raise RuntimeError(f"Failed to create coordinate transformation: {e}") from e

    try:
        transformed = transform.TransformPoints(pts)
    except Exception as e:
        raise ValueError(f"Batch coordinate transformation failed: {e}") from e

    x_out = [pt[0] for pt in transformed]
    y_out = [pt[1] for pt in transformed]

    return x_out, y_out
