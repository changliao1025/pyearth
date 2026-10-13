"""Raster spatial extent calculation utilities using GDAL."""

import os
from typing import Optional, Tuple, Union
from osgeo import gdal

from pyearth.gis.spatialref.is_wgs84_projection import is_wgs84_projection

PathLike = Union[str, os.PathLike]

__all__ = ["gdal_get_raster_extent"]


def gdal_get_raster_extent(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> Tuple[float, float, float, float]:
    """Compute the spatial extent of a raster dataset.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the raster file readable by GDAL.

    Returns
    -------
    tuple of float
        ``(min_x, max_x, min_y, max_y)`` in the dataset's coordinate system.
        For rasters in WGS84 (EPSG:4326), the result is clamped to global
        bounds ``[-180.0, 180.0]`` and ``[-90.0, 90.0]``.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    RuntimeError
        If GDAL cannot open the file.
    ValueError
        If the raster does not define a valid geotransform.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")

    if filename is None:
        raise ValueError("Filename must be provided.")

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist.")

    ds = gdal.Open(filepath, gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f"Unable to open raster {filepath}.")

    try:
        try:
            geotransform = ds.GetGeoTransform(can_return_null=True)
        except TypeError:
            geotransform = ds.GetGeoTransform()

        if not geotransform:
            raise ValueError(
                f"Raster {filepath} does not contain geotransform metadata."
            )

        raster_width = ds.RasterXSize
        raster_height = ds.RasterYSize

        corners = (
            gdal.ApplyGeoTransform(geotransform, 0, 0),
            gdal.ApplyGeoTransform(geotransform, raster_width, 0),
            gdal.ApplyGeoTransform(geotransform, raster_width, raster_height),
            gdal.ApplyGeoTransform(geotransform, 0, raster_height),
        )

        xs, ys = zip(*corners)
        min_x = min(xs)
        max_x = max(xs)
        min_y = min(ys)
        max_y = max(ys)

        projection = ds.GetProjection()
    finally:
        ds = None

    if projection and is_wgs84_projection(projection):
        min_x = max(min_x, -180.0)
        max_x = min(max_x, 180.0)
        min_y = max(min_y, -90.0)
        max_y = min(max_y, 90.0)

    return min_x, max_x, min_y, max_y
