import os
from typing import Tuple
from osgeo import gdal
from pyearth.gis.spatialref.is_wgs84_projection import is_wgs84_projection


def gdal_get_raster_extent(sFilename_in: str) -> Tuple[float, float, float, float]:
    """
    Compute the spatial extent of a raster dataset.

    Parameters
    ----------
    sFilename_in : str
        Path to the raster file readable by GDAL.

    Returns
    -------
    tuple of float
        ``(min_x, max_x, min_y, max_y)`` in the dataset's coordinate system.
        For rasters in EPSG:4326, the result is clamped to the global longitude
        and latitude bounds of ``[-180, 180]`` and ``[-90, 90]`` respectively.

    Raises
    ------
    FileNotFoundError
        If ``sFilename_in`` does not exist.
    RuntimeError
        If GDAL cannot open the file.
    ValueError
        If the raster does not define a valid geotransform.
    """

    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"File {sFilename_in} does not exist.")

    ds = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f"Unable to open raster {sFilename_in}.")

    try:
        try:
            geotransform = ds.GetGeoTransform(can_return_null=True)
        except TypeError:
            geotransform = ds.GetGeoTransform()

        if not geotransform:
            raise ValueError(
                f"Raster {sFilename_in} does not contain geotransform metadata."
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

    if is_wgs84_projection(projection):
        min_x = max(min_x, -180.0)
        max_x = min(max_x, 180.0)
        min_y = max(min_y, -90.0)
        max_y = min(max_y, 90.0)

    return min_x, max_x, min_y, max_y
