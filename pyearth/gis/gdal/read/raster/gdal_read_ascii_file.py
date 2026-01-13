import os
from typing import Any, Optional, Tuple
from osgeo import gdal


def gdal_read_ascii_file(
    sFilename_in: str,
) -> Tuple[Any, float, float, float, int, int, Optional[float], Tuple[float, ...], str]:
    """Read an ASCII/AAIGrid raster file and return data and metadata."""

    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"File {sFilename_in} does not exist.")

    dataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError(f"Unable to open raster {sFilename_in}.")

    band = None
    try:
        geotransform = dataset.GetGeoTransform()
        if not geotransform:
            raise ValueError(
                f"Raster {sFilename_in} does not provide geotransform metadata."
            )

        origin_x = geotransform[0]
        origin_y = geotransform[3]
        pixel_width = geotransform[1]

        ncolumn = dataset.RasterXSize
        nrow = dataset.RasterYSize

        band = dataset.GetRasterBand(1)
        if band is None:
            raise ValueError(f"Raster {sFilename_in} does not contain a raster band.")

        missing_value = band.GetNoDataValue()
        data = band.ReadAsArray(0, 0, ncolumn, nrow)
        if data is None:
            raise RuntimeError(f"Failed to read raster data from {sFilename_in}.")

        projection_wkt = dataset.GetProjectionRef() or dataset.GetProjection()
        if not projection_wkt:
            raise ValueError(
                f"Raster {sFilename_in} does not define spatial reference metadata."
            )

        return (
            data,
            pixel_width,
            origin_x,
            origin_y,
            nrow,
            ncolumn,
            missing_value,
            geotransform,
            projection_wkt,
        )
    finally:
        band = None
        dataset = None
