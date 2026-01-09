import os, sys
from typing import Any, Optional, Tuple
import numpy as np
from osgeo import gdal, osr


def gdal_read_envi_file(
    sFilename_in: str,
) -> Tuple[Any, float, float, float, int, int, Optional[float], Tuple[float, ...], str]:
    """Read a single-band ENVI raster and return data with metadata.

    Parameters
    ----------
    sFilename_in : str
        Path to the ENVI raster file.

    Returns
    -------
    tuple
        (data, pixel_width, origin_x, origin_y, nrow, ncolumn, missing_value, geotransform, projection_wkt)

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    RuntimeError
        If GDAL cannot open the file.
    ValueError
        If metadata is missing or invalid.
    """

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


def gdal_read_envi_file_multiple_band(
    sFilename_in: str,
) -> Tuple[
    np.ndarray,
    float,
    float,
    float,
    int,
    int,
    int,
    Optional[float],
    Tuple[float, ...],
    str,
]:
    """Read a multi-band ENVI raster and return stacked data with metadata.

    Parameters
    ----------
    sFilename_in : str
        Path to the ENVI raster file.

    Returns
    -------
    tuple
        (data, pixel_width, origin_x, origin_y, nband, nrow, ncolumn, missing_value, geotransform, projection_wkt)

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    RuntimeError
        If GDAL cannot open the file or read band data.
    ValueError
        If metadata is missing, invalid, or bands are missing.
    """

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
        nband = dataset.RasterCount

        if nband < 1:
            raise ValueError(f"Raster {sFilename_in} does not contain any bands.")

        band = dataset.GetRasterBand(1)
        if band is None:
            raise ValueError(f"Raster {sFilename_in} does not contain raster bands.")

        missing_value = band.GetNoDataValue()
        first_data = band.ReadAsArray(0, 0, ncolumn, nrow)
        if first_data is None:
            raise RuntimeError(f"Failed to read raster data from {sFilename_in}.")

        data = np.empty((nband, nrow, ncolumn), dtype=first_data.dtype)
        data[0, :, :] = first_data

        for i_band in range(1, nband):
            band = dataset.GetRasterBand(i_band + 1)
            if band is None:
                raise ValueError(f"Raster {sFilename_in} is missing band {i_band + 1}.")

            band_data = band.ReadAsArray(0, 0, ncolumn, nrow)
            if band_data is None:
                raise RuntimeError(
                    f"Failed to read raster data from band {i_band + 1} of {sFilename_in}."
                )

            data[i_band, :, :] = band_data

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
            nband,
            nrow,
            ncolumn,
            missing_value,
            geotransform,
            projection_wkt,
        )
    finally:
        band = None
        dataset = None
