import os
from typing import Any, Dict, Optional

import numpy as np
from osgeo import gdal, osr


def gdal_read_geotiff_file(
    sFilename_in: str, iFlag_metadata_only: int = 0
) -> Dict[str, Any]:
    """Read a single-band GeoTIFF raster and return data with metadata.

    Parameters
    ----------
    sFilename_in : str
        Path to the GeoTIFF raster file.
    iFlag_metadata_only : int, optional
        If 1, return only metadata without reading raster data. Default is 0.

    Returns
    -------
    dict
        Dictionary containing raster metadata and optionally data. Keys include:
        - 'pixelWidth', 'pixelHeight': pixel dimensions
        - 'originX', 'originY': origin coordinates
        - 'nrow', 'ncolumn': raster dimensions
        - 'geotransform': GDAL geotransform tuple
        - 'projection': spatial reference as WKT string
        - 'dataOut': raster data array (if iFlag_metadata_only=0)
        - 'dataType': GDAL data type constant (if iFlag_metadata_only=0)
        - 'missingValue': NoData value (if iFlag_metadata_only=0)

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
        pixel_height = geotransform[5]

        ncolumn = dataset.RasterXSize
        nrow = dataset.RasterYSize
        nband = dataset.RasterCount

        projection_wkt = dataset.GetProjectionRef() or dataset.GetProjection()
        if not projection_wkt:
            raise ValueError(
                f"Raster {sFilename_in} does not define spatial reference metadata."
            )

        if iFlag_metadata_only == 1:
            return {
                "pixelWidth": pixel_width,
                "pixelHeight": pixel_height,
                "originX": origin_x,
                "originY": origin_y,
                "nrow": nrow,
                "ncolumn": ncolumn,
                "geotransform": geotransform,
                "projection": projection_wkt,
            }

        band = dataset.GetRasterBand(1)
        if band is None:
            raise ValueError(f"Raster {sFilename_in} does not contain a raster band.")

        data_type = band.DataType
        missing_value = band.GetNoDataValue()

        data = band.ReadAsArray(0, 0, ncolumn, nrow)
        if data is None:
            raise RuntimeError(f"Failed to read raster data from {sFilename_in}.")

        return {
            "dataOut": data,
            "dataType": data_type,
            "pixelWidth": pixel_width,
            "pixelHeight": pixel_height,
            "originX": origin_x,
            "originY": origin_y,
            "nrow": nrow,
            "ncolumn": ncolumn,
            "missingValue": missing_value,
            "geotransform": geotransform,
            "projection": projection_wkt,
        }
    finally:
        band = None
        dataset = None


def gdal_read_geotiff_file_multiple_band(
    sFilename_in: str, iFlag_metadata_only: int = 0
) -> Dict[str, Any]:
    """Read a multi-band GeoTIFF raster and return stacked data with metadata.

    Parameters
    ----------
    sFilename_in : str
        Path to the GeoTIFF raster file.
    iFlag_metadata_only : int, optional
        If 1, return only metadata without reading raster data. Default is 0.

    Returns
    -------
    dict
        Dictionary containing raster metadata and optionally data. Keys include:
        - 'pixelWidth', 'pixelHeight': pixel dimensions
        - 'originX', 'originY': origin coordinates
        - 'nrow', 'ncolumn': raster dimensions
        - 'nband': number of bands (if iFlag_metadata_only=0)
        - 'geotransform': GDAL geotransform tuple
        - 'projection': spatial reference as WKT string
        - 'dataOut': 3D array (nband, nrow, ncolumn) (if iFlag_metadata_only=0)
        - 'dataType': GDAL data type constant (if iFlag_metadata_only=0)
        - 'missingValue': NoData value (if iFlag_metadata_only=0)

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
        pixel_height = geotransform[5]

        ncolumn = dataset.RasterXSize
        nrow = dataset.RasterYSize
        nband = dataset.RasterCount

        if nband < 1:
            raise ValueError(f"Raster {sFilename_in} does not contain any bands.")

        projection_wkt = dataset.GetProjectionRef() or dataset.GetProjection()
        if not projection_wkt:
            raise ValueError(
                f"Raster {sFilename_in} does not define spatial reference metadata."
            )

        if iFlag_metadata_only == 1:
            return {
                "pixelWidth": pixel_width,
                "pixelHeight": pixel_height,
                "originX": origin_x,
                "originY": origin_y,
                "nrow": nrow,
                "ncolumn": ncolumn,
                "geotransform": geotransform,
                "projection": projection_wkt,
            }

        band = dataset.GetRasterBand(1)
        if band is None:
            raise ValueError(f"Raster {sFilename_in} does not contain raster bands.")

        data_type = band.DataType
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

        return {
            "dataOut": data,
            "dataType": data_type,
            "pixelWidth": pixel_width,
            "pixelHeight": pixel_height,
            "originX": origin_x,
            "originY": origin_y,
            "nband": nband,
            "nrow": nrow,
            "ncolumn": ncolumn,
            "missingValue": missing_value,
            "geotransform": geotransform,
            "projection": projection_wkt,
        }
    finally:
        band = None
        dataset = None
