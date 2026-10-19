"""GeoTIFF raster reading utilities using GDAL."""

import os
from typing import Any, Dict, Optional, Union

import numpy as np
from osgeo import gdal

PathLike = Union[str, os.PathLike]

__all__ = [
    "gdal_read_geotiff_file",
    "gdal_read_geotiff_file_multiple_band",
]


def gdal_read_geotiff_file(
    filename: Optional[PathLike] = None,
    metadata_only: bool = False,
    **kwargs,
) -> Dict[str, Any]:
    """Read a single-band GeoTIFF raster and return data with metadata.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the GeoTIFF raster file.
    metadata_only : bool, default=False
        If True, return only metadata without reading raster data.

    Returns
    -------
    dict
        Dictionary containing raster metadata and optionally data:
        - 'pixelWidth', 'pixelHeight': pixel dimensions
        - 'originX', 'originY': origin coordinates
        - 'nrow', 'ncolumn': raster dimensions
        - 'geotransform': GDAL geotransform tuple
        - 'projection': spatial reference as WKT string
        - 'dataOut': raster data array (if metadata_only=False)
        - 'dataType': GDAL data type constant (if metadata_only=False)
        - 'missingValue': NoData value (if metadata_only=False)

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    RuntimeError
        If GDAL cannot open the file.
    ValueError
        If metadata is missing or invalid.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")
    if "iFlag_metadata_only" in kwargs:
        metadata_only = bool(kwargs["iFlag_metadata_only"])

    if filename is None:
        raise ValueError("Filename must be provided.")

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist.")

    dataset = gdal.Open(filepath, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError(f"Unable to open raster {filepath}.")

    band = None
    try:
        geotransform = dataset.GetGeoTransform()
        if not geotransform:
            raise ValueError(
                f"Raster {filepath} does not provide geotransform metadata."
            )

        origin_x = geotransform[0]
        origin_y = geotransform[3]
        pixel_width = geotransform[1]
        pixel_height = geotransform[5]

        ncolumn = dataset.RasterXSize
        nrow = dataset.RasterYSize

        projection_wkt = dataset.GetProjectionRef() or dataset.GetProjection()
        if not projection_wkt:
            raise ValueError(
                f"Raster {filepath} does not define spatial reference metadata."
            )

        if metadata_only:
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
            raise ValueError(f"Raster {filepath} does not contain a raster band.")

        data_type = band.DataType
        missing_value = band.GetNoDataValue()

        data = band.ReadAsArray(0, 0, ncolumn, nrow)
        if data is None:
            raise RuntimeError(f"Failed to read raster data from {filepath}.")

        return {
            "dataOut": data,
            "data": data,
            "dataType": data_type,
            "data_type": data_type,
            "pixelWidth": pixel_width,
            "pixel_width": pixel_width,
            "pixelHeight": pixel_height,
            "pixel_height": pixel_height,
            "originX": origin_x,
            "origin_x": origin_x,
            "originY": origin_y,
            "origin_y": origin_y,
            "nrow": nrow,
            "ncolumn": ncolumn,
            "missingValue": missing_value,
            "missing_value": missing_value,
            "geotransform": geotransform,
            "projection": projection_wkt,
            "projection_wkt": projection_wkt,
        }
    finally:
        band = None
        dataset = None


def gdal_read_geotiff_file_multiple_band(
    filename: Optional[PathLike] = None,
    metadata_only: bool = False,
    **kwargs,
) -> Dict[str, Any]:
    """Read a multi-band GeoTIFF raster and return stacked data with metadata.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the GeoTIFF raster file.
    metadata_only : bool, default=False
        If True, return only metadata without reading raster data.

    Returns
    -------
    dict
        Dictionary containing raster metadata and optionally data:
        - 'pixelWidth', 'pixelHeight': pixel dimensions
        - 'originX', 'originY': origin coordinates
        - 'nrow', 'ncolumn': raster dimensions
        - 'nband': number of bands
        - 'geotransform': GDAL geotransform tuple
        - 'projection': spatial reference as WKT string
        - 'dataOut': 3D array (nband, nrow, ncolumn) (if metadata_only=False)
        - 'dataType': GDAL data type constant (if metadata_only=False)
        - 'missingValue': NoData value (if metadata_only=False)
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")
    if "iFlag_metadata_only" in kwargs:
        metadata_only = bool(kwargs["iFlag_metadata_only"])

    if filename is None:
        raise ValueError("Filename must be provided.")

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist.")

    dataset = gdal.Open(filepath, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError(f"Unable to open raster {filepath}.")

    band = None
    try:
        geotransform = dataset.GetGeoTransform()
        if not geotransform:
            raise ValueError(
                f"Raster {filepath} does not provide geotransform metadata."
            )

        origin_x = geotransform[0]
        origin_y = geotransform[3]
        pixel_width = geotransform[1]
        pixel_height = geotransform[5]

        ncolumn = dataset.RasterXSize
        nrow = dataset.RasterYSize
        nband = dataset.RasterCount

        if nband < 1:
            raise ValueError(f"Raster {filepath} does not contain any bands.")

        projection_wkt = dataset.GetProjectionRef() or dataset.GetProjection()
        if not projection_wkt:
            raise ValueError(
                f"Raster {filepath} does not define spatial reference metadata."
            )

        if metadata_only:
            return {
                "pixelWidth": pixel_width,
                "pixelHeight": pixel_height,
                "originX": origin_x,
                "originY": origin_y,
                "nrow": nrow,
                "ncolumn": ncolumn,
                "nband": nband,
                "geotransform": geotransform,
                "projection": projection_wkt,
            }

        band = dataset.GetRasterBand(1)
        if band is None:
            raise ValueError(f"Raster {filepath} does not contain raster bands.")

        data_type = band.DataType
        missing_value = band.GetNoDataValue()

        first_data = band.ReadAsArray(0, 0, ncolumn, nrow)
        if first_data is None:
            raise RuntimeError(f"Failed to read raster data from {filepath}.")

        data = np.empty((nband, nrow, ncolumn), dtype=first_data.dtype)
        data[0, :, :] = first_data

        for i_band in range(1, nband):
            band = dataset.GetRasterBand(i_band + 1)
            if band is None:
                raise ValueError(f"Raster {filepath} is missing band {i_band + 1}.")

            band_data = band.ReadAsArray(0, 0, ncolumn, nrow)
            if band_data is None:
                raise RuntimeError(
                    f"Failed to read raster data from band {i_band + 1} of {filepath}."
                )

            data[i_band, :, :] = band_data

        return {
            "dataOut": data,
            "data": data,
            "dataType": data_type,
            "data_type": data_type,
            "pixelWidth": pixel_width,
            "pixel_width": pixel_width,
            "pixelHeight": pixel_height,
            "pixel_height": pixel_height,
            "originX": origin_x,
            "origin_x": origin_x,
            "originY": origin_y,
            "origin_y": origin_y,
            "nband": nband,
            "nrow": nrow,
            "ncolumn": ncolumn,
            "missingValue": missing_value,
            "missing_value": missing_value,
            "geotransform": geotransform,
            "projection": projection_wkt,
            "projection_wkt": projection_wkt,
        }
    finally:
        band = None
        dataset = None
