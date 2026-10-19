"""ASCII raster reading utilities using GDAL."""

import os
from typing import Optional, Union

from osgeo import gdal

from pyearth.gis.gdal.read.raster.gdal_read_envi_file import RasterDataset

gdal.UseExceptions()

PathLike = Union[str, os.PathLike]

__all__ = ["gdal_read_ascii_file"]


def gdal_read_ascii_file(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> RasterDataset:
    """Read an ASCII/AAIGrid raster file and return data and metadata.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the ASCII raster file.

    Returns
    -------
    RasterDataset
        Container unpackable as a 9-tuple:
        (data, pixel_width, origin_x, origin_y, nrow, ncolumn, missing_value, geotransform, projection_wkt)
        and accessible via dictionary keys or attributes.

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

        ncolumn = dataset.RasterXSize
        nrow = dataset.RasterYSize

        band = dataset.GetRasterBand(1)
        if band is None:
            raise ValueError(f"Raster {filepath} does not contain a raster band.")

        missing_value = band.GetNoDataValue()
        data = band.ReadAsArray(0, 0, ncolumn, nrow)
        if data is None:
            raise RuntimeError(f"Failed to read raster data from {filepath}.")

        projection_wkt = dataset.GetProjectionRef() or dataset.GetProjection() or ""

        return RasterDataset(
            data=data,
            pixel_width=pixel_width,
            origin_x=origin_x,
            origin_y=origin_y,
            nrow=nrow,
            ncolumn=ncolumn,
            missing_value=missing_value,
            geotransform=geotransform,
            projection_wkt=projection_wkt,
        )
    finally:
        band = None
        dataset = None
