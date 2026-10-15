"""Spatial dataset type identification utilities."""

import os
from typing import Literal, Optional, Union

from osgeo import gdal, ogr

gdal.UseExceptions()
ogr.UseExceptions()

PathLike = Union[str, os.PathLike]

__all__ = ["gdal_check_file_type"]


def gdal_check_file_type(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> Literal["raster", "vector", "unknown"]:
    """Determine whether a file is a raster, vector, or unknown spatial dataset.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the file to check.

    Returns
    -------
    str
        One of 'raster', 'vector', or 'unknown'.
        - 'raster': File can be opened as a GDAL raster dataset.
        - 'vector': File can be opened as an OGR vector dataset.
        - 'unknown': File does not exist or is not recognized.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")

    if filename is None:
        return "unknown"

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        return "unknown"

    gdal.PushErrorHandler("CPLQuietErrorHandler")
    try:
        # Try opening as raster
        try:
            raster = gdal.Open(filepath, gdal.GA_ReadOnly)
            if raster is not None:
                raster = None
                return "raster"
        except Exception:
            pass

        # Try opening as vector
        try:
            vector = ogr.Open(filepath, 0)
            if vector is not None:
                vector = None
                return "vector"
        except Exception:
            pass
    finally:
        gdal.PopErrorHandler()

    return "unknown"


