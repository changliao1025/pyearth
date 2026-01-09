import os
from typing import Literal

from osgeo import gdal, ogr


def gdal_check_file_type(sFilename_in: str) -> Literal["raster", "vector", "unknown"]:
    """Determine whether a file is a raster, vector, or unknown spatial data type.

    Parameters
    ----------
    filename : str
        Path to the file to check.

    Returns
    -------
    str
        One of 'raster', 'vector', or 'unknown'.
        - 'raster': File can be opened as a GDAL raster dataset.
        - 'vector': File can be opened as an OGR vector dataset.
        - 'unknown': File does not exist or is not a recognizable GDAL/OGR format.

    Notes
    -----
    This function attempts to open the file with GDAL first (raster), then OGR (vector).
    Some formats may be readable by both drivers; raster takes precedence.
    """

    if not os.path.exists(sFilename_in):
        return "unknown"

    # Try to open as raster
    raster = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    if raster is not None:
        raster = None  # Close dataset
        return "raster"

    # Try to open as vector
    vector = ogr.Open(sFilename_in, 0)  # 0 = read-only
    if vector is not None:
        vector = None  # Close dataset
        return "vector"

    return "unknown"
