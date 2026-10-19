"""GDAL raster format and driver support inspection utilities."""

import os
from functools import lru_cache
from typing import Dict, Union

from osgeo import gdal

from pyearth.system.filename import get_extension_from_path

PathLike = Union[str, os.PathLike]

SUPPORTED_RASTER_FORMATS: Dict[str, str] = {
    ".tif": "GTiff",
    ".tiff": "GTiff",
    ".jp2": "JP2OpenJPEG",
    ".j2k": "JP2OpenJPEG",
    ".img": "HFA",
    ".png": "PNG",
    ".jpg": "JPEG",
    ".jpeg": "JPEG",
    ".bmp": "BMP",
    ".gif": "GIF",
    ".nc": "netCDF",
    ".hdf": "HDF4",
    ".h5": "HDF5",
    ".hdf5": "HDF5",
    ".vrt": "VRT",
    ".asc": "AAIGrid",
    ".bil": "EHdr",
    ".bsq": "EHdr",
    ".bip": "EHdr",
    ".rst": "RST",
    ".grd": "GRD",
}

__all__ = [
    "SUPPORTED_RASTER_FORMATS",
    "gdal_raster_format_support",
    "get_available_raster_formats",
    "print_supported_raster_formats",
    "get_raster_format_from_extension",
    "get_raster_format_from_filename",
    "get_extension_from_raster_format",
    "get_raster_driver_from_format",
    "get_raster_driver_from_filename",
]


@lru_cache(maxsize=1)
def get_available_raster_formats() -> Dict[str, str]:
    """Get supported raster formats available in the active GDAL installation.

    Returns
    -------
    dict of str to str
        Mapping of file extensions to available GDAL driver names.
    """
    available_formats = {}
    for ext, driver_name in SUPPORTED_RASTER_FORMATS.items():
        if gdal.GetDriverByName(driver_name) is not None:
            available_formats[ext] = driver_name
    return available_formats


def gdal_raster_format_support() -> Dict[str, str]:
    """Return dictionary of supported raster formats available in GDAL.

    Returns
    -------
    dict of str to str
        Mapping of file extensions to available GDAL driver names.
    """
    return get_available_raster_formats()


def print_supported_raster_formats() -> None:
    """Print all supported and available raster formats in the current GDAL installation."""
    print("Supported raster formats (available in current GDAL installation):")
    available_formats = get_available_raster_formats()
    for ext, driver in available_formats.items():
        print(f"  {ext}: {driver}")


def get_raster_format_from_extension(extension: str) -> str:
    """Determine the GDAL raster format string from file extension.

    Parameters
    ----------
    extension : str
        File extension (e.g. '.tif').

    Returns
    -------
    str
        GDAL raster driver format name.

    Raises
    ------
    ValueError
        If the file extension is not supported.
    """
    ext = extension.strip().lower()
    if not ext.startswith("."):
        ext = "." + ext

    if ext in SUPPORTED_RASTER_FORMATS:
        return SUPPORTED_RASTER_FORMATS[ext]

    raise ValueError(f"Unsupported raster file extension: {ext}")


def get_raster_format_from_filename(filename: PathLike) -> str:
    """Determine the GDAL raster format string from a filename or path.

    Parameters
    ----------
    filename : str or os.PathLike
        Filename or path to inspect.

    Returns
    -------
    str
        GDAL raster driver format name.
    """
    if filename is None:
        raise ValueError("Filename cannot be None.")
    ext = get_extension_from_path(os.fspath(filename)).lower()
    if not ext:
        raise ValueError(
            f"Filename '{filename}' does not have an extension to determine raster format."
        )
    return get_raster_format_from_extension(ext)


def get_extension_from_raster_format(format_name: str) -> str:
    """Convert a GDAL raster format name to its standard file extension.

    Parameters
    ----------
    format_name : str
        GDAL driver name (e.g. 'GTiff', 'PNG', 'netCDF').

    Returns
    -------
    str
        Corresponding file extension (e.g. '.tif', '.png', '.nc').

    Raises
    ------
    ValueError
        If format_name is not found in supported raster formats.
    """
    for ext, fmt in SUPPORTED_RASTER_FORMATS.items():
        if fmt.lower() == format_name.strip().lower():
            return ext
    raise ValueError(f"Format '{format_name}' not found in supported raster formats.")


def get_raster_driver_from_format(file_format: str) -> gdal.Driver:
    """Get the GDAL driver object corresponding to a raster file format name.

    Parameters
    ----------
    file_format : str
        The raster file format name (e.g. 'GTiff', 'PNG').

    Returns
    -------
    gdal.Driver
        The GDAL driver object.

    Raises
    ------
    ValueError
        If the driver is not found.
    """
    driver = gdal.GetDriverByName(file_format.strip())
    if driver is None:
        raise ValueError(f"Driver for format '{file_format}' not found.")
    return driver


def get_raster_driver_from_filename(filename: PathLike) -> gdal.Driver:
    """Get the GDAL driver object corresponding to a filename.

    Parameters
    ----------
    filename : str or os.PathLike
        Filename or path to inspect.

    Returns
    -------
    gdal.Driver
        The GDAL driver object.
    """
    format_name = get_raster_format_from_filename(filename)
    return get_raster_driver_from_format(format_name)
