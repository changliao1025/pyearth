import os
from functools import lru_cache
from typing import Optional
from osgeo import ogr

from pyearth.system.filename import get_extension_from_path

SUPPORTED_VECTOR_FORMATS = {
    ".geojson": "GeoJSON",
    ".json": "GeoJSON",
    ".shp": "ESRI Shapefile",
    ".gpkg": "GPKG",
    ".kml": "KML",
    ".gml": "GML",
    ".sqlite": "SQLite",
    ".parquet": "Parquet",
    ".geoparquet": "Parquet",
    ".csv": "CSV",
    ".tab": "MapInfo File",
    ".mif": "MapInfo File",
    ".dxf": "DXF",
    ".gdb": "OpenFileGDB",
}


__all__ = [
    "SUPPORTED_VECTOR_FORMATS",
    "gdal_vector_format_support",
    "get_available_vector_formats",
    "print_supported_vector_formats",
    "get_vector_format_from_extension",
    "get_vector_driver_from_filename",
    "get_vector_driver_from_format",
    "get_vector_format_from_filename",
    "check_parquet_support",
    "has_parquet_support",
]


@lru_cache(maxsize=1)
def get_available_vector_formats() -> dict[str, str]:
    """
    Get supported vector formats that are actually available in the current GDAL/OGR installation.

    Returns:
    dict[str, str]: A dictionary of available vector formats mapping file extensions to OGR driver names.
    """
    available_formats = {}

    for ext, format_name in SUPPORTED_VECTOR_FORMATS.items():
        # Special handling for Parquet/GeoParquet
        if ext in (".parquet", ".geoparquet"):
            if check_parquet_support() is not None:
                available_formats[ext] = "Parquet"  # Use the canonical name
        else:
            # Check if the driver is available
            driver = ogr.GetDriverByName(format_name)
            if driver is not None:
                available_formats[ext] = format_name

    return available_formats


def gdal_vector_format_support() -> dict[str, str]:
    """
    Return a dictionary of supported vector formats that are available in the current GDAL/OGR installation.

    This function filters the static SUPPORTED_VECTOR_FORMATS to only include formats
    where the corresponding OGR driver is actually available at runtime.

    Returns:
    dict[str, str]: A dictionary of available vector formats mapping file extensions to OGR driver names.
    """
    return get_available_vector_formats()


def print_supported_vector_formats() -> None:
    """
    Print all supported vector formats that are available in the current GDAL/OGR installation.
    """
    print("Supported vector formats (available in current GDAL/OGR installation):")
    available_formats = get_available_vector_formats()
    for ext, driver in available_formats.items():
        print(f"  {ext}: {driver}")

    # Also show unavailable formats from the supported list
    unavailable_formats = {}
    for ext, format_name in SUPPORTED_VECTOR_FORMATS.items():
        if ext not in available_formats:
            unavailable_formats[ext] = format_name

    if unavailable_formats:
        print("\nFormats in supported list but not available in current installation:")
        for ext, driver in unavailable_formats.items():
            print(f"  {ext}: {driver} (driver not found)")


def get_vector_format_from_extension(sExtension: str) -> str:
    """
    Determine the OGR format string from file extension.

    Only accepts extensions that are both in the supported list and have available drivers.

    Parameters:
    sExtension (str): Input file extension (e.g., '.shp').

    Returns:
    str: Format string for OGR driver.

    Raises:
    ValueError: If the file extension is not supported or driver is not available.
    """
    ext = sExtension.lower()

    available_formats = get_available_vector_formats()
    if ext not in available_formats:
        if ext in SUPPORTED_VECTOR_FORMATS:
            raise ValueError(
                f"File extension '{ext}' is supported but driver '{SUPPORTED_VECTOR_FORMATS[ext]}' is not available in current GDAL/OGR installation"
            )
        else:
            raise ValueError(f"Unsupported vector file extension: {ext}")

    return available_formats[ext]


def get_vector_format_from_filename(filename: str) -> str:
    """
    Determine the OGR format string from a filename.

    Only accepts extensions that are both in the supported list and have available drivers.

    Parameters:
    filename (str): Input filename with extension.

    Returns:
    str: Format string for OGR driver.

    Raises:
    ValueError: If the file extension is not supported or driver is not available.
    """
    sExtension = get_extension_from_path(filename).lower()
    if sExtension == "":
        raise ValueError(
            f"Filename '{filename}' does not have an extension to determine vector format."
        )
    return get_vector_format_from_extension(sExtension)


def get_vector_driver_from_format(file_format: str) -> ogr.Driver:
    """
    Get the OGR driver based on the provided vector file format.

    Parameters:
    file_format (str): The vector file format (e.g., 'GeoJSON', 'ESRI Shapefile', 'GPKG').

    Returns:
    ogr.Driver: The OGR driver object corresponding to the file format.

    Raises:
    ValueError: If the driver is not found.
    """
    driver = ogr.GetDriverByName(file_format)
    if driver is None:
        raise ValueError(f"Driver for format '{file_format}' not found.")
    return driver


def get_vector_driver_from_filename(filename: str) -> ogr.Driver:
    """
    Get the OGR driver based on file extension.

    Parameters:
    filename (str): Path to the input file.

    Returns:
    ogr.Driver: The OGR driver object corresponding to the file format.

    Raises:
    ValueError: If the file extension is not supported or driver is not found.
    """
    sExtension = get_extension_from_path(filename).lower()
    if sExtension == "":
        raise ValueError(
            f"Filename '{filename}' does not have an extension to determine vector format."
        )
    format_name = get_vector_format_from_extension(sExtension)
    return get_vector_driver_from_format(format_name)


@lru_cache(maxsize=1)
def check_parquet_support() -> Optional[str]:
    """Check whether an OGR driver supporting Parquet/Arrow is available.

    Returns
    -------
    Optional[str]
        The first matching OGR driver name (for example 'Parquet' or 'Arrow')
        if available, otherwise ``None``.

    Notes
    -----
    - Matches driver names that contain 'parquet' or 'arrow' (case-insensitive).
    - Result is cached for the life of the process since available drivers do
      not change at runtime.
    """
    for i in range(ogr.GetDriverCount()):
        drv = ogr.GetDriver(i)
        if drv is None:
            continue
        name = drv.GetName()
        if not name:
            continue
        lname = name.lower()
        if "parquet" in lname or "arrow" in lname:
            return name
    return None


def has_parquet_support() -> bool:
    """Return True if a Parquet/Arrow OGR driver is available."""
    return check_parquet_support() is not None
