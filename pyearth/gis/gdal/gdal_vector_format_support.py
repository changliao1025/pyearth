import os
from functools import lru_cache
from typing import Optional
from osgeo import ogr

SUPPORTED_VECTOR_FORMATS = {
    '.geojson': 'GeoJSON',
    '.json': 'GeoJSON',
    '.shp': 'ESRI Shapefile',
    '.gpkg': 'GPKG',
    '.kml': 'KML',
    '.gml': 'GML',
    '.sqlite': 'SQLite',
    '.parquet': 'Parquet',
    '.geoparquet': 'Parquet',
    '.csv': 'CSV',
    '.tab': 'MapInfo File',
    '.mif': 'MapInfo File',
    '.dxf': 'DXF',
    '.gdb': 'OpenFileGDB'
}


__all__ = [
    'SUPPORTED_VECTOR_FORMATS',
    'gdal_vector_format_support',
    'print_supported_vector_formats',
    'get_vector_format_from_extension',
    'get_vector_driver',
    'get_vector_driver_from_extension',
    'check_parquet_support',
    'has_parquet_support',
    'get_preferred_vector_driver_for_extension',
]

def gdal_vector_format_support() -> dict[str, str]:
    """
    Return a dictionary of supported vector formats.

    Returns:
    dict[str, str]: A dictionary of supported vector formats mapping file extensions to OGR driver names.
    """
    return SUPPORTED_VECTOR_FORMATS

def print_supported_vector_formats() -> None:
    """
    Print all supported vector formats.
    """
    print("Supported vector formats:")
    for ext, driver in SUPPORTED_VECTOR_FORMATS.items():
        print(f"  {ext}: {driver}")

def get_vector_format_from_extension(filename: str) -> str:
    """
    Determine the OGR format string from file extension.

    Parameters:
    filename (str): Input filename with extension.

    Returns:
    str: Format string for OGR driver.

    Raises:
    ValueError: If the file extension is not supported.
    """
    _, ext = os.path.splitext(filename.lower())

    if ext not in SUPPORTED_VECTOR_FORMATS:
        raise ValueError(f"Unsupported vector file extension: {ext}")

    # For Parquet/GeoParquet we need to ensure a runtime driver is available.
    if ext in ('.parquet', '.geoparquet'):
        drv_name = check_parquet_support()
        if drv_name:
            return drv_name
        # No parquet driver available; raise a clear error so callers can
        # consider a fallback or instruct the user to install GDAL with
        # Parquet/Arrow support.
        raise ValueError(
            "Parquet/GeoParquet driver not available in OGR. "
            "Install GDAL with Parquet/Arrow support or use an alternative format (e.g. GPKG)."
        )

    return SUPPORTED_VECTOR_FORMATS[ext]

def get_vector_driver(file_format: str) -> ogr.Driver:
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

def get_vector_driver_from_extension(filename: str) -> ogr.Driver:
    """
    Get the OGR driver based on file extension.

    Parameters:
    filename (str): Input filename with extension.

    Returns:
    ogr.Driver: The OGR driver object corresponding to the file format.

    Raises:
    ValueError: If the file extension is not supported or driver is not found.
    """
    # Prefer a runtime-resolvable driver for extensions like Parquet.
    return get_preferred_vector_driver_for_extension(filename)


def get_preferred_vector_driver_for_extension(filename: str, fallback: Optional[str] = None) -> ogr.Driver:
    """Return an OGR driver suitable for the given filename extension.

    Parameters
    ----------
    filename : str
        File name (used to determine extension).
    fallback : str, optional
        Optional fallback driver name to use when the preferred driver is not
        available (for example, 'GPKG'). If provided and available, it will be
        returned instead of raising.

    Returns
    -------
    ogr.Driver
        The OGR driver instance corresponding to the filename extension.

    Raises
    ------
    ValueError
        If no suitable driver is available for the extension.
    """
    _, ext = os.path.splitext(filename.lower())
    if ext not in SUPPORTED_VECTOR_FORMATS:
        raise ValueError(f"Unsupported vector file extension: {ext}")

    # Special handling for Parquet/GeoParquet: resolve runtime driver name
    if ext in ('.parquet', '.geoparquet'):
        drv_name = check_parquet_support()
        if drv_name:
            drv = ogr.GetDriverByName(drv_name)
            if drv:
                return drv
        if fallback:
            drv = ogr.GetDriverByName(fallback)
            if drv:
                return drv
        raise ValueError(
            "Parquet/GeoParquet driver not available. Install GDAL with Parquet/Arrow support or provide a fallback driver."
        )

    # Normal path: use the mapping and get the driver
    format_name = SUPPORTED_VECTOR_FORMATS[ext]
    drv = ogr.GetDriverByName(format_name)
    if drv is None:
        # Allow fallback if provided
        if fallback:
            drv = ogr.GetDriverByName(fallback)
            if drv:
                return drv
        raise ValueError(f"Driver for format '{format_name}' not found.")
    return drv

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
        if 'parquet' in lname or 'arrow' in lname:
            return name
    return None


def has_parquet_support() -> bool:
    """Return True if a Parquet/Arrow OGR driver is available."""
    return check_parquet_support() is not None