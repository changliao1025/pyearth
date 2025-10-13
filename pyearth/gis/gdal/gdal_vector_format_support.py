import os
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

    if ext in SUPPORTED_VECTOR_FORMATS:
        return SUPPORTED_VECTOR_FORMATS[ext]
    else:
        raise ValueError(f"Unsupported vector file extension: {ext}")

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
    format_name = get_vector_format_from_extension(filename)
    return get_vector_driver(format_name)