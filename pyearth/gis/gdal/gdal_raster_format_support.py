import os
from osgeo import gdal

SUPPORTED_RASTER_FORMATS = {
    '.tif': 'GTiff',
    '.tiff': 'GTiff',
    '.jp2': 'JP2OpenJPEG',
    '.j2k': 'JP2OpenJPEG',
    '.img': 'HFA',
    '.png': 'PNG',
    '.jpg': 'JPEG',
    '.jpeg': 'JPEG',
    '.bmp': 'BMP',
    '.gif': 'GIF',
    '.nc': 'netCDF',
    '.hdf': 'HDF4',
    '.h5': 'HDF5',
    '.hdf5': 'HDF5',
    '.vrt': 'VRT',
    '.asc': 'AAIGrid',
    '.bil': 'EHdr',
    '.bsq': 'EHdr',
    '.bip': 'EHdr',
    '.rst': 'RST',
    '.grd': 'GRD'
}

def gdal_raster_format_support() -> dict[str, str]:
    """
    Return a dictionary of supported raster formats.

    Returns:
    dict[str, str]: A dictionary of supported raster formats mapping file extensions to GDAL driver names.
    """
    return SUPPORTED_RASTER_FORMATS

def print_supported_raster_formats() -> None:
    """
    Print all supported raster formats.
    """
    print("Supported raster formats:")
    for ext, driver in SUPPORTED_RASTER_FORMATS.items():
        print(f"  {ext}: {driver}")

def get_raster_format_from_extension(filename: str) -> str:
    """
    Determine the GDAL raster format string from file extension.

    Parameters:
    filename (str): Input filename with extension.

    Returns:
    str: Format string for GDAL driver.

    Raises:
    ValueError: If the file extension is not supported.
    """
    _, ext = os.path.splitext(filename.lower())

    if ext in SUPPORTED_RASTER_FORMATS:
        return SUPPORTED_RASTER_FORMATS[ext]
    else:
        raise ValueError(f"Unsupported raster file extension: {ext}")

def get_raster_driver(file_format: str) -> gdal.Driver:
    """
    Get the GDAL driver based on the provided raster file format.

    Parameters:
    file_format (str): The raster file format (e.g., 'GTiff', 'JP2OpenJPEG', 'HFA').

    Returns:
    gdal.Driver: The GDAL driver object corresponding to the file format.

    Raises:
    ValueError: If the driver is not found.
    """
    driver = gdal.GetDriverByName(file_format)
    if driver is None:
        raise ValueError(f"Driver for format '{file_format}' not found.")
    return driver

def get_raster_driver_from_extension(filename: str) -> gdal.Driver:
    """
    Get the GDAL driver based on file extension.

    Parameters:
    filename (str): Input filename with extension.

    Returns:
    gdal.Driver: The GDAL driver object corresponding to the file format.

    Raises:
    ValueError: If the file extension is not supported or driver is not found.
    """
    format_name = get_raster_format_from_extension(filename)
    return get_raster_driver(format_name)