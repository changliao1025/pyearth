import os
from osgeo import gdal

def get_gdal_driver_by_format(file_format):
    """
    Get the GDAL driver based on the provided raster file format.

    Args:
        file_format (str): The raster file format (e.g., 'GTiff', 'JP2', 'HFA').

    Returns:
        gdal.Driver: The GDAL driver object corresponding to the file format.
                      Returns None if the driver is not found.
    """
    # Get the list of GDAL drivers
    drivers = gdal.GetDriverCount()
    
    # Loop through the drivers and find the matching one
    for i in range(drivers):
        driver = gdal.GetDriver(i)
        if driver.GetMetadataItem(gdal.DCAP_RASTER) == 'YES':
            if file_format.upper() in [x.upper() for x in driver.GetMetadataItem(gdal.DMD_EXTENSION).split()]:
                return driver
    
    return None