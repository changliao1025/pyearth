
import os
import numpy as np
from osgeo import gdal
from scipy.ndimage import binary_dilation
from pyearth.gis.gdal.gdal_raster_format_support import get_raster_driver_from_filename
from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import gdal_write_geotiff_file

def create_raster_buffer_zone(sFilename_in,
                                  sFilename_out, target_value,
                                  buffer_pixel):
    """
    Create a buffer zone for a raster image based on a target value and buffer pixel size.

    Parameters:
    sFilename_in (str): The input raster file path.
    sFilename_out (str): The output raster file path.
    target_value (int or float): The pixel value to create the buffer zone around.
    buffer_pixel (int): The number of pixels to include in the buffer zone.

    Returns:
    None
    """
    # Open the input raster file
    dataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    if dataset is None:
        print("Could not open file: ", sFilename_in)
        return
    # Read the raster band
    band = dataset.GetRasterBand(1)
    if band is None:
        print("Could not read band 1 from: ", sFilename_in)
        dataset = None
        return

    geo_transform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()
    nodata_value = band.GetNoDataValue()
    output_datatype = band.DataType

    # Read the raster data as a numpy array
    raster_data = band.ReadAsArray()
    # Create a reusable buffer and build a binary mask where the target value is located
    buffer_array = np.zeros_like(raster_data, dtype=bool)
    np.equal(raster_data, target_value, out=buffer_array)
    # Perform binary dilation in-place on the buffer
    if buffer_pixel > 0:
        binary_dilation(buffer_array, iterations=buffer_pixel, output=buffer_array)
    # Overwrite only buffered cells with the target value
    raster_data[buffer_array] = target_value

    #there is an issue at the international dateline, we done want to make the IDL as coastline,


    # Create a new raster file to save the buffer zone
    if os.path.exists(sFilename_out):
        os.remove(sFilename_out)
    gdal_write_geotiff_file(
        sFilename_out,
        raster_data,
        geo_transform[1],
        geo_transform[5],
        geo_transform[0],
        geo_transform[3],
        nodata_value,
        projection,
        datatype=output_datatype,
    )

    print("Buffer zone created and saved to: ", sFilename_out)

def fix_raster_antimeridian_issue(sFilename_in,
                                  sFilename_out,
                                  target_value):

    #this function replace the raster value at the international dateline with the target value, to avoid the issue of creating a buffer zone at the IDL

    # Open the input raster file
    dataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    if dataset is None:
        print("Could not open file: ", sFilename_in)
        return
    # Read the raster band
    band = dataset.GetRasterBand(1)
    if band is None:
        print("Could not read band 1 from: ", sFilename_in)
        dataset = None
        return

    geo_transform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()
    nodata_value = band.GetNoDataValue()
    output_datatype = band.DataType
    # Read the raster data as a numpy array
    raster_data = band.ReadAsArray()
    # Identify columns corresponding to the antimeridian using normalized longitude
    ncols = raster_data.shape[1]
    lon = geo_transform[0] + np.arange(ncols) * geo_transform[1]
    lon180 = ((lon + 180.0) % 360.0) - 180.0
    tolerance = abs(geo_transform[1]) * 0.5
    dateline_columns = np.where(np.isclose(np.abs(lon180), 180.0, atol=tolerance))[0]

    # Replace values only if the antimeridian exists in this raster extent
    if dateline_columns.size > 0:
        raster_data[:, dateline_columns] = target_value
    # Create a new raster file to save the modified data
    if os.path.exists(sFilename_out):
        os.remove(sFilename_out)
    gdal_write_geotiff_file(
        sFilename_out,
        raster_data,
        geo_transform[1],
        geo_transform[5],
        geo_transform[0],
        geo_transform[3],
        nodata_value,
        projection,
        datatype=output_datatype,
    )
    print("International dateline fixed and saved to: ", sFilename_out)
