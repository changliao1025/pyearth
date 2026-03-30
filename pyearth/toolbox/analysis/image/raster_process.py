
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

def reset_antimeridian_column(column_data, source_value, target_value, iRaster_buffer_pixel):
    """
    Reset values in a column corresponding to the antimeridian to prevent buffer dilation artifacts.

    This function identifies contiguous segments of the target value in the column and preserves
    only a specified number of pixels at the top and bottom of each segment, while setting
    the rest to NaN. This helps to prevent buffer dilation from propagating across the dateline.

    Parameters:
    -----------
    column_data : 2D numpy array, shape (nrow, n_buffer_cols)
        The input column data from the raster. Shape should be (nrow, iRaster_buffer_pixel+1).
        Each column is processed independently.
    source_value : int or float
        The value to search for and partially replace (e.g., 1 for coast).
    target_value : int or float
        The value to assign to middle pixels of segments (e.g., 2 for land).
        #example,source_value is 1 (coast), target_value is 2(land), nan is ocean,
        # then we want to keep the top and bottom pixels of each segment of 1s, and set the middle as 2
        # Nan Nan Nan
        # Nan Nan Nan
        # Nan Nan Nan
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # Nan Nan Nan
        # Nan Nan Nan
        # Nan Nan Nan
        # we want to keep the (2*iRaster_buffer_pixel + 1) top 1s, and (2*iRaster_buffer_pixel + 1) bot 1s, for each segment shown in the column
        # then we can reset the data as
        # Nan Nan Nan
        # Nan Nan Nan
        # Nan Nan Nan
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 2 2 2
        # 1 1 1
        # 1 1 1
        # 1 1 1
        # Nan Nan Nan
        # Nan Nan Nan
        # Nan Nan Nan
        #remember that there could be multiple segments in the column,
        # for example, there could be multiple islands in the column, so we need to process each segment separately

    Returns:
    --------
    updated_columns : 2D numpy array, shape (nrow, n_buffer_cols)
        The modified column data with each column processed independently.
    """
    # Check if input is 2D
    if column_data.ndim != 2:
        raise ValueError(f"Expected 2D array, got {column_data.ndim}D array with shape {column_data.shape}")

    nrow, ncol = column_data.shape

    # Initialize output array
    updated_columns = np.copy(column_data)

    # Process each column independently
    for col_idx in range(ncol):
        current_column = column_data[:, col_idx]

        # Identify indices where the source_value is present
        source_indices = np.where(current_column == source_value)[0]

        if source_indices.size > 0:
            # Find contiguous segments of source_value
            segments = np.split(source_indices, np.where(np.diff(source_indices) != 1)[0] + 1)

            for segment in segments:
                if segment.size > 0:
                    # Preserve iRaster_buffer_pixel pixels at top and bottom, set middle to target_value
                    if len(segment) > 2 * iRaster_buffer_pixel:
                        # Only modify middle pixels if segment is large enough
                        middle_pixels = segment[iRaster_buffer_pixel+1:-iRaster_buffer_pixel-1] #when iRaster_buffer_pixel is 1, we want to keep the first 2 pixels and the last 2 pixels, and set the middle pixels as target_value
                        updated_columns[middle_pixels, col_idx] = target_value
                    # If segment is too small, keep all as source_value (no modification)

    return updated_columns

def fix_raster_antimeridian_issue(sFilename_in,
                                  sFilename_out,
                                  source_value,
                                  target_value,
                                  iRaster_buffer_pixel=1):
    """
    Replace raster values at the international dateline and surrounding buffer columns
    to avoid artifacts when creating buffer zones.

    This function identifies columns at the antimeridian (±180° longitude) and replaces
    values in those columns plus a configurable buffer on each side. This prevents
    buffer dilation operations from incorrectly propagating across the dateline.

    For global rasters (-180° to 180°), the antimeridian wraps around at the edges:
    - Left side (near -180°): last iRaster_buffer_pixel columns (right edge)
    - Right side (near +180°): first iRaster_buffer_pixel columns (left edge)

    Parameters:
    -----------
    sFilename_in : str
        Input raster file path
    sFilename_out : str
        Output raster file path
    target_value : int or float
        Value to assign to antimeridian buffer columns
    iRaster_buffer_pixel : int, optional (default=1)
        Number of buffer columns to include on each side of the dateline

    Returns:
    --------
    None

    Example:
    --------
    With iRaster_buffer_pixel=2 and dateline at edges (global raster with 360 columns):
    - Columns selected: [358, 359, 0, 1] (4 columns total)
    - All selected columns set to target_value
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

    # Identify columns corresponding to the antimeridian using normalized longitude
    ncols = raster_data.shape[1]
    lon = geo_transform[0] + np.arange(ncols) * geo_transform[1]
    lon180 = ((lon + 180.0) % 360.0) - 180.0

    # Check if this raster crosses the antimeridian
    tolerance = abs(geo_transform[1]) * 0.5
    dateline_columns = np.where(np.isclose(np.abs(lon180), 180.0, atol=tolerance))[0]

    # Process only if the antimeridian exists in this raster extent
    if dateline_columns.size > 0:
        # For global rasters (-180 to 180), the antimeridian wraps around:
        # - Left side (near -180°): last iRaster_buffer_pixel columns (right edge)
        # - Right side (near +180°): first iRaster_buffer_pixel columns (left edge)

        # Select buffer columns from the left edge (wrapping to left/west at -180°)
        left_start = 0
        left_buffer_cols = np.arange(left_start, iRaster_buffer_pixel+ 1)
        left_columns = raster_data[:, left_buffer_cols]
        right_end = ncols
        right_buffer_cols = np.arange(ncols - iRaster_buffer_pixel-1, right_end)
        right_columns = raster_data[:, right_buffer_cols]

        #now process the left, for example, when iRaster_buffer_pixel = 1, the data can be

        updated_columns_left = reset_antimeridian_column(left_columns, source_value, target_value, iRaster_buffer_pixel)
        # Select buffer columns from the right edge

        updated_columns_right = reset_antimeridian_column(right_columns, source_value, target_value, iRaster_buffer_pixel)
        #now replace the original columns with the updated columns
        raster_data[:, left_buffer_cols] = updated_columns_left
        raster_data[:, right_buffer_cols] = updated_columns_right
    else:
        print("No antimeridian detected in this raster extent")


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
