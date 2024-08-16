import os
import logging
import numpy as np
from osgeo import gdal, ogr, osr


def vectorize_raster(sFilename_raster_in, sFilename_vector_out, field_name='ID', field_type=ogr.OFTInteger64):
    """
    Converts a binary raster file to a vector shapefile.

    Parameters:
    sFilename_raster_in (str): Path to the input raster file.
    sFilename_vector_out (str): Path to the output vector file.
    field_name (str): Name of the field to create in the vector file.
    field_type (ogr.FieldType): Type of the field to create in the vector file.

    Returns:
    bool: True if successful, False otherwise.
    """
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    # Step 1: Open the binary raster file
    pDataset_raster = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)
    if pDataset_raster is None:
        logging.error('Could not open %s.', sFilename_raster_in)
        return False

    if os.path.exists(sFilename_vector_out):
        pDriver_shapefile.DeleteDataSource(sFilename_vector_out)

    pBand = pDataset_raster.GetRasterBand(1)  # Assuming the binary mask is in the first band

    # Get the projection of the raster
    pProjection_source = pDataset_raster.GetProjection()
    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromWkt(pProjection_source)
    nodata_value = pBand.GetNoDataValue()

    # Create a mask band
    mask_ds = gdal.GetDriverByName('MEM').Create('', pDataset_raster.RasterXSize, pDataset_raster.RasterYSize, 1, gdal.GDT_Byte)
    mask_band = mask_ds.GetRasterBand(1)

    # Populate mask band
    raster_array = pBand.ReadAsArray()
    mask_array = np.not_equal(raster_array, nodata_value).astype(np.uint8)
    mask_band.WriteArray(mask_array)
    mask_band.FlushCache()

    # Step 2: Create an output Shapefile to store the vector data

    pDataset_shapefile = pDriver_shapefile.CreateDataSource(sFilename_vector_out)
    pLayer_out = pDataset_shapefile.CreateLayer('polygonized', spatialRef, ogr.wkbPolygon)

    # Create a field
    pField_defn = ogr.FieldDefn(field_name, field_type)
    pLayer_out.CreateField(pField_defn)

    # Step 3: Polygonize
    if gdal.Polygonize(pBand, mask_band, pLayer_out, -1, [], callback=None) != 0:
        return False

    # Cleanup
    pDataset_shapefile = None

    print('Vectorization completed successfully.')
    return True