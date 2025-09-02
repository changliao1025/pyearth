import os
import logging
import numpy as np
from osgeo import gdal, ogr, osr
from pyearth.gis.gdal.gdal_to_numpy_datatype import numpy_to_gdal_type


def vectorize_raster(sFilename_raster_in, sFilename_vector_out, sFieldname='value',
                      sFieldtype=ogr.OFTInteger64):
    """
    Converts a binary raster file to a vector shapefile.

    Parameters:
    sFilename_raster_in (str): Path to the input raster file.
    sFilename_vector_out (str): Path to the output vector file.
    sFieldname (str): Name of the field to create in the vector file.
    sFieldtype (ogr.FieldType): Type of the field to create in the vector file.

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

    pGeotransform = pDataset_raster.GetGeoTransform()
    dOriginX = pGeotransform[0]
    dOriginY = pGeotransform[3]
    dPixelWidth = pGeotransform[1]
    dPixelHeight = pGeotransform[5]

    pBand = pDataset_raster.GetRasterBand(1)  # Assuming the binary mask is in the first band

    # Get the projection of the raster
    pProjection_source = pDataset_raster.GetProjection()
    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromWkt(pProjection_source)
    nodata_value = pBand.GetNoDataValue()
    eType_out = pBand.DataType

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

    # Create a field for id first
    # This field will be used to store the custom IDs for each polygon
    pField_defn = ogr.FieldDefn('id', ogr.OFTInteger64)
    pLayer_out.CreateField(pField_defn)
    # Create a field for the value

    #create another field for the value

    pField_defn = ogr.FieldDefn(sFieldname, sFieldtype)
    pLayer_out.CreateField(pField_defn)

    # Step 3: Polygonize
    #if gdal.Polygonize(pBand, mask_band, pLayer_out, -1, [], callback=None) != 0:
    #    return False

    # Step 3: Polygonize with custom pFeature IDs
    # Create a list of options for the polygonizer
    options = []

    # Set field index to 0 (the first field you created)
    field_index = 0

    # Use a callback function to assign custom IDs
    def polygonize_callback(progress, message, callback_data):
        # This function is called during the polygonization process
        # You can use it to monitor progress and assign custom IDs
        return 1  # Return 1 to continue, 0 to abort

    if gdal.Polygonize(pBand, mask_band, pLayer_out, field_index, options, callback=polygonize_callback) != 0:
        return False


    #reset the reading of the layer
    pLayer_out.ResetReading()
    pFeature = pLayer_out.GetNextFeature()
    lID = 1
    while pFeature:
        # Get the geometry

        # Get the original pixel value at this location if needed
        geom = pFeature.GetGeometryRef()
        pPolygonWKT = geom.ExportToWkt()
        aPolygon_extent= geom.GetEnvelope()
        minX, maxX, minY, maxY = aPolygon_extent
        aPolygon_extent = [minX, minY, maxX, maxY]
        #use the boundary to clip the raster
        pWrapOption = gdal.WarpOptions( cropToCutline=False,
                                           #cutlineDSName = sFilename_clip ,
                                           cutlineWKT=pPolygonWKT,
                                        xRes=dPixelWidth,
                                        yRes=abs(dPixelHeight),
                                        outputBounds=aPolygon_extent,
                                        dstSRS=spatialRef, format = 'MEM',
                                        resampleAlg='near',
                                         dstNodata=nodata_value,
                                         outputType=eType_out)
        pDataset_clip_warped = gdal.Warp('', sFilename_raster_in, options=pWrapOption)
        newGeoTransform = pDataset_clip_warped.GetGeoTransform()
        #convert the warped dataset to an array
        aData_clip = pDataset_clip_warped.ReadAsArray()
        #get the dominant value (most frequent) in the array
        if aData_clip.size > 0:
            # Handle potential NaN or nodata values
            valid_data = aData_clip.flatten()
            # Remove nodata values if they exist
            if nodata_value is not None:
                valid_data = valid_data[valid_data != nodata_value]

            # Find most frequent value if there's data left
            if valid_data.size > 0:
                if np.issubdtype(valid_data.dtype, np.integer):
                    # For integer arrays
                    dominant_value = np.bincount(valid_data).argmax()
                else:
                    # For float arrays
                    unique_values, counts = np.unique(valid_data, return_counts=True)
                    dominant_value = unique_values[np.argmax(counts)]
            else:
                dominant_value = nodata_value
        else:
            dominant_value = nodata_value

        # Set the new ID
        pFeature.SetField('id', lID)
        # Set additional attributes if needed
        dominant_value_gdal = numpy_to_gdal_type(dominant_value, sFieldtype )

        pFeature.SetField(sFieldname, dominant_value_gdal)

        # For example, set the original pixel value
        # Update the pFeature
        pLayer_out.SetFeature(pFeature)

        # Move to the next feature
        lID += 1
        pFeature = pLayer_out.GetNextFeature()

    # Cleanup
    pDataset_shapefile = None

    print('Vectorization completed successfully.')
    return True