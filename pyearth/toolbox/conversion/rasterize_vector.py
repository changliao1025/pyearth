import os
from osgeo import gdal, ogr, osr

def rasterize_vector(sFilename_vector_in, sFilename_raster_out, sAttribute_name_in,
                     dResolution_x, dResolution_y, 
                     dMissing_value_in=None,
                     iDataType_out=None,
                     dMin_x_in =None, dMax_x_in=None, dMin_y_in=None, dMax_y_in=None, 
                      pSpatialRef_in = None ):
    """
    Converts any shapefile to a raster
    :param sFilename_vector_in: STR of a shapefile name (with directory e.g., "C:/temp/poly.shp")
    :param sFilename_raster_out: STR of target file name, including directory; must end on ".tif"
    :param pixel_size: INT of pixel size (default: 10)
    :param no_data_value: Numeric (INT/FLOAT) for no-data pixels (default: -9999)
    :param iDataType_out: gdal.GDALDataType raster data type - default=gdal.GDT_Float32 (32 bit floating point)
    :kwarg field_name: name of the shapefile's field with values to burn to the raster
    :return: produces the shapefile defined with sFilename_vector_in
    """

    # open data source
    try:
        pDatasource_vector = ogr.Open(sFilename_vector_in)
    except RuntimeError as e:
        print("Error: Could not open %s." % str(sFilename_vector_in))
        return None
    
    if os.path.exists(sFilename_raster_out):
        os.remove(sFilename_raster_out)
    
    if pSpatialRef_in is None:
        pSpatialRef = pDatasource_vector.GetLayer().GetSpatialRef()
        if pSpatialRef is None:
            pSpatialRef = osr.SpatialReference()
            pSpatialRef.ImportFromEPSG(4326)
            
    else:
        pSpatialRef = pSpatialRef_in
    
    if iDataType_out is None:
        iDataType = gdal.GDT_Float32 #gdal.GDT_Int16
    else:
        iDataType = iDataType_out

    if dMissing_value_in is None:
        dMissing_value = -9999
    else:
        dMissing_value = dMissing_value_in

    pLayer_vector = pDatasource_vector.GetLayer()

    if dMin_x_in is None or dMax_x_in is None or dMin_y_in is None or dMax_y_in is None:
        # read extent
        dMin_x, dMax_x, dMin_y, dMax_y = pLayer_vector.GetExtent()
    else:
        dMin_x = dMin_x_in
        dMax_x = dMax_x_in
        dMin_y = dMin_y_in
        dMax_y = dMax_y_in

    # get x and y resolution
    nrow = int((dMax_y - dMin_y) / dResolution_y)    
    ncolumn = int((dMax_x - dMin_x) / dResolution_x)

    #get raster drive
    pRaster_Driver = gdal.GetDriverByName('GTiff')
    pDatasource_raster = pRaster_Driver.Create(sFilename_raster_out, ncolumn, nrow, 1, eType=iDataType)
    pDatasource_raster.SetGeoTransform((dMin_x, dResolution_x, 0, dMax_y, 0, -dResolution_y))
    band = pDatasource_raster.GetRasterBand(1)
    band.Fill(dMissing_value)
    band.SetNoDataValue(dMissing_value)

    # get spatial reference system and assign to raster
    pDatasource_raster.SetProjection(pSpatialRef.ExportToWkt())

    # RasterizeLayer(Dataset dataset, int bands, Layer layer, pfnTransformer=None, pTransformArg=None,
    # int burn_values=0, options=None, GDALProgressFunc callback=0, callback_data=None)
    gdal.RasterizeLayer(pDatasource_raster, [1], pLayer_vector, None, None, burn_values=[8],
                        options=["ALL_TOUCHED=TRUE", "ATTRIBUTE=" + sAttribute_name_in])

    # release raster band
    band.FlushCache()