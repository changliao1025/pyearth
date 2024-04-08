import os
import numpy as np
from osgeo import gdal, ogr, osr

def rasterize_vector(sFilename_vector_in, sFilename_raster_out, 
                     dResolution_x, dResolution_y, 
                     iFlag_boundary_only_in = None,
                     iFlag_use_field_value_in=None,
                     sAttribute_name_in=None,
                     dMissing_value_in=None,
                     dField_value_in=None,
                     dFill_value_in=None,
                     iDataType_out=None,
                     dMin_x_in =None, dMax_x_in=None, dMin_y_in=None, dMax_y_in=None, 
                     nRow_in=None, nColumn_in=None,
                     pProjection_target_in = None ):
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
    
    # Create a new vector layer for the boundary
    pDatasource_boundary = ogr.GetDriverByName('Memory').CreateDataSource('out')
    
    if os.path.exists(sFilename_raster_out):
        os.remove(sFilename_raster_out)
    
    if sAttribute_name_in is not None:        
        sAttribute_name = sAttribute_name_in

    if iFlag_use_field_value_in is None:
        iFlag_use_field_value = 0
    else:
        iFlag_use_field_value = 1
    
    if dField_value_in is None:
        dField_value = 1
    else:
        dField_value = dField_value_in

    if iFlag_boundary_only_in is None:
        iFlag_boundary_only = 1
    else:
        iFlag_boundary_only = iFlag_boundary_only_in
    
    if dFill_value_in is None:
        dFill_value = 1
    else:
        dFill_value = dFill_value_in
    
    pLayer_vector = pDatasource_vector.GetLayer()
    pSpatialRef_source = pLayer_vector.GetSpatialRef()
    pProjection_source = pSpatialRef_source.ExportToWkt()
    if pProjection_target_in is None:
        pProjection_target = pProjection_source
    else:
        pProjection_target = pProjection_target_in  
    
    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(pProjection_target)
    
    if iDataType_out is None:
        iDataType = gdal.GDT_Float32 #gdal.GDT_Int16
    else:
        iDataType = iDataType_out

    if dMissing_value_in is None:
        dMissing_value = -9999
    else:
        dMissing_value = dMissing_value_in    

    if dMin_x_in is None or dMax_x_in is None or dMin_y_in is None or dMax_y_in is None:
        # read extent
        dMin_x, dMax_x, dMin_y, dMax_y = pLayer_vector.GetExtent()
    else:
        dMin_x = dMin_x_in
        dMax_x = dMax_x_in
        dMin_y = dMin_y_in
        dMax_y = dMax_y_in       

    #check the geometry type is polygon or not
    pFeature = pLayer_vector.GetNextFeature()
    pGeometry = pFeature.GetGeometryRef()
    if pGeometry.GetGeometryName() == 'POLYGON':
        pass
    else:
        if pGeometry.GetGeometryName() == 'MULTIPOLYGON':
            print('This is multipolygon')
        else:
            iFlag_boundary_only = 1    
     
    if nRow_in is not None:
        nrow = nRow_in
    else:
        nrow = int((dMax_y - dMin_y) / dResolution_y)    
    if nColumn_in is not None:
        ncolumn = nColumn_in
    else:
        ncolumn = int((dMax_x - dMin_x) / dResolution_x)

    if nrow == 0 or ncolumn == 0:
        #print('Error: resolution is too high')
        return

    #get raster drive
    pRaster_Driver = gdal.GetDriverByName('GTiff')
    pDatasource_raster = pRaster_Driver.Create(sFilename_raster_out, ncolumn, nrow, 1, eType=iDataType)
    pDatasource_raster.SetGeoTransform((dMin_x, dResolution_x, 0, dMax_y, 0, -dResolution_y))
    band = pDatasource_raster.GetRasterBand(1)    
    band.Fill(dMissing_value)
    band.SetNoDataValue(dMissing_value)

    # get spatial reference system and assign to raster
    pDatasource_raster.SetProjection(pSpatialRef_target.ExportToWkt())

    if iFlag_use_field_value == 1:
        if iFlag_boundary_only == 0:
            pLayer_boundary = pDatasource_boundary.CreateLayer('boundary', srs=pSpatialRef_source)
            geometry = pFeature.GetGeometryRef()        
            pFeature_boundary = ogr.Feature(pLayer_boundary.GetLayerDefn()) 
            pFeature_boundary.SetGeometry(geometry.Boundary())
            pLayer_boundary.CreateFeature(pFeature_boundary)
            gdal.RasterizeLayer(pDatasource_raster, [1], pLayer_vector, None, None, burn_values=[dField_value],
                    options=["ALL_TOUCHED=TRUE"])                
            gdal.RasterizeLayer(pDatasource_raster, [1], pLayer_boundary, None, None,
                options=["ALL_TOUCHED=TRUE"])
        else:
            gdal.RasterizeLayer(pDatasource_raster, [1], pLayer_vector, None, None, 
                    options=["ALL_TOUCHED=TRUE", "ATTRIBUTE=" + sAttribute_name])
    else: #use user provided value for field, 
        if iFlag_boundary_only == 0: #fill the polygon with a different value
            pLayer_boundary = pDatasource_boundary.CreateLayer('boundary', srs=pSpatialRef_source)
            geometry = pFeature.GetGeometryRef()        
            pFeature_boundary = ogr.Feature(pLayer_boundary.GetLayerDefn()) 
            pFeature_boundary.SetGeometry(geometry.Boundary())
            pLayer_boundary.CreateFeature(pFeature_boundary)
            gdal.RasterizeLayer(pDatasource_raster, [1], pLayer_vector, None, None, burn_values=[dFill_value],
                    options=["ALL_TOUCHED=TRUE"])                
            gdal.RasterizeLayer(pDatasource_raster, [1], pLayer_boundary, None, None, burn_values=[dField_value],
                options=["ALL_TOUCHED=TRUE"])
        else:
            gdal.RasterizeLayer(pDatasource_raster, [1], pLayer_vector, None, None, burn_values=[dField_value],
                    options=["ALL_TOUCHED=TRUE"])

    # release raster band
    band.FlushCache()
    pDatasource_boundary.Destroy()
    pDatasource_raster = None
    pSpatialRef_source = None
    pSpatialRef_target = None

    return