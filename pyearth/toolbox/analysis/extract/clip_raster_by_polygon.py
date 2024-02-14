import os, sys
import numpy as np
from osgeo import gdal, osr, ogr, gdalconst
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file


def clip_raster_by_geojson(sFilename_raster_in, sFilename_geojson_in, sFilename_raster_out, sFormat='GTiff'):
    """
    Clip a raster by a polygon
    :param sFilename_raster_in: input raster filename
    :param sFilename_polygon_in: input polygon filename
    :param sFilename_raster_out: output raster filename
    :param sFormat: output format
    :return: None
    """

    #check input files
    if os.path.exists(sFilename_raster_in):
        pass
    else:
        print('The raster file does not exist!')
        return
    
    if os.path.exists(sFilename_geojson_in):
        pass
    else:   
        print('The geojson does not exist!')
        return
    
    sDriverName = 'GTiff'
    pDriver = gdal.GetDriverByName(sDriverName)
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    pDataset_elevation = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)
   
    dummy = gdal_read_geotiff_file(sFilename_raster_in)
    aDem_in= dummy['dataOut']    
    eType = dummy['dataType']  
    dPixelWidth = dummy['pixelWidth']                        
    pPixelHeight = dummy['pixelHeight']
    dOriginX = dummy['originX']
    dOriginY = dummy['originY']
    nrow = dummy['nrow']
    ncolumn = dummy['ncolumn']
    dMissing_value= dummy['missingValue']
    pProjection = dummy['projection']
    pSpatialRef_target = dummy['spatialReference']
    wkt1 = pSpatialRef_target.ExportToWkt()
    print(wkt1)
    dX_left=dOriginX
    dX_right = dOriginX + ncolumn * dPixelWidth
    dY_top = dOriginY
    dY_bot = dOriginY + nrow * pPixelHeight

    #get the spatial reference of the shapefile
    shapefile_ds = ogr.Open(sFilename_geojson_in)    
    # Get the first layer in the shapefile
    pLayer_clip = shapefile_ds.GetLayer(0)   
    # Count the number of features (polygons)
    num_polygons = pLayer_clip.GetFeatureCount()
    # Get the spatial reference of the layer
    pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
    wkt2 = pSpatial_reference_clip.ExportToWkt()
    print(wkt2)
    #comparison = pSpatialRef_target.IsSame(pSpatial_reference_clip)

    if( wkt1 != wkt2):
    #if(comparison != 1):
        iFlag_transform = 1
        transform = osr.CoordinateTransformation(pSpatial_reference_clip, pSpatialRef_target)
        #in this case, we can reproject the shapefile to the same spatial reference as the raster
        #get the folder that contains the shapefile
        sFolder = os.path.dirname(sFilename_geojson_in)
        #get the name of the shapefile
        sName = os.path.basename(sFilename_geojson_in)
        #get the name of the shapefile without extension
        sName_no_extension = os.path.splitext(sName)[0]
        #create a new shapefile
        sFilename_shapefile_out = sFolder + '/' + sName_no_extension + '_transformed.shp'
        #create a new shapefile
        
        #create a new shapefile
        pDataset_transform = pDriver_shapefile.CreateDataSource(sFilename_shapefile_out)
        #create a new shapefile layer
        pLayer_transform = pDataset_transform.CreateLayer('boundary', pSpatialRef_target, geom_type=ogr.wkbPolygon)
        #create a new feature
        pFeature_transform = ogr.Feature(pLayer_transform.GetLayerDefn())
        if num_polygons ==1:
            pFeature_clip = pLayer_clip.GetNextFeature()
            pPolygon = pFeature_clip.GetGeometryRef()
            pPolygon.Transform(transform)            
            pass
        else:
            print('The shapefile contains more than one polygon!')
            #merge the polygons
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            # Iterate over all features in the layer
            for i in range(num_polygons):
                pFeature_clip = pLayer_clip.GetFeature(i)
                pGeometry_clip = pFeature_clip.GetGeometryRef()                
                if pGeometry_clip is not None:
                    pGeometry_clip.Transform(transform)     
                    # Union the geometry of each feature with the merged polygon
                    pPolygon = pPolygon.Union(pGeometry_clip)
             
        
        #set the geometry
        pFeature_transform.SetGeometry(pPolygon)
        #add the feature to the layer
        pLayer_transform.CreateFeature(pFeature_transform)        
        pDataset_transform.FlushCache()
        #use the new shapefile to clip the raster
        sFilename_clip = sFilename_shapefile_out

    else:
        sFilename_clip = sFilename_geojson_in        
        if num_polygons ==1:
           #read the first polygon 
            pFeature_clip = pLayer_clip.GetNextFeature()
            pPolygon = pFeature_clip.GetGeometryRef() 
        else:
            print('The shapefile contains more than one polygon!')
            #merge the polygons
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            # Iterate over all features in the layer
            for i in range(num_polygons):
                pFeature_clip = pLayer_clip.GetFeature(i)
                geometry = pFeature_clip.GetGeometryRef()                
                if geometry is not None:                    
                    # Union the geometry of each feature with the merged polygon
                    pPolygon = pPolygon.Union(geometry)
             

        #get the envelope of the polygon

        iFlag_transform = 0

    #use the gdal warp function to clip the raster
    minX, maxX, minY, maxY = pPolygon.GetEnvelope()
    iNewWidth = int( (maxX - minX) / abs(dPixelWidth)  )
    iNewHeigh = int( (maxY - minY) / abs(dPixelWidth) )
    newGeoTransform = (minX, dPixelWidth, 0,  maxY, 0, -dPixelWidth)  
    if minX > dX_right or maxX < dX_left or minY > dY_top or maxY < dY_bot:        
        #this polygon is out of bound            
        pass
    else:       
        pDataset_clip = pDriver.Create(sFilename_raster_out, iNewWidth, iNewHeigh, 1, eType)
        pDataset_clip.SetGeoTransform( newGeoTransform )
        pDataset_clip.SetProjection( pProjection)   
        pWrapOption = gdal.WarpOptions( cropToCutline=True,cutlineDSName = sFilename_clip , 
                width=iNewWidth,   
                    height=iNewHeigh,      
                        dstSRS=pProjection , format = 'MEM' )
        pDataset_clip_warped = gdal.Warp('', pDataset_elevation, options=pWrapOption)

        #convert the warped dataset to an array
        aDem_clip = pDataset_clip_warped.ReadAsArray()
        #change the missing value to the original missing value
        aDem_clip[aDem_clip == dMissing_value] = -9999
        # Write the warped dataset to the output raster
        pDataset_clip.GetRasterBand(1).WriteArray(aDem_clip)

        
        #close the dataset
        pDataset_clip = None
        #close the dataset
        pDataset_elevation = None
        print('The raster has been clipped by the geojson!')
    


def clip_raster_by_shapefile(sFilename_raster_in, sFilename_shapefile_in, sFilename_raster_out, sFormat='GTiff'):
    """
    Clip a raster by a polygon
    :param sFilename_raster_in: input raster filename
    :param sFilename_polygon_in: input polygon filename
    :param sFilename_raster_out: output raster filename
    :param sFormat: output format
    :return: None
    """

    #check input files
    if os.path.exists(sFilename_raster_in):
        pass
    else:
        print('The raster file does not exist!')
        return
    
    if os.path.exists(sFilename_shapefile_in):
        pass
    else:   
        print('The shapefile does not exist!')
        return
    
    sDriverName = 'GTiff'
    pDriver = gdal.GetDriverByName(sDriverName)
    pDataset_elevation = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)
   
    dummy = gdal_read_geotiff_file(sFilename_raster_in)
    aDem_in= dummy['dataOut']
    dPixelWidth = dummy['pixelWidth']                        
    pPixelHeight = dummy['pixelHeight']
    dOriginX = dummy['originX']
    dOriginY = dummy['originY']
    nrow = dummy['nrow']
    ncolumn = dummy['ncolumn']
    dMissing_value= dummy['missingValue']
    pProjection = dummy['projection']
    pSpatialRef_target = dummy['spatialReference']
    dX_left=dOriginX
    dX_right = dOriginX + ncolumn * dPixelWidth
    dY_top = dOriginY
    dY_bot = dOriginY + nrow * pPixelHeight

    #get the spatial reference of the shapefile
    shapefile_ds = ogr.Open(sFilename_shapefile_in)    
    # Get the first layer in the shapefile
    layer = shapefile_ds.GetLayer(0)   

    # Count the number of features (polygons)
    num_polygons = layer.GetFeatureCount()
    # Get the spatial reference of the layer
    pSpatial_reference_b = layer.GetSpatialRef()

    comparison = pSpatialRef_target.IsSame(pSpatial_reference_b)
    
    if(comparison != 1):
        iFlag_transform = 1
        transform = osr.CoordinateTransformation(pSpatialRef_target, pSpatial_reference_b)
        #in this case, we can reproject the shapefile to the same spatial reference as the raster
        #get the folder that contains the shapefile
        sFolder = os.path.dirname(sFilename_shapefile_in)
        #get the name of the shapefile
        sName = os.path.basename(sFilename_shapefile_in)
        #get the name of the shapefile without extension
        sName_no_extension = os.path.splitext(sName)[0]
        #create a new shapefile
        sFilename_shapefile_out = sFolder + '/' + sName_no_extension + '_transformed.shp'
        #create a new shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')
        #create a new shapefile
        pDataset3 = driver.CreateDataSource(sFilename_shapefile_out)
        #create a new shapefile layer
        layer_out = pDataset3.CreateLayer('boundary', pSpatialRef_target, geom_type=ogr.wkbPolygon)
        #create a new feature
        feature_out = ogr.Feature(layer_out.GetLayerDefn())
        if num_polygons ==1:
            feature = layer.GetNextFeature()
            pPolygon = feature.GetGeometryRef()
            pPolygon.Transform(transform)            
            pass
        else:
            print('The shapefile contains more than one polygon!')
            #merge the polygons
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            # Iterate over all features in the layer
            for feature in layer:
                geometry = feature.GetGeometryRef()                
                if geometry is not None:
                    geometry.Transform(transform)     
                    # Union the geometry of each feature with the merged polygon
                    pPolygon = pPolygon.Union(geometry)
             
        
        #set the geometry
        feature_out.SetGeometry(pPolygon)
        #add the feature to the layer
        layer_out.CreateFeature(feature_out)
        #destroy the feature
        feature_out.Destroy()
        #destroy the feature
        feature = None    
        pDataset3.FlushCache()
        #use the new shapefile to clip the raster
        sFilename_clip = sFilename_shapefile_out

    else:
        sFilename_clip = sFilename_shapefile_in
        #read the first polygon 
        pFeature = layer.GetNextFeature()
        pPolygon = pFeature.GetGeometryRef()
        #get the envelope of the polygon

        iFlag_transform = 0

    #use the gdal warp function to clip the raster
    minX, maxX, minY, maxY = pPolygon.GetEnvelope()
    iNewWidth = int( (maxX - minX) / abs(dPixelWidth)  )
    iNewHeigh = int( (maxY - minY) / abs(dPixelWidth) )
    newGeoTransform = (minX, dPixelWidth, 0,    maxY, 0, -dPixelWidth)  
    if minX > dX_right or maxX < dX_left    or minY > dY_top or maxY < dY_bot:        
        #this polygon is out of bound            
        pass
    else:       
        pDataset_clip = pDriver.Create(sFilename_raster_out, iNewWidth, iNewHeigh, 1, gdalconst.GDT_Float32)
        pDataset_clip.SetGeoTransform( newGeoTransform )
        pDataset_clip.SetProjection( pProjection)   
        pWrapOption = gdal.WarpOptions( cropToCutline=True,cutlineDSName = sFilename_clip , \
                width=iNewWidth,   \
                    height=iNewHeigh,      \
                        dstSRS=pProjection , format = sDriverName )
        pDataset_clip = gdal.Warp(sFilename_raster_out, pDataset_elevation, options=pWrapOption)
        
        #close the dataset
        pDataset_clip = None
        #close the dataset
        pDataset_elevation = None
        print('The raster has been clipped by the shapefile!')
      

    