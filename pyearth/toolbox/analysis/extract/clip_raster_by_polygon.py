import os, sys
import numpy as np
from osgeo import gdal, ogr
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.management.vector.reproject import reproject_vector
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file


def clip_raster_by_polygon_file(sFilename_raster_in, sFilename_polygon_in, sFilename_raster_out, sFormat_in='GTiff'):
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
    
    if os.path.exists(sFilename_polygon_in):
        pass
    else:   
        print('The polygon file does not exist!')
        return
    
    #get the extension of raster file
    if sFormat_in is not None:
        sDriverName = sFormat_in
    else: 
        sDriverName = 'GTiff'
        sExtension_raster = os.path.splitext(sFilename_raster_in)[1]
    #get the driver for the extension    
    pDriver_raster = gdal.GetDriverByName(sDriverName)

    #get the extension of polygon file  
    sExtension_vector = os.path.splitext(sFilename_polygon_in)[1]
    #get the driver for the extension
    if sExtension_vector == '.geojson':
        pDriver_vector = ogr.GetDriverByName('GeoJSON')
    else:
        pDriver_vector = ogr.GetDriverByName(sExtension_vector)
   

    pDataset_data = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)   
    dummy = gdal_read_geotiff_file(sFilename_raster_in)
    aData= dummy['dataOut']    
    eType = dummy['dataType']  
    dPixelWidth = dummy['pixelWidth']                        
    pPixelHeight = dummy['pixelHeight']
    dOriginX = dummy['originX']
    dOriginY = dummy['originY']
    nrow = dummy['nrow']
    ncolumn = dummy['ncolumn']
    dMissing_value= dummy['missingValue']
    pProjection = dummy['projection']
    pSpatial_reference_target = dummy['spatialReference']
    pProjection_target = pSpatial_reference_target.ExportToWkt()
    dX_left=dOriginX
    dX_right = dOriginX + ncolumn * dPixelWidth
    dY_top = dOriginY
    dY_bot = dOriginY + nrow * pPixelHeight

    #get the spatial reference of the shapefile
    pDataset_clip = ogr.Open(sFilename_polygon_in)    
    # Get the first layer in the shapefile
    pLayer_clip = pDataset_clip.GetLayer(0)   
    # Count the number of features (polygons)
    nPolygon = pLayer_clip.GetFeatureCount()
    if nPolygon == 0:
        print('The shapefile does not contain any polygon!')
        return
    else:
        if nPolygon > 1:
            pDataset_clip = None
            pLayer_clip = None
            print('The polygon contains more than one polygon, the program will attempt to merge them as one!')
            #obtain the file extension
            sFilename_clip_new = sFilename_polygon_in.replace(sExtension_vector, '_merged' + sExtension_vector)
            merge_features(sFilename_polygon_in, sFilename_clip_new)
            sFilename_polygon_in = sFilename_clip_new
            #open the new file 
            pDataset_clip = ogr.Open(sFilename_polygon_in)    
            # Get the first layer in the shapefile
            pLayer_clip = pDataset_clip.GetLayer(0) 
            pass

    # Get the spatial reference of the layer
    pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
    pProjection_clip = pSpatial_reference_clip.ExportToWkt()
    print(pProjection_clip)   

    if( pProjection_target != pProjection_clip):        
        pDataset_clip = None
        pLayer_clip = None
        #in this case, we can reproject the shapefile to the same spatial reference as the raster
        #get the folder that contains the shapefile
        sFolder = os.path.dirname(sFilename_polygon_in)
        #get the name of the shapefile
        sName = os.path.basename(sFilename_polygon_in)
        #get the name of the shapefile without extension
        sName_no_extension = os.path.splitext(sName)[0]
        #create a new shapefile
        sFilename_clip_out = sFolder + '/' + sName_no_extension + '_transformed' + sExtension_vector
        reproject_vector(sFilename_polygon_in, sFilename_clip_out, pProjection_target)        
        #use the new shapefile to clip the raster        
        sFilename_clip = sFilename_clip_out
        pDataset_clip = ogr.Open(sFilename_clip_out)
        pLayer_clip = pDataset_clip.GetLayer(0) 
        pFeature_clip = pLayer_clip.GetNextFeature()
        pPolygon = pFeature_clip.GetGeometryRef() 

    else:
        sFilename_clip = sFilename_polygon_in              
        #read the first polygon 
        pFeature_clip = pLayer_clip.GetNextFeature()
        pPolygon = pFeature_clip.GetGeometryRef() 
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
        pDataset_clip = pDriver_raster.Create(sFilename_raster_out, iNewWidth, iNewHeigh, 1, eType)
        pDataset_clip.SetGeoTransform( newGeoTransform )
        pDataset_clip.SetProjection( pProjection)   
        pWrapOption = gdal.WarpOptions( cropToCutline=True,cutlineDSName = sFilename_clip , 
                width=iNewWidth,   
                    height=iNewHeigh,      
                        dstSRS=pSpatial_reference_target , format = 'MEM' )
        pDataset_clip_warped = gdal.Warp('', pDataset_data, options=pWrapOption)

        #convert the warped dataset to an array
        aData_clip = pDataset_clip_warped.ReadAsArray()
        #change the missing value to the original missing value
        aData_clip[aData_clip == dMissing_value] = -9999
        # Write the warped dataset to the output raster
        pDataset_clip.GetRasterBand(1).WriteArray(aData_clip)        
        #close the dataset
        pDataset_clip = None
        #close the dataset
        pDataset_data = None
        print('The raster has been clipped by the polygon file!')

    pSpatial_reference_target = None
    pSpatial_reference_clip = None
  
    

      

    