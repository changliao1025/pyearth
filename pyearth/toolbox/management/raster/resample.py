import os
import numpy as np
from osgeo import gdal, osr, gdalconst

from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.gdal.read.raster.gdal_get_raster_spatial_ref import gdal_get_raster_spatial_ref_wkt
from pyearth.gis.gdal.gdal_to_numpy_datatype import gdal_to_numpy_datatype
gdal.UseExceptions()    
#set default spatial reference
pSpatialRef_default = osr.SpatialReference()
pSpatialRef_default.ImportFromEPSG(4326)  # WGS84
sProjection_default = pSpatialRef_default.ExportToWkt()
pSpatialRef_default = None

def resample_raster(sFilename_in, sFilename_out, dResolution_x, dResolution_y,
                     pProjection_target_in = None, iData_type = gdalconst.GDT_Int16, 
                     sResampleAlg = 'MODE', dMissing_value_source= 255, dMissing_value_target = -9999):

    #check if the input raster exists
    if not os.path.exists(sFilename_in):
        print('Error: the input raster does not exist')
        return
    
    pDriver_tiff = gdal.GetDriverByName('GTiff')
    
    #get the spatial from the input raster
    pProjection_source = gdal_get_raster_spatial_ref_wkt(sFilename_in)

    #check whether it is wgs84 or not using project string
    pSpatialRef_source = osr.SpatialReference()
    pSpatialRef_source.ImportFromWkt(pProjection_source)
    if pProjection_source == sProjection_default:
        iFlag_wgs84_source = 1
    else:
        iFlag_wgs84_source = 0

    if pProjection_target_in is None:
        pProjection_target = pProjection_source
    else:
        pProjection_target = pProjection_target_in  

    
    if pProjection_target == sProjection_default:
        iFlag_wgs84_target = 1
    else:   
        iFlag_wgs84_target = 0
    
    if os.path.exists(sFilename_out):
        #check file size and remove if it is empty
        if os.path.getsize(sFilename_out) == 0:
            print('remove empty file: ' + sFilename_out)
            os.remove(sFilename_out) 
    else:
        print('create file: ' + sFilename_out)

    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(pProjection_target)
    #use the same approach gire to define the extent
    #get the extent of raster geotiff file
    
    pDataset_in = gdal.Open(sFilename_in, gdal.GA_ReadOnly)

    if iFlag_wgs84_source == 1:       #from wgs84
        dLon_min, dLon_max, dLat_min, dLat_max = gdal_get_raster_extent(sFilename_in)
        if iFlag_wgs84_target == 1: #to wgs84
            #if the target is wgs84, we can use the same approach
            nleft  = np.floor(  (dLon_min - (-180)) /(dResolution_x)  )
            nright = np.ceil(  (dLon_max - (-180)) /(dResolution_x)  )
            ntop  = np.floor(  (90 - dLat_max) /(dResolution_y)  )
            nbot = np.ceil(  (90 - dLat_min) /(dResolution_y)  )
            nrow = int(nbot-ntop)
            ncolumn = int(nright - nleft)

            dMin_x = -180 + nleft*dResolution_x
            dMax_y = 90 - ntop*dResolution_y
        else: #to something else
            #if the target is not wgs84, we need to reproject the extent
            pTransform = osr.CoordinateTransformation(pSpatialRef_source, pSpatialRef_target)
            dX_min, dY_min, _ = pTransform.TransformPoint(dLon_min, dLat_min)
            dX_max, dY_max, _ = pTransform.TransformPoint(dLon_max, dLat_max)
            nrow = int((dY_max - dY_min)/dResolution_y)
            ncolumn = int((dX_max - dX_min)/dResolution_x)
            dMin_x = dX_min
            dMax_y = dY_max
    else:
        #if it is not wgs84, it has a projection already
        #get the extent of raster geotiff file
        dX_min, dX_max, dY_min, dY_max = gdal_get_raster_extent(sFilename_in)
        #if the projection is not the same, we need to reproject the extent
        if iFlag_wgs84_target == 1: #to wgs84
            pTransform = osr.CoordinateTransformation(pSpatialRef_source, pSpatialRef_target)
            dLon_min, dLat_min, _ = pTransform.TransformPoint(dX_min, dY_min)
            dLon_max, dLat_max, _ = pTransform.TransformPoint(dX_max, dY_max)
            nleft  = np.floor(  (dLon_min - (-180)) /(dResolution_x)  )
            nright = np.ceil(  (dLon_max - (-180)) /(dResolution_x)  )
            ntop  = np.floor(  (90 - dLat_max) /(dResolution_y)  )
            nbot = np.ceil(  (90 - dLat_min) /(dResolution_y)  )
            nrow = int(nbot-ntop)
            ncolumn = int(nright - nleft)
            dMin_x = -180 + nleft*dResolution_x
            dMax_y = 90 - ntop*dResolution_y
        else: #to something else, but could be the same or different
            if pProjection_source == pProjection_target:
                nrow = int((dY_max - dY_min)/dResolution_y)  
                ncolumn = int((dX_max - dX_min)/dResolution_x)   
                dMin_x = dX_min
                dMax_y = dY_max
            else: #not the same
                pTransform = osr.CoordinateTransformation(pSpatialRef_source, pSpatialRef_target)
                dX_min, dY_min, _ = pTransform.TransformPoint(dX_min, dY_min)
                dX_max, dY_max, _ = pTransform.TransformPoint(dX_max, dY_max)
                nrow = int((dY_max - dY_min)/dResolution_y) 
                ncolumn = int((dX_max - dX_min)/dResolution_x)    
                dMin_x = dX_min
                dMax_y = dY_max            

    iNewWidth = ncolumn
    iNewHeigh = nrow

    newGeoTransform = (dMin_x, dResolution_x, 0, dMax_y, 0, -dResolution_y)  
    pDataset_clip = pDriver_tiff.Create(sFilename_out, iNewWidth, iNewHeigh, 1, iData_type) #this data type may be provided as an input argument
    pDataset_clip.SetGeoTransform( newGeoTransform )
    pDataset_clip.SetProjection( pProjection_target) 
    
    pWrapOption = gdal.WarpOptions( cropToCutline=False, #could be true if vector file is provided
                                width=iNewWidth,   
                                    height=iNewHeigh,      
                                        dstSRS=pSpatialRef_target , format = 'MEM',
                                        resampleAlg=sResampleAlg ) #this resample algorithm may be provided as an input argument

    pDataset_clip_warped = gdal.Warp('', pDataset_in, options=pWrapOption)#gdal.Warp(sFilename_out, pDataset_in, options=pWrapOption)

    #convert the warped dataset to an array
    aData_clip = pDataset_clip_warped.ReadAsArray()
    #change the gdal data type to numpy data type
    iData_type_numpy = gdal_to_numpy_datatype(iData_type)

    aData_clip = aData_clip.astype(iData_type_numpy)
    aData_clip[aData_clip == dMissing_value_source] = dMissing_value_target
    pDataset_clip.GetRasterBand(1).WriteArray(aData_clip)        
    pDataset_clip = None
    pSpatialRef_source = None
    pSpatialRef_target = None
    return