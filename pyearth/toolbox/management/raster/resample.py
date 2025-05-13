import os
import numpy as np
from osgeo import gdal, osr, gdalconst

from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.gdal.read.raster.gdal_get_raster_spatial_reference import gdal_get_raster_spatial_reference_wkt
from pyearth.gis.gdal.gdal_to_numpy_datatype import gdal_to_numpy_datatype
gdal.UseExceptions()
#set default spatial reference


def resample_raster(sFilename_in, sFilename_out, dResolution_x, dResolution_y,
                     pProjection_target_in = None, iData_type = gdalconst.GDT_Int16,
                     sResampleAlg = 'MODE', dMissing_value_source= 255, dMissing_value_target = -9999):

    pSpatialRef_default = osr.SpatialReference()
    pSpatialRef_default.ImportFromEPSG(4326)  # WGS84
    sProjection_default = pSpatialRef_default.ExportToWkt()
    pSpatialRef_default = None
    #check if the input raster exists


    if not os.path.exists(sFilename_in):
        print('Error: the input raster does not exist')
        return

    pDriver_tiff = gdal.GetDriverByName('GTiff')

    #get the spatial from the input raster
    pProjection_source = gdal_get_raster_spatial_reference_wkt(sFilename_in)

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
    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']
    pWrapOption = gdal.WarpOptions( cropToCutline=False, #could be true if vector file is provided
                                     xRes=dResolution_x,
                                     yRes=dResolution_y,
                                        dstSRS=pSpatialRef_target , format = 'MEM',
                                        resampleAlg=sResampleAlg ) #this resample algorithm may be provided as an input argument

    pDataset_clip_warped = gdal.Warp('', pDataset_in, options=pWrapOption)#gdal.Warp(sFilename_out, pDataset_in, options=pWrapOption)
    newGeoTransform = pDataset_clip_warped.GetGeoTransform()
    #convert the warped dataset to an array
    aData_clip = pDataset_clip_warped.ReadAsArray()
    iNewWidth = aData_clip.shape[1]
    iNewHeigh = aData_clip.shape[0]
    pDataset_clip = pDriver_tiff.Create(sFilename_out, iNewWidth, iNewHeigh, 1, iData_type, options= options) #this data type may be provided as an input argument
    pDataset_clip.SetGeoTransform( newGeoTransform )
    pDataset_clip.SetProjection( pProjection_target)
    pDataset_clip.GetRasterBand(1).SetNoDataValue(dMissing_value_target)
    #change the gdal data type to numpy data type
    iData_type_numpy = gdal_to_numpy_datatype(iData_type)
    aData_clip = aData_clip.astype(iData_type_numpy)
    aData_clip[aData_clip == dMissing_value_source] = dMissing_value_target
    pDataset_clip.GetRasterBand(1).WriteArray(aData_clip)
    pDataset_clip.GetRasterBand(1).FlushCache()  # Corrected method name to FlushCache()
    pDataset_clip.FlushCache()
    pDataset_clip = None
    pSpatialRef_source = None
    pSpatialRef_target = None
    return