
import os
import numpy as np
from osgeo import gdal, osr, gdalconst
from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.gdal.gdal_to_numpy_datatype import gdal_to_numpy_datatype
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file

def merge_rasters(aFilename_rasters, sFilename_merge_out, dResolution_x, dResolution_y,
                           sFilename_mask,sResampleAlg = 'NEAREST',
                          dMissing_value_source=-9999, dMissing_value_target=-9999, iData_type = gdal.GDT_Int16):


    pDriver_tiff = gdal.GetDriverByName('GTiff')
    # List of raster files to merge
    iData_type_numpy = gdal_to_numpy_datatype( iData_type )
    # Output file

    if os.path.exists(sFilename_merge_out):
        os.remove(sFilename_merge_out)
    #set up the width et al
    #how to set the extent?
    #we can use the global land-ocean mask as extent?
    #get the spatial extent
    dLon_min, dLon_max, dLat_min, dLat_max = gdal_get_raster_extent(sFilename_mask)
    #nleft  = np.floor(  (dLon_min - (-180)) /(dResolution_x)  )
    #nright = np.ceil(  (dLon_max - (-180)) /(dResolution_x)  )
    #ntop  = np.ceil(  (90 - dLat_max) /(dResolution_y)  )
    #nbot = np.ceil(  (90 - dLat_min) /(dResolution_y)  )
    #nrow = int(nbot-ntop)
    #ncolumn = int(nright - nleft)

    dummy = gdal_read_geotiff_file(sFilename_mask)
    pGeoTransform = dummy['geotransform']
    pProjection_target = dummy['projection']

    ncolumn = dummy['ncolumn']
    nrow = dummy['nrow']

    iNewWidth = ncolumn
    iNewHeigh = nrow

    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(pProjection_target)
    #pSpatialRef_target.ImportFromEPSG(4326)  # WGS84
    #pProjection_target = pSpatialRef_target.ExportToWkt()

    #newGeoTransform = (dLon_min, dResolution_x, 0, dLat_max, 0, -dResolution_y)
    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']
    pDataset_clip = pDriver_tiff.Create(sFilename_merge_out, iNewWidth, iNewHeigh, 1, eType=iData_type, options = options) #this data type may be provided as an input argument
    pDataset_clip.SetGeoTransform( pGeoTransform )
    pDataset_clip.SetProjection( pProjection_target)
    pDataset_clip.GetRasterBand(1).SetNoDataValue(dMissing_value_target)

    pWrapOption = gdal.WarpOptions( cropToCutline=False, #could be true if vector file is provided
                                    width=iNewWidth,
                                        height=iNewHeigh,
                                        outputBounds = [dLon_min, dLat_min, dLon_max, dLat_max],
                                            dstSRS=pSpatialRef_target , format = 'MEM',
                                            resampleAlg=sResampleAlg ) #this resample algorithm may be provided as an input argument
    # Use gdal.Warp to merge the raster files
    pDataset_clip_warped = gdal.Warp('', aFilename_rasters, options=pWrapOption)#gdal.Warp(sFilename_merge, aFilename_river_network, options=pWrapOption)

    #convert the warped dataset to an array
    aData_clip = pDataset_clip_warped.ReadAsArray()
    #change the gdal data type to numpy data type
    aData_clip[np.where(aData_clip == dMissing_value_source)] = dMissing_value_target
    aData_clip = aData_clip.astype(iData_type_numpy)
    pDataset_clip.GetRasterBand(1).WriteArray(aData_clip)
    pDataset_clip.GetRasterBand(1).FlushCache()
    pDataset_clip = None
    pSpatialRef_target = None
    return

