#we will use gdal api for most operations
import os, sys
import numpy as np
from osgeo import ogr, osr, gdal
from pyearth.toolbox.conversion.rasterize_vector import rasterize_vector

def convert_vector_to_global_raster(sFilename_vector_in, sFilename_tif_out,
                                       dResolution_x_in, dResolution_y_in,
                                       iFlag_boundary_only_in = 0, dFill_value_in = 2 ):

    #get the folder that contains the output file
    sFolder = os.path.dirname(sFilename_tif_out)
    if not os.path.exists(sFolder):
        os.makedirs(sFolder)

    if os.path.exists(sFilename_tif_out):
        #check file size and remove if it is empty
        if os.path.getsize(sFilename_tif_out) == 0:
            print('remove empty file: ' + sFilename_tif_out)
            os.remove(sFilename_tif_out)
        else:
            os.remove(sFilename_tif_out)

    pDataSource_clip = ogr.Open(sFilename_vector_in)
    if pDataSource_clip is None:
        print("Could not open clip polygon")
        return

    pLayer_clip = pDataSource_clip.GetLayer()
    dLon_min, dLon_max, dLat_min, dLat_max = pLayer_clip.GetExtent()

    dLon_min = -180
    dLon_max = 180
    dLat_min = -90
    dLat_max = 90

    nleft  = np.floor(  (dLon_min - (-180)) /(dResolution_x_in)  )
    nright = np.ceil(  (dLon_max - (-180)) /(dResolution_x_in)  )
    ntop  = np.floor(  (90 - dLat_max) /(dResolution_y_in)  )
    nbot = np.ceil(  (90 - dLat_min) /(dResolution_y_in)  )
    dMin_x = -180 + nleft*dResolution_x_in
    dMax_x = -180 + nright*dResolution_x_in
    dMin_y = 90 - nbot*dResolution_y_in
    dMax_y = 90 - ntop*dResolution_y_in
    if dMin_x < -180:
        dMin_x = -180
    if dMax_x > 180:
        dMax_x = 180
    if dMin_y < -90:
        dMin_y = -90
    if dMax_y > 90:
        dMax_y = 90

    nrow = int(nbot-ntop)
    ncolumn = int(nright - nleft)

    rasterize_vector(sFilename_vector_in, sFilename_tif_out,
                      dResolution_x_in, dResolution_y_in,
                         dMissing_value_in=0,
                         iDataType_out = gdal.GDT_Byte ,
                         dMin_x_in = dMin_x, dMax_x_in = dMax_x, dMin_y_in = dMin_y, dMax_y_in = dMax_y,
                         nRow_in=nrow, nColumn_in=ncolumn,
                             iFlag_boundary_only_in = iFlag_boundary_only_in,
                             dFill_value_in = dFill_value_in )
    return