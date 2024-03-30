import os
import numpy as np
from osgeo import gdal, ogr, osr, gdalconst

def gdal_get_raster_spatial_ref_wkt(sFilename_in):
    ds = gdal.Open(sFilename_in, gdal.GA_ReadOnly)    
    pProjection =ds.GetProjection()
    return pProjection