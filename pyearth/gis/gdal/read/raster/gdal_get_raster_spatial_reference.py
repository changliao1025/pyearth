import os
from osgeo import gdal

def gdal_get_raster_spatial_reference_wkt(sFilename_in):
    ds = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    pProjection =ds.GetProjection()
    return pProjection