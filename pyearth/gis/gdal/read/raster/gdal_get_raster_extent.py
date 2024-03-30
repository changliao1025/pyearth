
from osgeo import gdal, ogr, osr, gdalconst
def gdal_get_raster_extent(sFilename_in):
    ds = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    gt = ds.GetGeoTransform()
    min_x = gt[0]
    max_x = gt[0] + gt[1]*ds.RasterXSize
    min_y = gt[3] + gt[5]*ds.RasterYSize
    max_y = gt[3]
    ds = None
    return min_x, max_x, min_y, max_y