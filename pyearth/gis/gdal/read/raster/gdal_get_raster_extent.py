
from osgeo import gdal, osr
def gdal_get_raster_extent(sFilename_in):
    ds = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    gt = ds.GetGeoTransform()
    min_x = gt[0]
    max_x = gt[0] + gt[1]*ds.RasterXSize
    min_y = gt[3] + gt[5]*ds.RasterYSize
    max_y = gt[3]

    #check the srs of the raster
    pProjection = ds.GetProjection()

    pSpatialRef_WGS84 = osr.SpatialReference()
    pSpatialRef_WGS84.ImportFromEPSG(4326)
    pProjection_WGS84 = pSpatialRef_WGS84.ExportToWkt()

    #check whether it is a lat/lon projection
    if pProjection == pProjection_WGS84:
        #check whether the extent is in the range of -180, 180, -90, 90
        if min_x < -180:
            min_x = -180
        if max_x > 180:
            max_x = 180
        if min_y < -90:
            min_y = -90
        if max_y > 90:
            max_y = 90
    else:
        pass

    ds = None

    return min_x, max_x, min_y, max_y

