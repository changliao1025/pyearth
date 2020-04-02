import sys
from osgeo import gdal, osr, gdalconst
def gdal_write_geotiff(sFilename_in, aData_in,\
                       dPixelWidth_in, \
                       dOriginX_in, \
                       dOriginY_in, \
                       dMissing_value_in,\
                       pSpatialRef_in):

    pDriver = gdal.GetDriverByName('GTiff')
    pDriver.Register()
    nrow, ncolumn = aData_in.shape
    nband = 1

    pDataset = pDriver.Create(
        sFilename_in,
        ncolumn,
        nrow,
        nband,
        gdal.GDT_Float32, )

    pDataset.SetGeoTransform([
        dOriginX_in,    # 0
        dPixelWidth_in,  # 1
        0,                      # 2
        dOriginY_in,    # 3
        0,                      # 4
        -dPixelWidth_in])

    pProjection= pSpatialRef_in.ExportToPrettyWkt()
    pDataset.SetProjection(pProjection)

    pBand = pDataset.GetRasterBand(1)
    pBand.WriteArray(aData_in)
    pBand.SetNoDataValue(dMissing_value_in)

    pGeotransform_out = pDataset.GetGeoTransform()
    pProjection_out = pDataset.GetProjection()

    pDataset.FlushCache()  # Write to disk.
    pDriver = None
    pDataset = None
    pBand=None
    return  sFilename_in, pGeotransform_out, pProjection_out
