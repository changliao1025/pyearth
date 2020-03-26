import sys
from osgeo import gdal, osr, gdalconst

def gdal_write_geotiff_multiple_band(sFilename_in, aData_in, \
                                     dPixelWidth_in,\
                                     dOriginX_in, \
                                     dOriginY_in,  \
                                     dMissing_value_in,\
                                     pSpatialRef_in ):

    pDriver = gdal.GetDriverByName('GTiff')
    pDriver.Register()
    nband, nrow, ncolumn = aData_in.shape


    # Creates a new raster data source
    pDataset = pDriver.Create(sFilename_in, \
                              ncolumn, \
                              nrow, \
                              nband,\
                              gdal.GDT_Float32)

    # Write metadata

    pDataset.SetGeoTransform([dOriginX_in, \
                              dPixelWidth_in, \
                              0.0, \
                              dOriginY_in, \
                              0.0, \
                              -dPixelWidth_in])

    pProjection= pSpatialRef_in.ExportToPrettyWkt()
    pDataset.SetProjection(pProjection)

    #Write raster datasets
    for iBand in range(nband):
        pBand = pDataset.GetRasterBand(iBand + 1)
        pBand.WriteArray(aData_in[iBand, :, :])
        pBand.SetNoDataValue(dMissing_value_in)

    pGeotransform_out = pDataset.GetGeoTransform()
    pProjection_out = pDataset.GetProjection()

    pDataset.FlushCache()  # Write to disk.
    pDriver =None
    pDataset=None
    pBand=None

    return  sFilename_in, pGeotransform_out, pProjection_out
