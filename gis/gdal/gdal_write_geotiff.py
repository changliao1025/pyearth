
import sys
from osgeo import gdal, osr
def gdal_write_geotiff(sFilename_in, ncolumn_in, nrow_in, dX_origin_in, dY_origin_in, dPixelWidth_in, dMissing_value_in,\
     aData_image_in, pSpatialRef_in):

    """Array > Raster
    Save a raster from a C order aData_image_in.

    :param aData_image_in: ndarray    
    """
    pDriver = gdal.GetDriverByName('GTiff')

    pDataset = pDriver.Create(
        sFilename_in,
        ncolumn_in,
        nrow_in,
        1,
        gdal.GDT_Float32, )

    pDataset.SetGeoTransform((
        dX_origin_in,    # 0
        dPixelWidth_in,  # 1
        0,                      # 2
        dY_origin_in,    # 3
        0,                      # 4
        -dPixelWidth_in))  
    pProjection= pSpatialRef_in.ExportToPrettyWkt()
    pDataset.SetProjection(pProjection)
    pDataset.GetRasterBand(1).WriteArray(aData_image_in)
    #pDataset.GetRasterBand(1).
    pDataset.FlushCache()  # Write to disk.
    pDriver = None
    pDataset = None
    #return pDataset, pDataset.GetRasterBand(1) 