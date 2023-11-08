
import os, sys

from osgeo import gdal, osr


def gdal_read_ascii_file(sFilename_in):
    """Read a ENVI standard format raster file.

    Args:
        sFilename_in (string): The filename

    Returns:
        tuple: aData_out, pPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value , pGeotransform, pProjection,  pSpatial_reference
    """
    
    if os.path.exists(sFilename_in):
        pass
    else:
        print('The file does not exist!')
        return

    sDriverName='AAIGrid'
    pDriver = gdal.GetDriverByName(sDriverName)  

    if pDriver is None:
        print ("%s pDriver not available.\n" % sDriverName)
    else:
        print  ("%s pDriver IS available.\n" % sDriverName)
         
    pDataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)

    if pDataset is None:
        print("Couldn't open this file: " + sFilename_in)
        print('Perhaps you need an ENVI .hdr file?')
        
        sys.exit("Try again!")
    else:
        
        pProjection = pDataset.GetProjection()
       
        ncolumn = pDataset.RasterXSize
        nrow = pDataset.RasterYSize
        nband = pDataset.RasterCount

        pGeotransform = pDataset.GetGeoTransform()
        dOriginX = pGeotransform[0]
        dOriginY = pGeotransform[3]
        pPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]
       
        pBand = pDataset.GetRasterBand(1)
        dMissing_value = pBand.GetNoDataValue()
        aData_out = pBand.ReadAsArray(0, 0, ncolumn, nrow)
        pSpatial_reference = osr.SpatialReference(wkt=pProjection)

        pDriver = None
        pDataset = None
        pBand = None

        return aData_out, pPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value , pGeotransform, pProjection,  pSpatial_reference
