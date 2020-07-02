import sys
import numpy as np
from osgeo import gdal, osr,  gdalconst
def gdal_read_envi_file_multiple_band(sFilename_in):
    pDriver = gdal.GetDriverByName('ENVI')
    pDriver.Register()

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

        aData_out = np.full( (nband, nrow, ncolumn) , -9999.0, dtype= float )
        for iBand in range(nband):
            pBand = pDataset.GetRasterBand( iBand + 1)
            dMissing_value = pBand.GetNoDataValue()
            aData_out[iBand, :, :] = pBand.ReadAsArray(0, 0, ncolumn, nrow)

        pSpatialRef = osr.SpatialReference(wkt=pProjection)

        pDriver = None
        pDataset = None
        pBand = None

        return aData_out, pPixelWidth, dOriginX, dOriginY, nband, nrow, ncolumn, dMissing_value,  pSpatialRef, pGeotransform