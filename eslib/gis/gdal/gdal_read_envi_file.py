

import sys
from osgeo import gdal,  gdalconst, osr
def gdal_read_envi_file(sFilename_in):
    '''
	Converts a binary file of ENVI type to a numpy array.
	Lack of an ENVI .hdr file will cause this to crash.
	'''
    pDriverName = gdal.GetDriverByName('ENVI')
    pDriverName.Register()

    pFile = gdal.Open(sFilename_in, gdal.GA_ReadOnly)

    if pFile is None:
        print("Couldn't open this file: " + sFilename_in)
        print('Perhaps you need an ENVI .hdr file?')
        
        sys.exit("Try again!")
    else:
        print("%s opened successfully" % sFilename_in)
        pProjection = pFile.GetProjection()
        print('~~~~~~~~~~~~~~')
        print('Get image size')
        print('~~~~~~~~~~~~~~')
        ncolumn = pFile.RasterXSize
        nrow = pFile.RasterYSize
        nband = pFile.RasterCount

        print("columns: %i" % ncolumn)
        print("nrow: %i" % nrow)
        print("nband: %i" % nband)

        print('~~~~~~~~~~~~~~')
        print('Get georeference information')
        print('~~~~~~~~~~~~~~')
        pGeotransform = pFile.GetGeoTransform()
        dOriginX = pGeotransform[0]
        dOriginY = pGeotransform[3]
        pPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]

        print("origin x: %i" % dOriginX)
        print("origin y: %i" % dOriginY)
        print("width: %2.2f" % pPixelWidth)
        print("height: %2.2f" % pPixelHeight)

        # Set pixel offset.....
        print('~~~~~~~~~~~~~~')
        print('Convert image to 2D array')
        print('~~~~~~~~~~~~~~')
        pBand = pFile.GetRasterBand(1)
        aImage_array = pBand.ReadAsArray(0, 0, ncolumn, nrow)
        image_array_name = sFilename_in
        print(type(aImage_array))
        print(aImage_array.shape)
        pSpatialRef = osr.SpatialReference(wkt=pProjection)
        return aImage_array, pPixelWidth, dOriginX, dOriginY, ncolumn, nrow, pSpatialRef, pGeotransform