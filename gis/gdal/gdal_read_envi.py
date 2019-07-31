

import sys
from osgeo import gdal,  gdalconst
def gdal_read_envi(sFilename_in):
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
        pSpatialRef = pFile.GetProjection()
        print('~~~~~~~~~~~~~~')
        print('Get image size')
        print('~~~~~~~~~~~~~~')
        cols = pFile.RasterXSize
        rows = pFile.RasterYSize
        bands = pFile.RasterCount

        print("columns: %i" % cols)
        print("rows: %i" % rows)
        print("bands: %i" % bands)

        print('~~~~~~~~~~~~~~')
        print('Get georeference information')
        print('~~~~~~~~~~~~~~')
        geotransform = pFile.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]

        print("origin x: %i" % originX)
        print("origin y: %i" % originY)
        print("width: %2.2f" % pixelWidth)
        print("height: %2.2f" % pixelHeight)

        # Set pixel offset.....
        print('~~~~~~~~~~~~~~')
        print('Convert image to 2D array')
        print('~~~~~~~~~~~~~~')
        band = pFile.GetRasterBand(1)
        image_array = band.ReadAsArray(0, 0, cols, rows)
        image_array_name = sFilename_in
        print(type(image_array))
        print(image_array.shape)
        pSpatialRef = osr.SpatialReference(wkt=pSpatialRef)
        return image_array, pixelWidth, originX, originY, cols, rows, pSpatialRef, geotransform