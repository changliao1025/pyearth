
import sys
from osgeo import gdal,osr
def gdal_read_tiff(sFilename_in):
    '''
	Converts a binary file of ENVI type to a numpy array.
	Lack of an ENVI .hdr file will cause this to crash.
	'''
    

    pFile = gdal.Open(sFilename_in, gdal.GA_ReadOnly)

    if pFile is None:
        print("Couldn't open this file: " + sFilename_in)
        print('Perhaps you need an ENVI .hdr file?')
        
        sys.exit("Try again!")
    else:
        # Check type of the variable 'pFile'
        print("%s opened successfully" % sFilename_in)
        type(pFile)
        # Projection
        proj = pFile.GetProjection()

        # Dimensions       

        # Metadata for the pFile dataset
        pFile.GetMetadata()
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

        # Check type of the variable 'band'
        type(band)

        # Data type of the values
        gdal.GetDataTypeName(band.DataType)
        # Compute statistics if needed
        if band.GetMinimum() is None or band.GetMaximum()is None:
            band.ComputeStatistics(0)
            print("Statistics computed.")

        # Fetch metadata for the band
        

        # Print only selected metadata:
        print ("[ NO DATA VALUE ] = ", band.GetNoDataValue()) # none
        print ("[ MIN ] = ", band.GetMinimum())
        print ("[ MAX ] = ", band.GetMaximum())
        image_array = band.ReadAsArray(0, 0, cols, rows)
        image_array_name = sFilename_in
        print(type(image_array))
        print(image_array.shape)

        
        proj = osr.SpatialReference(wkt=proj)
        return image_array, pixelWidth,originX, originY, cols,rows,proj, geotransform