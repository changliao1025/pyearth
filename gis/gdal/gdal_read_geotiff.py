
import sys
from osgeo import gdal, osr
def gdal_read_geotiff(sFilename_in):
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
        #type(pFile)
        # Projection
        pSpatialRef = pFile.GetProjection()

        # Dimensions       

        # Metadata for the pFile dataset
        pFile.GetMetadata()
        print('~~~~~~~~~~~~~~')
        print('Get image size')
        print('~~~~~~~~~~~~~~')
        ncolumn = pFile.RasterXSize
        nrow = pFile.RasterYSize
        nband = pFile.RasterCount

        #print("columns: %i" % ncolumn)
        #print("nrow: %i" % nrow)
        #print("nband: %i" % nband)

        print('~~~~~~~~~~~~~~')
        print('Get georeference information')
        print('~~~~~~~~~~~~~~')
        pGeotransform = pFile.GetGeoTransform()
        dX_origin = pGeotransform[0]
        dY_origin = pGeotransform[3]
        dPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]

        print("origin x: %i" % dX_origin)
        print("origin y: %i" % dY_origin)
        print("width: %2.2f" % dPixelWidth)
        print("height: %2.2f" % pPixelHeight)

        # Set pixel offset.....
        print('~~~~~~~~~~~~~~')
        print('Convert image to 2D array')
        print('~~~~~~~~~~~~~~')
        
        pBand = pFile.GetRasterBand(1)

        # Check type of the variable 'pBand'
        #type(pBand)

        # Data type of the values
        gdal.GetDataTypeName(pBand.DataType)
        # Compute statistics if needed
        if pBand.GetMinimum() is None or pBand.GetMaximum()is None:
            pBand.ComputeStatistics(0)
            print("Statistics computed.")

        # Fetch metadata for the pBand
        

        # Print only selected metadata:
        
        dMissing_value = pBand.GetNoDataValue()
        print ("[ NO DATA VALUE ] = ", dMissing_value) # none
        print ("[ MIN ] = ", pBand.GetMinimum())
        print ("[ MAX ] = ", pBand.GetMaximum())
        aData_image = pBand.ReadAsArray(0, 0, ncolumn, nrow)
        sFilename_image = sFilename_in
        #print(type(aData_image))
        #print(aData_image.shape)

        #there are different types of spatial reference systems
        #we will use one of them to keep the consistency
        pSpatialRef = osr.SpatialReference(wkt=pSpatialRef)
        #there are many information in a raster data, we will use some standard way to output them
        #beblow is an example for ArcGIS ASCII file
        #NCOLS xxx
        #NROWS xxx
        #XLLCENTER xxx
        #YLLCENTER xxx
        #CELLSIZE xxx
        #NODATA_VALUE xxx
        #row 1
        #row 2
        #...
        #row n
        return ncolumn, nrow, dX_origin, dY_origin, dPixelWidth, dMissing_value, aData_image, pSpatialRef