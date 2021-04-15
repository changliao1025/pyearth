import os, sys
import numpy as np
from osgeo import gdal, osr

def gdal_read_geotiff_file(sFilename_in):
    pDriver = gdal.GetDriverByName('GTiff')
    pDriver.Register()
    pDataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)

    if pDataset is None:
        print("Couldn't open this file: " + sFilename_in)
        sys.exit("Try again!")
    else:       
        pProjection = pDataset.GetProjection()

        pDataset.GetMetadata()
       
        ncolumn = pDataset.RasterXSize
        nrow = pDataset.RasterYSize
        nband = pDataset.RasterCount

        pGeotransform = pDataset.GetGeoTransform()
        dOriginX = pGeotransform[0]
        dOriginY = pGeotransform[3]
        dPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]

        pBand = pDataset.GetRasterBand(1)

        # Data type of the values
        gdal.GetDataTypeName(pBand.DataType)
        # Compute statistics if needed
        if pBand.GetMinimum() is None or pBand.GetMaximum()is None:
            pBand.ComputeStatistics(0)

        dMissing_value = pBand.GetNoDataValue()
       
        aData_out = pBand.ReadAsArray(0, 0, ncolumn, nrow)
    
        #we will use one of them to keep the consistency
        pSpatialRef = osr.SpatialReference(wkt=pProjection)
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
        #close file

        pDataset = None
        pBand = None
        
        

        return aData_out, dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value,  pSpatialRef, pGeotransform


def gdal_read_geotiff_file_multiple_band(sFilename_in):
    pDriver = gdal.GetDriverByName('GTiff')
    pDriver.Register()
    

    pDataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)

    if pDataset is None:
        print("Couldn't open this file: " + sFilename_in)
        print('Perhaps you need an ENVI .hdr file?')        
        sys.exit("Try again!")
    else:        
        pProjection = pDataset.GetProjection()
        pDataset.GetMetadata()
       
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
        #close file

        pDriver = None
        pDataset = None
        pBand = None
        return aData_out, pPixelWidth, dOriginX, dOriginY, nband, nrow, ncolumn, dMissing_value,  pSpatialRef, pGeotransform