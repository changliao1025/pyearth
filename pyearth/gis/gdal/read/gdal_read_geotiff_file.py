import os, sys
import numpy as np
from osgeo import gdal, osr,ogr, gdalconst

def gdal_read_geotiff_file(sFilename_in):
    """Read a Geotiff raster file."""
    sDriverName='GTiff'
    pDriver = gdal.GetDriverByName(sDriverName)  

    if pDriver is None:
        print ("%s pDriver not available.\n" % sDriverName)
    else:
        print  ("%s pDriver IS available.\n" % sDriverName)  
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
        if pBand.GetMinimum() is None or pBand.GetMaximum() is None:
            pBand.ComputeStatistics(0)

        dMissing_value = pBand.GetNoDataValue()
       
        aData_out = pBand.ReadAsArray(0, 0, ncolumn, nrow)
    
        #we will use one of them to keep the consistency
        pSpatialRef = osr.SpatialReference(wkt=pProjection)
       

        pDataset = None
        pBand = None      
        pBand = None

        return aData_out, dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value, pGeotransform, pProjection,  pSpatialRef


def gdal_read_geotiff_file_multiple_band(sFilename_in):
    """Read a Geotiff file with multiple bands."""

   
    sDriverName='GTiff'
    pDriver = gdal.GetDriverByName(sDriverName)  

    if pDriver is None:
        print ("%s pDriver not available.\n" % sDriverName)
    else:
        print  ("%s pDriver IS available.\n" % sDriverName) 
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
        pPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]

        pBand = pDataset.GetRasterBand(1)
        
        dummy_data = pBand.ReadAsArray(0, 0, ncolumn, nrow)
        dt = type(dummy_data)
        #there is a chance that GDAL datetype is not compatiable with numpy datatype.
        dMissing_value = pBand.GetNoDataValue()
        aData_out = np.full( (nband, nrow, ncolumn) , -9999.0, dtype= dt )
        for iBand in range(nband):
            pBand = pDataset.GetRasterBand( iBand + 1)
            
            aData_out[iBand, :, :] = pBand.ReadAsArray(0, 0, ncolumn, nrow)
       

      
        pSpatialRef = osr.SpatialReference(wkt=pProjection)
        

        pDriver = None
        pDataset = None
        pBand = None
        return aData_out, pPixelWidth, dOriginX, dOriginY, nband, nrow, ncolumn, dMissing_value, pGeotransform, pProjection,  pSpatialRef