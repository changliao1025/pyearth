import os, sys
import numpy as np
from osgeo import gdal, osr

def gdal_read_geotiff_file(sFilename_in, iFlag_metadata_only = 0):
    """Read a Geotiff format raster file.

    Args:
        sFilename_in (string): The file name

    Returns:
        tuple: aData_out, dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value , pGeotransform, pProjection,  pSpatial_reference
    """
    
    if os.path.exists(sFilename_in):
        pass
    else:
        print('The file does not exist!')
        return
    
    

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
        dMissing_value = -9999.0
        pProjection = pDataset.GetProjection()

        pDataset.GetMetadata()
       
        ncolumn = pDataset.RasterXSize
        nrow = pDataset.RasterYSize
        nband = pDataset.RasterCount

        pGeotransform = pDataset.GetGeoTransform()
        dOriginX = pGeotransform[0]
        dOriginY = pGeotransform[3]
        dPixelWidth = pGeotransform[1]
        dPixelHeight = pGeotransform[5]
        #we will use one of them to keep the consistency
        pSpatial_reference = osr.SpatialReference(wkt=pProjection)  

        if iFlag_metadata_only ==1:
            pass
        else:    
            if nband >1: 
                print('This is a multi-band raster file!')
                print('We will only read the first band!')
                pass

            pBand = pDataset.GetRasterBand(1)

            # Data type of the values
            
            iData_type = pBand.DataType
            sDate_type = gdal.GetDataTypeName(iData_type)
            # Compute statistics if needed
            if pBand.GetMinimum() is None or pBand.GetMaximum() is None:
                pBand.ComputeStatistics(0)

            dMissing_value = pBand.GetNoDataValue()
            aData_out = pBand.ReadAsArray(0, 0, ncolumn, nrow)                

        pDataset = None
        pBand = None             

        if iFlag_metadata_only == 1:
            pData = {
                'pixelWidth': dPixelWidth,
                'pixelHeight': dPixelHeight,
                'originX': dOriginX,
                'originY': dOriginY,
                'nrow': nrow,
                'ncolumn': ncolumn,
                #'missingValue': dMissing_value,
                'geotransform': pGeotransform,
                'projection': pProjection,
                'spatialReference': pSpatial_reference
            }
            return pData
           
        else:   
            pData = {
                'dataOut': aData_out,
                'dataType': iData_type, # 'Float32', 'Float64', 'UInt16', 'Int16', 'UInt32', 'Int32', 'Byte', 'UInt32', 'Int32', 'CInt16', 'CInt32', 'CFloat32', 'CFloat64
                'pixelWidth': dPixelWidth,
                'pixelHeight': dPixelHeight,
                'originX': dOriginX,
                'originY': dOriginY,
                'nrow': nrow,
                'ncolumn': ncolumn,
                'missingValue': dMissing_value,
                'geotransform': pGeotransform,
                'projection': pProjection,
                'spatialReference': pSpatial_reference
            }
            return pData

def gdal_read_geotiff_file_multiple_band(sFilename_in, iFlag_metadata_only = 0):
    
    """Read a Geotiff format raster file with multiple bands.

    Args:
        sFilename_in (string): The filename

    Returns:
        tuple: aData_out, dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value , pGeotransform, pProjection,  pSpatial_reference
    """
    
    if os.path.exists(sFilename_in):
        pass
    else:
        print('The file does not exist!')
        return
   
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
        dPixelHeight = pGeotransform[5]

        pBand = pDataset.GetRasterBand(1)
        
        dummy_data = pBand.ReadAsArray(0, 0, ncolumn, nrow)
        dt = type(dummy_data)
        #there is a chance that GDAL datetype is not compatiable with numpy datatype.
        dMissing_value = pBand.GetNoDataValue()
        pSpatial_reference = osr.SpatialReference(wkt=pProjection)

        if iFlag_metadata_only ==1:
            pass
        else:    
            aData_out = np.full( (nband, nrow, ncolumn) , -9999.0, dtype= dt )
            for iBand in range(nband):
                pBand = pDataset.GetRasterBand( iBand + 1)
                iData_type = pBand.DataType
                aData_out[iBand, :, :] = pBand.ReadAsArray(0, 0, ncolumn, nrow)
                pass

        pDriver = None
        pDataset = None
        pBand = None
        if iFlag_metadata_only == 1:
            pData = {
                'pixelWidth': dPixelWidth,
                'pixelHeight': dPixelHeight,
                'originX': dOriginX,
                'originY': dOriginY,
                'nrow': nrow,
                'ncolumn': ncolumn,
                #'missingValue': dMissing_value,
                'geotransform': pGeotransform,
                'projection': pProjection,
                'spatialReference': pSpatial_reference
            }
            return pData
        else:
            pData = {
                'dataOut': aData_out,
                'dataType': iData_type,
                'pixelWidth': dPixelWidth,
                'pixelHeight': dPixelHeight,
                'originX': dOriginX,
                'originY': dOriginY,
                'nrow': nrow,
                'ncolumn': ncolumn,
                'missingValue': dMissing_value,
                'geotransform': pGeotransform,
                'projection': pProjection,
                'spatialReference': pSpatial_reference
            }
            return pData