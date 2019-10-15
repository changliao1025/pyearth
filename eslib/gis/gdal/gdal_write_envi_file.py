import sys
from osgeo import gdal,  gdalconst

def gdal_write_envi_file(sFilename_out, aData_in, dPixelWidth,\
     dOriginX, dOriginY, pProjection ):
    
    pDriverName = gdal.GetDriverByName('ENVI')
    pDriverName.Register()
	nrow, ncolumn = aData_in.shape
	nband = 1

	# Creates a new raster data source
	pFile = driver.Create(sFilename_out, ncolumn, nrow, nband, gdal.GDT_Float32)
	
	# Write metadata	

	pFile.SetGeoTransform([dOriginX, dPixelWidth, 0.0, dOriginY, 0.0, -dPixelWidth])
	pFile.SetProjection(pProjection)

	#Write raster datasets
	outBand = pFile.GetRasterBand(1)
	outBand.WriteArray(aData_in)
	
	pGeotransform_out = pFile.GetGeoTransform()
	pProjection_out = pFile.GetProjection()
	
	print ("Output binary saved: ", sFilename_out)
	
	return  sFilename_out, pGeotransform_out, pProjection_out