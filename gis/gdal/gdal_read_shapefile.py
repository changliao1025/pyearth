import os
import sys
from osgeo import gdal, gdalconst, osr, ogr
def gdal_read_shapefile(sFilename_in):
    driverName = "ESRI Shapefile"
    driver = ogr.GetDriverByName( driverName )
    if driver is None:
        print ("%s driver not available.\n" % driverName)
    else:
        print  ("%s driver IS available.\n" % driverName)
    dataSource = driver.Open(sFilename_in, 0) # 0 means      read-only. 1 means writeable.
    daShapefile = sFilename_in
# Check to see if shapefile is found.
    if dataSource is None:
        print ('Could not open %s' % (daShapefile))
    else:
        print ('Opened %s' % (daShapefile))
        layer = dataSource.GetLayer()
        featureCount = layer.GetFeatureCount()
        print ("Number of features in %s: %d" %  (os.path.basename(daShapefile),featureCount))

        layer = dataSource.GetLayer()
        sr = layer.GetSpatialRef() 

        
        return sr
        #for feature in layer:
            

            #geom = feature.GetGeometryRef()
            #print (geom.Centroid().ExportToWkt())

        #return layer