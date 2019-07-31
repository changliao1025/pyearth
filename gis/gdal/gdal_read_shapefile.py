import os
import sys
from osgeo import gdalconst, ogr, osr
def gdal_read_shapefile(sFilename_in):
    sDriverName = "ESRI Shapefile"
    pDriver = ogr.GetDriverByName( sDriverName )
    if pDriver is None:
        print ("%s pDriver not available.\n" % sDriverName)
    else:
        print  ("%s pDriver IS available.\n" % sDriverName)
    pDataSource = pDriver.Open(sFilename_in, 0) # 0 means      read-only. 1 means writeable.
    
    # Check to see if shapefile is found.
    if pDataSource is None:
        print ('Could not open %s' % (sFilename_in))
    else:
        print ('Opened %s' % (sFilename_in))
        pLayer = pDataSource.GetLayer()
        lFeatureCount = pLayer.GetFeatureCount()
        print ("Number of features in %s: %d" %  (os.path.basename(sFilename_in),lFeatureCount))

        
        pSpatailRef = pLayer.GetSpatialRef() 

        
        aFeature=[]
        for pFeature in pLayer:         

            pGeometry = pFeature.GetGeometryRef()
            print (pGeometry.Centroid().ExportToWkt())
            aFeature.append(pGeometry)

        return aFeature, pSpatailRef