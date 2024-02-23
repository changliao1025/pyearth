import os, sys
import numpy as np
from osgeo import gdal, osr, ogr, gdalconst

def merge_features(sFilename_in, sFilename_out, sFormat):


    pDataset_in = ogr.Open(sFilename_in)    
    # Get the first layer in the file
    pLayer_in = pDataset_in.GetLayer(0)   
    # Count the number of features (polygons)
    nFeature = pLayer_in.GetFeatureCount()
    # Get the spatial reference of the layer
    pSpatial_reference = pLayer_in.GetSpatialRef()
    wkt2 = pSpatial_reference.ExportToWkt()

    if nFeature == 0:
        print('No feature found in the input shapefile')
        return
    else:
        print('Number of features found in the input shapefile: ' + str(nFeature))

    # Create a new dataset using the output filename
    pDriver = ogr.GetDriverByName(sFormat)
    if pDriver is None:
        print('Driver not found')
        return

    if os.path.exists(sFilename_out):
        pDriver.DeleteDataSource(sFilename_out)

    pDataset_out = pDriver.CreateDataSource(sFilename_out)
    if pDataset_out is None:
        print('Dataset not created')
        return    

    #obtain the geotype 
    iGeomType = pLayer_in.GetGeomType()
    if iGeomType == ogr.wkbPoint:
        #create the layer
        pLayer_out = pDataset_out.CreateLayer('layer', pSpatial_reference, geom_type=ogr.wkbPoint)
        pGeometry_merge = ogr.Geometry(ogr.wkbPoint)
    else:
        if iGeomType == ogr.wkbLineString:
            #create the layer
            pLayer_out = pDataset_out.CreateLayer('layer', pSpatial_reference, geom_type=ogr.wkbLineString)
            pGeometry_merge = ogr.Geometry(ogr.wkbLineString)
        else:
            if iGeomType == ogr.wkbPolygon:
                #create the layer
                pLayer_out = pDataset_out.CreateLayer('layer', pSpatial_reference, geom_type=ogr.wkbPolygon)
                pGeometry_merge = ogr.Geometry(ogr.wkbPolygon)
            else:
                print('Geometry type not supported')
                return
            pass
        pass
    # Create a new layer in the output shapefile
    
    # Loop through the input features and merge them into the output layer
    for i in range(nFeature):
        pFeature = pLayer_in.GetFeature(i)
       
        pGeometry = pFeature.GetGeometryRef()                
        if pGeometry is not None:
          
            # Union the geometry of each feature with the merged polygon
            pGeometry_merge = pGeometry_merge.Union(pGeometry)

    # Create a new feature in the output layer
    pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
    pFeature_out.SetGeometry(pGeometry_merge)
    pLayer_out.CreateFeature(pFeature_out)

    #close the dataset
    pDataset_out.Destroy()      
        


    return