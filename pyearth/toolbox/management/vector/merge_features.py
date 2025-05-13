import os, sys
import numpy as np
from osgeo import gdal, osr, ogr, gdalconst

def merge_features(sFilename_in, sFilename_out, sFormat='GeoJSON'):

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

    #obtain the geotype of first layer and
    iGeomType = pLayer_in.GetGeomType()
    #obtain the geotype of first geometry
    pLayer_in.ResetReading()

    # Obtain the first feature
    pFeature = pLayer_in.GetNextFeature()
    pGeometry = pFeature.GetGeometryRef()
    if pGeometry is None:
        print('Geometry not found')
        return
    #get the geometry type
    iGeomType = pGeometry.GetGeometryType()
    #get geometry type name
    sGeomType = ogr.GeometryTypeToName(iGeomType)
    print('Geometry type: ' + sGeomType)
    #check whether it is a multi-geometry
    if iGeomType == ogr.wkbMultiPoint or iGeomType == ogr.wkbMultiLineString or iGeomType == ogr.wkbMultiPolygon:
        #get the number of geometries
        nGeom = pGeometry.GetGeometryCount()
        #get the first geometry
        pGeometry_single = pGeometry.GetGeometryRef(0)
        iGeomType = pGeometry_single.GetGeometryType()

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
    pLayer_in.ResetReading()  # Reset reading to start from the first feature again
    pFeature = pLayer_in.GetNextFeature()
    while pFeature:
        pGeometry = pFeature.GetGeometryRef()
        if pGeometry is not None:
            #check geotype again
            iGeomType_new = pGeometry.GetGeometryType()
            if iGeomType_new == iGeomType:
                # Union the geometry of each feature with the merged polygon
                pGeometry_merge = pGeometry_merge.Union(pGeometry)
            else:
                #check whether the geometry type is a multi-geometry
                if iGeomType_new == ogr.wkbMultiPoint or iGeomType_new == ogr.wkbMultiLineString or iGeomType_new == ogr.wkbMultiPolygon:
                    #get the number of geometries
                    nGeom = pGeometry.GetGeometryCount()
                    for i in range(nGeom):
                        pGeometry_single = pGeometry.GetGeometryRef(i)
                        #check again its geometry type
                        iGeomType_single = pGeometry_single.GetGeometryType()
                        if iGeomType_single == iGeomType:
                            pGeometry_merge = pGeometry_merge.Union(pGeometry_single)
                        else:
                            print('Geometry type not supported')
                else:
                    print('Geometry type not supported')

            pFeature = pLayer_in.GetNextFeature()

    # Create a new feature in the output layer
    pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
    pFeature_out.SetGeometry(pGeometry_merge)
    pLayer_out.CreateFeature(pFeature_out)

    #close the dataset
    pDataset_out.Destroy()
    pDataset_in.Destroy()
    print('Finished merging features.')

    return