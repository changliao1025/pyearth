import os, sys
import numpy as np
from osgeo import gdal, osr, ogr, gdalconst

def merge_files(aFilename_in, sFilename_out, sFormat='GeoJSON'):

    # Register all drivers
    ogr.RegisterAll()

    # Get the driver
    pDriver = ogr.GetDriverByName('GeoJSON')
    if pDriver is None:
        print('GeoJSON driver not available.')
        return

    # Check if the output file exists and delete it if it does
    if os.path.exists(sFilename_out):
        pDriver.DeleteDataSource(sFilename_out)

    # Create the output data source
    pDataset_out = pDriver.CreateDataSource(sFilename_out)
    if pDataset_out is None:
        print('Dataset not created')
        return

    # Loop through each input file
    iFlag_first = 1
    lid = 0
    for sFilename_in in aFilename_in:
        # Open the input file
        pDataset_in = ogr.Open(sFilename_in)
        if pDataset_in is None:
            print(f'Failed to open file: {sFilename_in}')
            continue

        # Get the input layer
        pLayer_in = pDataset_in.GetLayer()
        if pLayer_in is None:
            print(f'Failed to get layer from file: {sFilename_in}')
            continue

        # Obtain the geometry type
        iGeomType = pLayer_in.GetGeomType()
        print(f'Geometry type for {sFilename_in}: {iGeomType}')

        # Create the output layer based on the geometry type
        if iFlag_first == 1:
            if iGeomType == ogr.wkbPoint:
                pLayer_out = pDataset_out.CreateLayer('layer', geom_type=ogr.wkbPoint)
            elif iGeomType == ogr.wkbLineString:
                pLayer_out = pDataset_out.CreateLayer('layer', geom_type=ogr.wkbLineString)
            elif iGeomType == ogr.wkbPolygon:
                pLayer_out = pDataset_out.CreateLayer('layer', geom_type=ogr.wkbPolygon)
            else:
                print(f'Unsupported geometry type: {iGeomType}')
                sGeomType = ogr.GeometryTypeToName(iGeomType)
                print('Geometry type not supported:', sGeomType)
                continue
            pField = ogr.FieldDefn('id', ogr.OFTInteger)
            if pLayer_out.CreateField(pField) != 0:
                print('Failed to create id field in output layer')
                return
            iFlag_first = 0
        else:
            pass

        # Copy features from the input layer to the output layer
        for feature in pLayer_in:
            #copy the feature
            pFeature_new = ogr.Feature(pLayer_out.GetLayerDefn())
            # Set the geometry
            pFeature_new.SetGeometry(feature.GetGeometryRef())
            #set id for each feature
            pFeature_new.SetField('id', lid)
            pLayer_out.CreateFeature(pFeature_new)
            lid = lid + 1

        # Clean up
        pDataset_in = None

    # Clean up
    pDataset_out = None
    print('Merge completed.')

    return