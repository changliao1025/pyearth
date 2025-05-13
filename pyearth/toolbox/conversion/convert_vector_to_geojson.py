import os
import osgeo
from osgeo import ogr, osr, gdal
from pyearth.system.define_global_variables import *
def convert_vector_to_geojson(sFilename_vector_in, sFilename_geojson_out):
    #check input file
    if not os.path.exists(sFilename_vector_in):
        print('Input file not found:' , sFilename_vector_in)
        return
    #check if the file exists
    if os.path.exists(sFilename_geojson_out):
        os.remove(sFilename_geojson_out)
        pass

    #open the vector file
    pDataset_in = ogr.Open(sFilename_vector_in, 0)
    if pDataset_in is None:
        print('Could not open file')
        return
    #create the driver
    pDriver = ogr.GetDriverByName('GeoJSON')
    if pDriver is None:
        print('Driver not available')
        return
    #create the output file
    pDataset_out = pDriver.CreateDataSource(sFilename_geojson_out)
    if pDataset_out is None:
        print('Could not create file')
        return


    pSrs_out = osr.SpatialReference()
    pSrs_out.ImportFromEPSG(4326)    # WGS84 lat/lon
    wkt2 = pSrs_out.ExportToWkt()
    print(wkt2)

    #get the layer
    pLayer_in = pDataset_in.GetLayer()
    #obtain the geotype
    #get the first feature
    pLayer_in.ResetReading()
    pFeature = pLayer_in.GetNextFeature()
    #get the geometry
    pGeometry = pFeature.GetGeometryRef()
    #get geometry type
    iGeomType = pGeometry.GetGeometryType()
    #get the spatial reference
    pSrs_in = pLayer_in.GetSpatialRef()
    wkt1 = pSrs_in.ExportToWkt()
    print(wkt1)

    comparison = pSrs_in.IsSame(pSrs_out)
    #if(comparison != 1):
    if( wkt1 != wkt2):
        iFlag_transform = 1
        if int(osgeo.__version__[0]) >= 3:
            # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
            pSrs_in.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
            pSrs_out.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

        transform = osr.CoordinateTransformation(pSrs_in, pSrs_out)
        pass
    else:
        iFlag_transform = 0
        pass

    if iGeomType == ogr.wkbPoint:
        #create the layer
        pLayer_out = pDataset_out.CreateLayer('layer', pSrs_out, geom_type=ogr.wkbPoint)
    else:
        if iGeomType == ogr.wkbLineString:
            #create the layer
            pLayer_out = pDataset_out.CreateLayer('layer', pSrs_out, geom_type=ogr.wkbLineString)
        else:
            if iGeomType == ogr.wkbPolygon:
                #create the layer
                pLayer_out = pDataset_out.CreateLayer('layer', pSrs_out, geom_type=ogr.wkbPolygon)
            else:
                print('Geometry type not supported')
                return
            pass
        pass

    pFeatureDefn = pLayer_in.GetLayerDefn()
    pFeature_out = ogr.Feature(pFeatureDefn)
    #copy the fields
    for i in range(pFeatureDefn.GetFieldCount()):
        pFieldDefn = pFeatureDefn.GetFieldDefn(i)
        pLayer_out.CreateField(pFieldDefn)
        pass

    if iFlag_transform == 1:
        #get the feature count
        iFeatureCount = pLayer_in.GetFeatureCount()
        pLayer_in.ResetReading()
        #for i in range(iFeatureCount):
        pFeature_in = pLayer_in.GetNextFeature()
        while pFeature_in:
            #pFeature_in = pLayer_in.GetFeature(i)
            pGeometry = pFeature_in.GetGeometryRef()
            pGeometry.Transform(transform)
            pFeature_out = ogr.Feature(pFeatureDefn)
            pFeature_out.SetGeometry(pGeometry)
            #copy the attributes
            for i in range(pFeatureDefn.GetFieldCount()):
                pFeature_out.SetField(i, pFeature_in.GetField(i))
                pass

            pLayer_out.CreateFeature(pFeature_out)
            pFeature_in = pLayer_in.GetNextFeature()
        pass
    else:
        iFeatureCount = pLayer_in.GetFeatureCount()
        pLayer_in.ResetReading()
        #for i in range(iFeatureCount):
        pFeature_in = pLayer_in.GetNextFeature()
        while pFeature_in:
            #pFeature_in = pLayer_in.GetFeature(i)
            pGeometry = pFeature_in.GetGeometryRef()
            pFeature_out.SetGeometry(pGeometry)
            #copy the attributes
            for i in range(pFeatureDefn.GetFieldCount()):
                pFeature_out.SetField(i, pFeature_in.GetField(i))
            pLayer_out.CreateFeature(pFeature_out)
            pFeature_in = pLayer_in.GetNextFeature()
        pass

    pDataset_in = None
    pLayer_out = pFeature_out = None
    pDataset_out = None

    # Verify the output geometry type
    #dataset_out = ogr.Open(sFilename_geojson_out)
    #if dataset_out is None:
    #    print(f'Failed to open output file: {sFilename_geojson_out}')
    #    return
#
    #layer_out = dataset_out.GetLayer()
    #if layer_out is None:
    #    print(f'Failed to get layer from output file: {sFilename_geojson_out}')
    #    return
#
    #iGeomType0 = layer_out.GetGeomType()
    #sGeomType0 = ogr.GeometryTypeToName(iGeomType0)
    #for feature in layer_out:
    #    geometry = feature.GetGeometryRef()
    #    geom_type_out = geometry.GetGeometryType()
    #    sGeomType1 = geometry.GetGeometryName() #ogr.GeometryTypeToName(geom_type_out)
    #    print(f'Output geometry type: {sGeomType1}')

    print('Conversion completed')


    return
