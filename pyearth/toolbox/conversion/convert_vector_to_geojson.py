import os
from osgeo import ogr, osr, gdal
def convert_vector_to_geojson(sFilename_vector_in, sFilename_geojson_out):
    #check input file
    if not os.path.exists(sFilename_vector_in):
        print('Input file not found')
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
    iGeomType = pLayer_in.GetGeomType()
    #get the spatial reference
    pSrs_in = pLayer_in.GetSpatialRef()
    wkt1 = pSrs_in.ExportToWkt()
    print(wkt1)

    comparison = pSrs_in.IsSame(pSrs_out)    
    #if(comparison != 1):
    if( wkt1 != wkt2):
        iFlag_transform = 1
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
        
    if iFlag_transform == 1:
        
        #get the feature count 
        iFeatureCount = pLayer_in.GetFeatureCount()

        for i in range(iFeatureCount):
            pFeature_in = pLayer_in.GetFeature(i)    
            pGeometry = pFeature_in.GetGeometryRef()
            pGeometry.Transform(transform)
            pFeatureDefn = pLayer_in.GetLayerDefn()
            pFeature_out = ogr.Feature(pFeatureDefn)
            pFeature_out.SetGeometry(pGeometry)
            pLayer_out.CreateFeature(pFeature_out)            
        pass
    else:
        iFeatureCount = pLayer_in.GetFeatureCount()
        for i in range(iFeatureCount):
            pFeature_in = pLayer_in.GetFeature(i)    
            pGeometry = pFeature_in.GetGeometryRef()
            pFeatureDefn = pLayer_in.GetLayerDefn()
            pFeature_out = ogr.Feature(pFeatureDefn)
            pFeature_out.SetGeometry(pGeometry)
            pLayer_out.CreateFeature(pFeature_out)            
        pass

    pDataset_in = None
    pDataset_out = None
    
    print('Conversion completed')


    return
