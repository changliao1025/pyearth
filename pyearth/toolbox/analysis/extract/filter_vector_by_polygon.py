import os
from osgeo import osr, ogr

def filter_vector_by_polygon(sFilename_vector_in, sFilename_polygon_in, sFilename_vector_out):
    # Open the input vector layer
    pDriver_parquet = ogr.GetDriverByName('Parquet')
    pDataSource_in = ogr.Open(sFilename_vector_in)
    if pDataSource_in is None:
        print("Could not open input vector")
        return

    pLayer_in = pDataSource_in.GetLayer()

    # Open the clip polygon layer
    pDataSource_clip = ogr.Open(sFilename_polygon_in)
    if pDataSource_clip is None:
        print("Could not open clip polygon")
        pDataSource_in = None
        return

    pLayer_clip = pDataSource_clip.GetLayer()

    # Create a new shapefile for the clipped data
    
    pDataSource_out = pDriver_parquet.CreateDataSource(sFilename_vector_out)
    pLayer_out = pDataSource_out.CreateLayer('clipped', geom_type=ogr.wkbPolygon)

    # Copy the fields from the input layer to the output layer
    pLayer_defn_in = pLayer_in.GetLayerDefn()
    for i in range(pLayer_defn_in.GetFieldCount()):
        field_defn = pLayer_defn_in.GetFieldDefn(i)
        pLayer_out.CreateField(field_defn)   
    
    iFlag_option = 1
    if iFlag_option == 1:
        # Get the extent (bounding box) of the clip layer
        aExtent_clip = pLayer_clip.GetExtent()
        print(aExtent_clip)
        pLayer_in.SetSpatialFilterRect(aExtent_clip[0], aExtent_clip[2], aExtent_clip[1], aExtent_clip[3])
        pass
    else:
        if iFlag_option == 2: #use the outline of the clip layer, this is more accurate
            pass

    # Loop through the input layer and create clipped features
    for feature in pLayer_in:
        pGeometry = feature.GetGeometryRef()        
        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        pFeature_out.SetGeometry(pGeometry)
        for i in range(pLayer_out.GetLayerDefn().GetFieldCount()):
            pFeature_out.SetField(pLayer_out.GetLayerDefn().GetFieldDefn(i).GetNameRef(),
                                feature.GetField(i))
        pLayer_out.CreateFeature(pFeature_out)
        pFeature_out = None

    # Reset spatial filter and close data sources
    pLayer_in.SetSpatialFilter(None)
    pDataSource_in = pDataSource_clip = pDataSource_out = None