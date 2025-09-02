import os
from osgeo import ogr, osr
def convert_polygon_to_polyline_file(sFilename_polygon_in, sFilename_polyline_out):
    if os.path.exists(sFilename_polyline_out):
        os.remove(sFilename_polyline_out)

    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_polyline_out)
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)  # WGS84 lat/lon
    pLayerOut = pDataset_out.CreateLayer('polyline', pSpatial_reference_gcs, ogr.wkbLineString)
    # Add one attribute
    pLayerOut.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64))  # long type for high resolution
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)
    # read the polygon file
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    pDataset_polygon = pDriver_geojson.Open(sFilename_polygon_in, 0)
    # find the number of geometries in the polygon file
    pLayer_polygon = pDataset_polygon.GetLayer()
    nFeature_polygon = pLayer_polygon.GetFeatureCount()
    # check whether the polygon file has the attribute ID
    for pFeature_polygon in pLayer_polygon:
        fid = pFeature_polygon.GetFID()
        pGeometry_polygon = pFeature_polygon.GetGeometryRef()
        if pGeometry_polygon.GetGeometryType() == ogr.wkbPolygon:
            # get the exterior ring
            pRing = pGeometry_polygon.GetGeometryRef(0)
            # create a new line string geometry
            pLineString = ogr.Geometry(ogr.wkbLineString)
            # add the points from the polygon to the line string
            for i in range(pRing.GetPointCount()):
                x, y, z = pRing.GetPoint(i)
                pLineString.AddPoint(x, y, z)
            # create a new feature for the line string
            pFeatureOut.SetGeometry(pLineString)
            # set the ID field
            pFeatureOut.SetField("id", fid)
            # add the feature to the output layer
            pLayerOut.CreateFeature(pFeatureOut)
        else:
            pass
    # clean up
    pFeatureOut.Destroy()
    pDataset_polygon = None
    pDataset_out = None
    print("Conversion completed successfully.")

