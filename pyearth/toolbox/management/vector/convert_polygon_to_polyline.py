import os, sys
from osgeo import ogr, osr

def convert_polygon_to_polyline(sFilename_polygon_in, sFilename_polyline_out, sFormat='GeoJSON'):
    # Open the input polygon file
    if sFormat == 'GeoJSON':
        pDriver_polygon = ogr.GetDriverByName('GeoJSON')
    else:
        pDriver_polygon = ogr.GetDriverByName('ESRI Shapefile')

    pDriver_polyline = pDriver_polygon

    if not os.path.exists(sFilename_polygon_in):
        print(f"Error: file {sFilename_polygon_in} does not exist!")
        return

    pDataset_polygon = pDriver_polygon.Open(sFilename_polygon_in, 0)  # 0 means read-only

    if pDataset_polygon is None:
        raise FileNotFoundError(f"Could not open {sFilename_polygon_in}")

    #check whether the output is already exists
    if os.path.exists(sFilename_polyline_out):
        pDriver_polyline.DeleteDataSource(sFilename_polyline_out)

    # Get the polygon layer
    pLayer_polygon = pDataset_polygon.GetLayer()

    # Create the output polyline file
    pDataset_polyline = pDriver_polyline.CreateDataSource(sFilename_polyline_out)
    if pDataset_polyline is None:
        raise RuntimeError(f"Could not create {sFilename_polyline_out}")

    # Create the spatial reference for the output file
    pSpatialRef = pLayer_polygon.GetSpatialRef()

    # Create the polyline layer
    pLayer_polyline = pDataset_polyline.CreateLayer('polyline', pSpatialRef, ogr.wkbLineString)

    # Copy the fields from the polygon layer to the polyline layer
    pLayerDefn_polygon = pLayer_polygon.GetLayerDefn()
    for i in range(pLayerDefn_polygon.GetFieldCount()):
        fieldDefn = pLayerDefn_polygon.GetFieldDefn(i)
        pLayer_polyline.CreateField(fieldDefn)

    # Get the feature definition for the polyline layer
    pLayerDefn_polyline = pLayer_polyline.GetLayerDefn()

    # Iterate over the features in the polygon layer
    for pFeature_polygon in pLayer_polygon:
        # Get the geometry of the polygon feature
        pGeometry_polygon = pFeature_polygon.GetGeometryRef()
        if pGeometry_polygon.GetGeometryName() == 'POLYGON':
            pGeometry_polyline = ogr.Geometry(ogr.wkbLineString)
            ring = pGeometry_polygon.GetGeometryRef(0)  # Get the exterior ring
            npoints = ring.GetPointCount()
            for i in range(npoints):
                point = ring.GetPoint(i)
                pGeometry_polyline.AddPoint(point[0], point[1])

            pFeature_polyline = ogr.Feature(pLayerDefn_polyline)
            pFeature_polyline.SetGeometry(pGeometry_polyline)
            pLayer_polyline.CreateFeature(pFeature_polyline)
            pFeature_polyline = None
        else:
            if pGeometry_polygon.GetGeometryName() == 'MULTIPOLYGON':
                for i in range(pGeometry_polygon.GetGeometryCount()):
                    pGeometry_polyline = ogr.Geometry(ogr.wkbLineString)
                    polygon = pGeometry_polygon.GetGeometryRef(i)
                    ring = polygon.GetGeometryRef(0)  # Get the exterior ring of each polygon
                    npoints = ring.GetPointCount()
                    for j in range(npoints):
                        point = ring.GetPoint(j)
                        pGeometry_polyline.AddPoint(point[0], point[1])

                    # Create a new feature for the polyline layer
                    pFeature_polyline = ogr.Feature(pLayerDefn_polyline)
                    pFeature_polyline.SetGeometry(pGeometry_polyline)
                    # Add the polyline feature to the polyline layer
                    pLayer_polyline.CreateFeature(pFeature_polyline)
                    pFeature_polyline = None

    # Close the datasets
    pDataset_polygon = None
    pDataset_polyline = None
    pSpatialRef = None

    print(f"Conversion from polygon to polyline completed. Output saved to {sFilename_polyline_out}")

