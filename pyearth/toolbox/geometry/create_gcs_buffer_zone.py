
import os, sys
import numpy as np
from osgeo import ogr, gdal, osr
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area  import calculate_polygon_area

from gcsbuffer.classes.vertex import pyvertex
from gcsbuffer.classes.edge import pyedge
from gcsbuffer.classes.polygon import pypolygon

pDriver = ogr.GetDriverByName('GeoJSON')

def create_gcs_buffer_zone_polygon(sFilename_polygon_in, sFilename_polygon_out,
                                    dThreshold_in = 1.0E9, #m2 to filter out small polygons
                                      dBuffer_distance_in = 5000): #m

    pDataSource = pDriver.Open(sFilename_polygon_in, 0)
    pLayer = pDataSource.GetLayer()
    pFeature = pLayer.GetNextFeature()
    pGeometry = pFeature.GetGeometryRef()

    sGeometry_type = pGeometry.GetGeometryName()
    if sGeometry_type =='MULTIPOLYGON':
        aaCoords_gcs = get_geometry_coordinates(pGeometry)
        nPart = len(aaCoords_gcs)
    else:
        aCoords_gcs = get_geometry_coordinates(pGeometry)
        aaCoords_gcs= [aCoords_gcs]
        nPart = 1
    for i in range(nPart):
        aCoords_gcs = aaCoords_gcs[i]
        aCoords_gcs = np.array(aCoords_gcs)
        #calculate the area of the polygon
        dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1], iFlag_algorithm=2)
        if (dArea) < dThreshold_in:
            continue


        nPoint = len(aCoords_gcs)-1 #remove the last point which is the same as the first point
        aVertex = list()
        point= dict()
        for i in range(0, nPoint):
            point['dLongitude_degree'] = aCoords_gcs[i,0]
            point['dLatitude_degree'] =  aCoords_gcs[i,1]
            pVertex = pyvertex(point)
            aVertex.append(pVertex)

    nVertex = len(aVertex)
    aEdge = list()
    for i in range(0, nVertex-1):
        if aVertex[i] != aVertex[i+1]:
            pEdge = pyedge(aVertex[i], aVertex[i+1])
            aEdge.append(pEdge)
        else:
            pass

    pEdge = pyedge(aVertex[nVertex-1], aVertex[0])
    aEdge.append(pEdge)
    pPolygon = pypolygon(aEdge)
    pPolygon.calculate_buffer_zone_polygon(dBuffer_distance_in, sFilename_out = sFilename_polygon_out)
    print('create_gcs_buffer_zone_polygon is done!')
    return