
import os, sys
import numpy as np
from osgeo import ogr, gdal, osr
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area  import calculate_polygon_area

from gcsbuffer.classes.vertex import pyvertex
from gcsbuffer.classes.edge import pyedge
from gcsbuffer.classes.polyline import pypolyline
from gcsbuffer.classes.polygon import pypolygon

pDriver = ogr.GetDriverByName('GeoJSON')

def create_point_buffer_zone(sWkt, dBuffer_distance_in):
    #create a point geometry from WKT
    pGeometry = ogr.CreateGeometryFromWkt(sWkt)
    if pGeometry is None:
        print('Error: Invalid WKT input for point buffer creation.')
        return None
    sGeometry_type = pGeometry.GetGeometryName()
    if sGeometry_type != 'POINT':
        print('Error: Input geometry must be a POINT for buffer creation.')
        return None

    aCoords_gcs = get_geometry_coordinates(pGeometry)
    if len(aCoords_gcs) != 1:
        print('Error: Input geometry must be a single point for buffer creation.')
        return None

    point = dict()
    point['dLongitude_degree'] = aCoords_gcs[0][0]
    point['dLatitude_degree'] = aCoords_gcs[0][1]
    pVertex = pyvertex(point)

    sWkt_buffer_polygon = pVertex.calculate_buffer_zone(dBuffer_distance_in)

    return sWkt_buffer_polygon

def create_polyline_buffer_zone(sWkt, dBuffer_distance_in):

    #create a polyline geometry from WKT
    pGeometry = ogr.CreateGeometryFromWkt(sWkt)
    if pGeometry is None:
        print('Error: Invalid WKT input for polyline buffer creation.')
        return None
    sGeometry_type = pGeometry.GetGeometryName()
    if sGeometry_type != 'LINESTRING':
        print('Error: Input geometry must be a LINESTRING for buffer creation.')
        return None

    aEdge = list()
    aCoords_gcs = get_geometry_coordinates(pGeometry)
    nPoint = len(aCoords_gcs) #remove the last point which is the same as the first point
    aVertex = list()
    point= dict()
    for i in range(0, nPoint):
        point['dLongitude_degree'] = aCoords_gcs[i][0]
        point['dLatitude_degree'] =  aCoords_gcs[i][1]
        pVertex = pyvertex(point)
        aVertex.append(pVertex)

    nVertex = len(aVertex)
    for i in range(0, nVertex-1):
        if aVertex[i] != aVertex[i+1]:
            pEdge = pyedge(aVertex[i], aVertex[i+1])
            aEdge.append(pEdge)
        else:
            pass

    ppolyline = pypolyline(aEdge)
    sWkt_buffer_polygon = ppolyline.calculate_buffer_zone(dBuffer_distance_in)


    return sWkt_buffer_polygon

def create_buffer_zone_polygon_file(sFilename_polygon_in, sFilename_polygon_out,
                                    dThreshold_in = 1.0E9, #m2 to filter out small polygons
                                      dBuffer_distance_in = 5000): #m

    pDataSource = pDriver.Open(sFilename_polygon_in, 0)
    pLayer = pDataSource.GetLayer()

    # Prepare output (overwrite if exists)
    if os.path.exists(sFilename_polygon_out):
        pDriver.DeleteDataSource(sFilename_polygon_out)
    pOutDataSource = pDriver.CreateDataSource(sFilename_polygon_out)
    pOutLayer = pOutDataSource.CreateLayer("buffer", geom_type=ogr.wkbPolygon)

    pFeature = pLayer.GetNextFeature()

    while pFeature:
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
            if dThreshold_in is None:
                pass
            else:
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
            sWkt_buffer = pPolygon.calculate_buffer_zone(dBuffer_distance_in)

            if sWkt_buffer:
                buffer_geom = ogr.CreateGeometryFromWkt(sWkt_buffer)
                outFeature = ogr.Feature(pOutLayer.GetLayerDefn())
                outFeature.SetGeometry(buffer_geom)
                pOutLayer.CreateFeature(outFeature)
                outFeature = None  # Free feature

        pFeature = pLayer.GetNextFeature()

    pDataSource = None
    pOutDataSource = None
    print('create_gcs_buffer_zone_polygon is done!')
    return



