import os, sys
import numpy as np
from osgeo import ogr, osr, gdal
gdal.UseExceptions()
from pyearth.system.define_global_variables import *
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.douglas_peucker_geodetic import douglas_peucker_geodetic
from pyearth.gis.geometry.visvalingam_whyatt_geodetic import  visvalingam_whyatt_geodetic

def remove_small_polygon(sFilename_vector_in, sFilename_vector_out, dThreshold_in):
    """
    This function is used to remove small polygons from a vector file
    :param sFilename_vector_in: input vector file
    :param sFilename_vector_out: output vector file
    :param dThreshold_in: threshold for removing small polygons, in square kilo meters
    :param dSimplify_tolerance: tolerance for simplifying the geometry, it depends on the algorithm used
    :return: None
    """

    if not os.path.exists(sFilename_vector_in):
        print('Error: file %s does not exist!' % sFilename_vector_in)
        return

    if os.path.exists(sFilename_vector_out):
        os.remove(sFilename_vector_out)

    dThreshold = float(dThreshold_in)
    pDriver = ogr.GetDriverByName('GeoJSON')
    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)

    pDataSource_in = pDriver.Open(sFilename_vector_in, 0)
    pLayer_in = pDataSource_in.GetLayer()

    pDataSource_out = pDriver.CreateDataSource(sFilename_vector_out)
    pLayer_out = pDataSource_out.CreateLayer('layer', pSrs, ogr.wkbPolygon)
    #create fields for id and area
    pFieldDefn = ogr.FieldDefn('id', ogr.OFTInteger)
    pLayer_out.CreateField(pFieldDefn)
    pFieldDefn = ogr.FieldDefn('area', ogr.OFTReal)
    pLayer_out.CreateField(pFieldDefn)
    pLayerDefn_in = pLayer_in.GetLayerDefn()
    nFieldCount = pLayerDefn_in.GetFieldCount()
    for i in range(nFieldCount):
        pFieldDefn = pLayerDefn_in.GetFieldDefn(i)
        if pFieldDefn.GetName() == 'id':
            continue
        if pFieldDefn.GetName() == 'area':
            continue
        pLayer_out.CreateField(pFieldDefn)

    pLayerDefn_out = pLayer_out.GetLayerDefn()
    pFeature_in = pLayer_in.GetNextFeature()
    lID = 1
    while pFeature_in:
        pGeometry = pFeature_in.GetGeometryRef()
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
            dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1], iFlag_algorithm=2)
            dAreakm = dArea/1.0E6

            if dAreakm > dThreshold:
                pGeometry_partial = ogr.Geometry(ogr.wkbPolygon)
                pRing = ogr.Geometry(ogr.wkbLinearRing)
                for aCoord in aCoords_gcs:
                    pRing.AddPoint(aCoord[0], aCoord[1])
                pRing.CloseRings()
                pGeometry_partial.AddGeometry(pRing)
                pGeometry_partial.AssignSpatialReference(pSrs)
                pFeature_out = ogr.Feature(pLayerDefn_out)
                pFeature_out.SetGeometry(pGeometry_partial)
                pFeature_out.SetField('id', lID)
                pFeature_out.SetField('area', dAreakm)
                for i in range(nFieldCount):
                    field_name = pLayerDefn_in.GetFieldDefn(i).GetName()
                    if field_name == 'id' or field_name == 'area':
                        continue
                    pFeature_out.SetField(field_name, pFeature_in.GetField(field_name))

                pLayer_out.CreateFeature(pFeature_out)
                pFeature_out = None
                lID = lID + 1

        pFeature_in = pLayer_in.GetNextFeature()

    pDataSource_in.Destroy()
    pDataSource_out.Destroy()
    print('The small polygons have been removed!')
    return