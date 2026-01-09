# import gdal
import os
from osgeo import ogr, osr
from datetime import datetime
import numpy as np
from rtree.index import Index as RTreeindex
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area

osr.UseExceptions()


# these functions assume that the input vector files have an attribute called 'id'
def polygon_difference(
    sFilename_base,
    sFilename_new,
    sAttribute_name_base,
    sAttribute_name_new,
    sFilename_output_in,
):
    start_time = datetime.now()
    # read the base file
    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    # delete if file exist
    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)
    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_output_in)
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)  # WGS84 lat/lon

    pLayerOut = pDataset_out.CreateLayer("diff", pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayerOut.CreateField(
        ogr.FieldDefn("polygonid", ogr.OFTInteger64)
    )  # long type for high resolution
    pLayerOut.CreateField(
        ogr.FieldDefn("diff", ogr.OFTReal)
    )  # long type for high resolution

    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)

    pDataset_base = pDriver_geojson.Open(sFilename_base, 0)
    # read the new file
    pDataset_new = pDriver_geojson.Open(sFilename_new, 0)

    # find the numebr of geometries in the base file
    pLayer_base = pDataset_base.GetLayer()
    nFeature_base = pLayer_base.GetFeatureCount()

    # do the same for the new file
    pLayer_new = pDataset_new.GetLayer()
    nFeature_new = pLayer_new.GetFeatureCount()
    lID_polygon = 1
    for i in range(nFeature_base):
        pFeature_base = pLayer_base.GetFeature(i)
        pGeometry_base = pFeature_base.GetGeometryRef()
        dRunoff_base = pFeature_base.GetField(sAttribute_name_base)

        for j in range(nFeature_new):
            pFeature_new = pLayer_new.GetFeature(j)
            pGeometry_new = pFeature_new.GetGeometryRef()
            iFlag_intersect = pGeometry_new.Intersects(pGeometry_base)
            if iFlag_intersect == True:
                iFlag_intersected = 1
                pGeometry_intersect = pGeometry_new.Intersection(pGeometry_base)
                pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                if pGeometrytype_intersect == "POLYGON":
                    dRunoff_new = pFeature_new.GetField(sAttribute_name_new)
                    dData_diff = dRunoff_new - dRunoff_base
                    pFeatureOut.SetGeometry(pGeometry_intersect)
                    pFeatureOut.SetField("polygonid", lID_polygon)
                    pFeatureOut.SetField("diff", dData_diff)

                    pLayerOut.CreateFeature(pFeatureOut)
                    lID_polygon = lID_polygon + 1

    # close files
    pDataset_base = None
    pDataset_new = None
    end_time = datetime.now()
    # get difference
    delta = end_time - start_time
    sec = delta.total_seconds()
    print("difference in seconds:", sec)
    min = sec / 60
    print("difference in minutes:", min)
    return None


def polygon_difference_rtree(
    sFilename_base,
    sFilename_new,
    sAttribute_name_base,
    sAttribute_name_new,
    sFilename_output_in,
):
    import rtree

    start_time = datetime.now()
    # read the base file
    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    # delete if file exist
    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)

    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_output_in)
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)  # WGS84 lat/lon
    pLayerOut = pDataset_out.CreateLayer("diff", pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayerOut.CreateField(
        ogr.FieldDefn("polygonid", ogr.OFTInteger64)
    )  # long type for high resolution
    pLayerOut.CreateField(
        ogr.FieldDefn("diff", ogr.OFTReal)
    )  # long type for high resolution
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)
    # read the base file
    pDataset_base = pDriver_geojson.Open(sFilename_base, 0)
    # find the numebr of geometries in the base file
    pLayer_base = pDataset_base.GetLayer()
    interleaved = True
    index_base = RTreeindex()
    aData_base = list()
    lID = 0
    for pFeature_base in pLayer_base:
        pGeometry_base = pFeature_base.GetGeometryRef()
        left, right, bottom, top = pGeometry_base.GetEnvelope()
        pBound = (left, bottom, right, top)
        index_base.insert(lID, pBound)  #
        dData_base = pFeature_base.GetField(sAttribute_name_base)
        aData_base.append(dData_base)
        lID = lID + 1

    # read the new file
    pDataset_new = pDriver_geojson.Open(sFilename_new, 0)
    pLayer_new = pDataset_new.GetLayer()
    lID_polygon = 1
    for pFeature_new in pLayer_new:
        pGeometry_new = pFeature_new.GetGeometryRef()
        left, right, bottom, top = pGeometry_new.GetEnvelope()
        pBound = (left, bottom, right, top)
        aIntersect = list(index_base.intersection(pBound))
        dData_new = pFeature_new.GetField(sAttribute_name_new)
        for k in aIntersect:
            pFeature_base = pLayer_base.GetFeature(k)
            pGeometry_base = pFeature_base.GetGeometryRef()
            dData_base = aData_base[k]
            dData_diff = dData_new - dData_base
            iFlag_intersect = pGeometry_new.Intersects(pGeometry_base)
            if iFlag_intersect == True:
                pGeometry_intersect = pGeometry_new.Intersection(pGeometry_base)
                pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                if pGeometrytype_intersect == "POLYGON":
                    pFeatureOut.SetGeometry(pGeometry_intersect)
                    pFeatureOut.SetField("polygonid", lID_polygon)
                    pFeatureOut.SetField("diff", dData_diff)
                    pLayerOut.CreateFeature(pFeatureOut)
                    lID_polygon = lID_polygon + 1
                else:
                    print("pGeometrytype_intersect:", pGeometrytype_intersect)

    # close files
    pDataset_base = None
    pDataset_new = None
    end_time = datetime.now()
    # get difference
    delta = end_time - start_time
    sec = delta.total_seconds()
    print("difference in seconds:", sec)
    min = sec / 60
    print("difference in minutes:", min)
    return None


def polygon_difference_cython(
    sFilename_base,
    sFilename_new,
    sAttribute_name_base,
    sAttribute_name_new,
    sFilename_output_in,
    dArea_threshold_in=1000.0,
    dBuffer_threshold_in=-0.0001,
):

    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)  # WGS84 lat/lon
    start_time = datetime.now()
    # read the base file
    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    # delete if file exist
    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)

    # create a temp file first by replacing the
    sFilename_tmp = sFilename_output_in.replace(".geojson", "_tmp.geojson")
    if os.path.exists(sFilename_tmp):
        os.remove(sFilename_tmp)

    pDataset_merge_out = pDriver_geojson.CreateDataSource(sFilename_tmp)
    pLayerOut_merge = pDataset_merge_out.CreateLayer(
        "merge", pSpatial_reference_gcs, ogr.wkbPolygon
    )
    pLayerOut_merge.CreateField(
        ogr.FieldDefn("id", ogr.OFTInteger64)
    )  # long type for high resolution
    # Add one attribute
    pLayerDefn_merge = pLayerOut_merge.GetLayerDefn()
    pFeature_merge = ogr.Feature(pLayerDefn_merge)

    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_output_in)
    pLayerOut = pDataset_out.CreateLayer("diff", pSpatial_reference_gcs, ogr.wkbPolygon)
    pLayerOut.CreateField(
        ogr.FieldDefn("polygonid", ogr.OFTInteger64)
    )  # long type for high resolution
    pLayerOut.CreateField(
        ogr.FieldDefn("diff", ogr.OFTReal)
    )  # long type for high resolution
    pLayerOut.CreateField(
        ogr.FieldDefn("area", ogr.OFTReal)
    )  # long type for high resolution
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)
    # read the base file
    pDataset_base = pDriver_geojson.Open(sFilename_base, 0)
    # find the numebr of geometries in the base file
    pLayer_base = pDataset_base.GetLayer()
    index_base = RTreeindex()
    aData_base = list()
    lID = 0
    for pFeature_base in pLayer_base:
        pGeometry_base = pFeature_base.GetGeometryRef()
        left, right, bottom, top = pGeometry_base.GetEnvelope()
        pBound = (left, bottom, right, top)
        index_base.insert(lID, pBound)  #
        dData_base = pFeature_base.GetField(sAttribute_name_base)
        aData_base.append(dData_base)
        lID = lID + 1

    # read the new file
    pDataset_new = pDriver_geojson.Open(sFilename_new, 0)
    pLayer_new = pDataset_new.GetLayer()

    # create a union geometery for the all the intersection polygon
    pGeometry_union = ogr.Geometry(ogr.wkbPolygon)
    lID_polygon = 1
    for pFeature_new in pLayer_new:
        pGeometry_new = pFeature_new.GetGeometryRef()
        left, right, bottom, top = pGeometry_new.GetEnvelope()
        pBound = (left, bottom, right, top)
        aIntersect = list(index_base.search(pBound))
        dData_new = pFeature_new.GetField(sAttribute_name_new)
        nIntersect = len(aIntersect)
        for k in aIntersect:
            pFeature_base = pLayer_base.GetFeature(k)
            pGeometry_base = pFeature_base.GetGeometryRef()
            dData_base = aData_base[k]
            dData_diff = dData_new - dData_base

            if dData_new == -9999 or dData_base == -9999:
                dData_percent = 0.0
            else:
                dData_percent = (dData_diff) / np.max([dData_new, dData_base]) * 100
            if dData_percent > 100:
                dData_percent = 100
            if dData_percent < -100:
                dData_percent = -100

            iFlag_intersect = pGeometry_new.Intersects(pGeometry_base)
            if iFlag_intersect == True:
                pGeometry_intersect = pGeometry_new.Intersection(pGeometry_base)
                pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                iCount = pGeometry_intersect.GetGeometryCount()
                if pGeometrytype_intersect == "POLYGON":
                    aCoords_gcs = pGeometry_intersect
                    dArea = calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1])
                    if dArea > dArea_threshold_in:
                        pFeatureOut.SetGeometry(pGeometry_intersect)
                        pFeatureOut.SetField("polygonid", lID_polygon)
                        pFeatureOut.SetField("diff", dData_diff)
                        pFeatureOut.SetField("perc", dData_percent)
                        pFeatureOut.SetField("area", dArea)
                        pLayerOut.CreateFeature(pFeatureOut)
                        lID_polygon = lID_polygon + 1

                    pGeometry_union = pGeometry_union.Union(pGeometry_intersect)
                else:
                    print("pGeometrytype_intersect:", pGeometrytype_intersect)
                    if pGeometrytype_intersect == "MULTIPOLYGON":
                        for i in range(iCount):
                            pPolygon = pGeometry_intersect.GetGeometryRef(i)
                            aCoords_gcs = get_geometry_coordinates(pPolygon)
                            dArea = calculate_polygon_area(
                                aCoords_gcs[:, 0], aCoords_gcs[:, 1]
                            )
                            if dArea > dArea_threshold_in:
                                pFeatureOut.SetGeometry(pPolygon)
                                pFeatureOut.SetField("polygonid", lID_polygon)
                                pFeatureOut.SetField("diff", dData_diff)
                                pFeatureOut.SetField("perc", dData_percent)
                                pFeatureOut.SetField("area", dArea)
                                pLayerOut.CreateFeature(pFeatureOut)
                                lID_polygon = lID_polygon + 1

                            pGeometry_union = pGeometry_union.Union(pPolygon)

    # get the coordinates of the union
    pPolygon = pGeometry_union.GetGeometryRef(0)
    aCoords_gcs = get_geometry_coordinates(pPolygon)
    pFeature_merge.SetGeometry(pGeometry_union)
    pFeature_merge.SetField("id", 1)
    pLayerOut_merge.CreateFeature(pFeature_merge)
    pDataset_merge_out = None

    iFlag_remaining_new = 0
    iFlag_remaining_base = 0
    # remaining news
    if iFlag_remaining_new == 1:
        nFeature_new = pLayer_new.GetFeatureCount()
        for j in range(nFeature_new):
            pFeature_new = pLayer_new.GetFeature(j)
            pGeometry_new = pFeature_new.GetGeometryRef()
            iFlag_within = pGeometry_new.Within(pGeometry_union)
            iFlag_intersect = pGeometry_new.Intersects(pGeometry_union)
            dData_new = pFeature_new.GetField(sAttribute_name_new)
            dData_percent = 100
            aCoords_gcs = get_geometry_coordinates(pGeometry_new)
            dArea = calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1])
            if iFlag_intersect == True:
                if iFlag_within == True:
                    continue
                # not within as well

                pGeometry_difference = pGeometry_new.Difference(pGeometry_union)
                pGeometrytype_difference = pGeometry_difference.GetGeometryName()
                iCount = pGeometry_difference.GetGeometryCount()
                if pGeometrytype_difference == "POLYGON":
                    for i in range(iCount):
                        pPolygon0 = pGeometry_difference.GetGeometryRef(i)
                        pGeometerytype = pPolygon0.GetGeometryName()
                        if pGeometerytype == "POLYGON":
                            pPolygon = pPolygon0
                            pass
                        else:
                            if pGeometerytype == "LINEARRING":
                                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                                pPolygon.AddGeometry(pPolygon0)
                            else:
                                print(pGeometerytype)
                                continue
                        aCoords_gcs = get_geometry_coordinates(pPolygon)
                        dArea = calculate_polygon_area(
                            aCoords_gcs[:, 0], aCoords_gcs[:, 1]
                        )
                        if dArea > dArea_threshold_in:
                            pPolygon = pPolygon.Buffer(dBuffer_threshold_in)
                            pFeatureOut.SetGeometry(pPolygon)
                            pFeatureOut.SetField("polygonid", lID_polygon)
                            pFeatureOut.SetField("diff", dData_new)
                            pFeatureOut.SetField("perc", dData_percent)
                            pFeatureOut.SetField("area", dArea)
                            pLayerOut.CreateFeature(pFeatureOut)
                            lID_polygon = lID_polygon + 1
                        else:
                            print("Unwanted polygon:")
                else:
                    print("pGeometrytype_difference:", pGeometrytype_difference)
                    if pGeometrytype_difference == "MULTIPOLYGON":
                        for i in range(iCount):
                            pPolygon0 = pGeometry_difference.GetGeometryRef(i)
                            pGeometerytype = pPolygon0.GetGeometryName()
                            if pGeometerytype == "POLYGON":
                                pass
                            else:
                                print(pGeometerytype)
                                continue
                            aCoords_gcs = get_geometry_coordinates(pPolygon0)
                            dArea = calculate_polygon_area(
                                aCoords_gcs[:, 0], aCoords_gcs[:, 1]
                            )
                            if dArea > dArea_threshold_in:
                                pPolygon = pPolygon0.Buffer(dBuffer_threshold_in)
                                pFeatureOut.SetGeometry(pPolygon)
                                pFeatureOut.SetField("polygonid", lID_polygon)
                                pFeatureOut.SetField("diff", dData_new)
                                pFeatureOut.SetField("perc", dData_percent)
                                pFeatureOut.SetField("area", dArea)
                                pLayerOut.CreateFeature(pFeatureOut)
                                lID_polygon = lID_polygon + 1
                            else:
                                print("Unwanted polygon:")

            else:
                if dArea > dArea_threshold_in:
                    pFeatureOut.SetGeometry(pGeometry_new)
                    pFeatureOut.SetField("polygonid", lID_polygon)
                    pFeatureOut.SetField("diff", dData_new)
                    pFeatureOut.SetField("perc", dData_percent)
                    pFeatureOut.SetField("area", dArea)
                    pLayerOut.CreateFeature(pFeatureOut)
                    lID_polygon = lID_polygon + 1

    # remaining base
    if iFlag_remaining_base == 1:
        nFeature_base = pLayer_base.GetFeatureCount()
        for j in range(nFeature_base):
            pFeature_base = pLayer_base.GetFeature(j)
            pGeometry_base = pFeature_base.GetGeometryRef()
            iFlag_within = pGeometry_base.Within(pGeometry_union)
            iFlag_intersect = pGeometry_base.Intersects(pGeometry_union)
            dData_base = pFeature_base.GetField(sAttribute_name_base)

            dData_percent = -100
            aCoords_gcs = get_geometry_coordinates(pGeometry_base)
            dArea = calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1])
            if iFlag_intersect == True:
                if iFlag_within == True:
                    continue
                pGeometry_difference = pGeometry_base.Difference(pGeometry_union)
                pGeometrytype_difference = pGeometry_difference.GetGeometryName()
                iCount = pGeometry_difference.GetGeometryCount()
                if pGeometrytype_difference == "POLYGON":
                    for i in range(iCount):
                        pPolygon0 = pGeometry_difference.GetGeometryRef(i)
                        pGeometerytype = pPolygon0.GetGeometryName()
                        if pGeometerytype == "POLYGON":
                            pPolygon = pPolygon0
                            pass
                        else:
                            if pGeometerytype == "LINEARRING":
                                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                                pPolygon.AddGeometry(pPolygon0)
                            else:
                                print(pGeometerytype)
                                continue
                        aCoords_gcs = get_geometry_coordinates(pPolygon0)
                        dArea = calculate_polygon_area(
                            aCoords_gcs[:, 0], aCoords_gcs[:, 1]
                        )
                        if dArea > dArea_threshold_in:
                            pPolygon = pPolygon.Buffer(dBuffer_threshold_in)
                            pFeatureOut.SetGeometry(pPolygon)
                            pFeatureOut.SetField("polygonid", lID_polygon)
                            pFeatureOut.SetField("diff", -1 * dData_base)
                            pFeatureOut.SetField("perc", dData_percent)
                            pFeatureOut.SetField("area", dArea)
                            pLayerOut.CreateFeature(pFeatureOut)
                            lID_polygon = lID_polygon + 1
                        else:
                            print("Unwanted polygon:")
                else:
                    print("pGeometrytype_difference:", pGeometrytype_difference)
                    if pGeometrytype_difference == "MULTIPOLYGON":
                        for i in range(iCount):
                            pPolygon0 = pGeometry_difference.GetGeometryRef(i)
                            pGeometerytype = pPolygon0.GetGeometryName()
                            if pGeometerytype == "POLYGON":
                                pass
                            else:
                                print(pGeometerytype)
                                continue

                            aCoords_gcs = get_geometry_coordinates(pPolygon0)
                            dArea = calculate_polygon_area(
                                aCoords_gcs[:, 0], aCoords_gcs[:, 1]
                            )
                            if dArea > dArea_threshold_in:
                                pPolygon = pPolygon0.Buffer(dBuffer_threshold_in)
                                pFeatureOut.SetGeometry(pPolygon)
                                pFeatureOut.SetField("polygonid", lID_polygon)
                                pFeatureOut.SetField("diff", -1 * dData_base)
                                pFeatureOut.SetField("perc", dData_percent)
                                pFeatureOut.SetField("area", dArea)
                                pLayerOut.CreateFeature(pFeatureOut)
                                lID_polygon = lID_polygon + 1
                            else:
                                print("Unwanted polygon:")

            else:
                if dArea > dArea_threshold_in:
                    pFeatureOut.SetGeometry(pGeometry_base)
                    pFeatureOut.SetField("polygonid", lID_polygon)
                    pFeatureOut.SetField("diff", -1 * dData_base)
                    pFeatureOut.SetField("perc", dData_percent)
                    pFeatureOut.SetField("area", dArea)
                    pLayerOut.CreateFeature(pFeatureOut)
                    lID_polygon = lID_polygon + 1

    # close files
    pDataset_base = None
    pDataset_new = None
    pSpatial_reference_gcs = None
    end_time = datetime.now()
    # get difference
    delta = end_time - start_time
    sec = delta.total_seconds()
    print("difference in seconds:", sec)
    min = sec / 60
    print("difference in minutes:", min)

    return None


def polygon_difference_cython_channel(
    sFilename_base,
    sFilename_new,
    sAttribute_name_base,
    sAttribute_name_new,
    aChannel_base,
    aChannel_new,
    sFilename_output_in,
    dArea_threshold_in=1000.0,
    dBuffer_threshold_in=-0.0001,
):

    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)  # WGS84 lat/lon
    start_time = datetime.now()
    # read the base file
    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    # delete if file exist
    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)

    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_output_in)
    pLayerOut = pDataset_out.CreateLayer("diff", pSpatial_reference_gcs, ogr.wkbPolygon)
    pLayerOut.CreateField(
        ogr.FieldDefn("polygonid", ogr.OFTInteger64)
    )  # long type for high resolution
    pLayerOut.CreateField(ogr.FieldDefn("diff", ogr.OFTReal))
    pLayerOut.CreateField(ogr.FieldDefn("perc", ogr.OFTReal))
    pLayerOut.CreateField(ogr.FieldDefn("area", ogr.OFTReal))
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)
    # read the base file
    pDataset_base = pDriver_geojson.Open(sFilename_base, 0)
    # find the numebr of geometries in the base file
    pLayer_base = pDataset_base.GetLayer()
    index_base = RTreeindex()
    aData_base = list()
    lID = 0
    nFeature_base = pLayer_base.GetFeatureCount()
    for j in range(nFeature_base):
        # for pFeature_base in pLayer_base:
        pFeature_base = pLayer_base.GetFeature(j)
        iFlag_channel = aChannel_base[j]
        # if iFlag_channel ==1:
        pGeometry_base = pFeature_base.GetGeometryRef()
        left, right, bottom, top = pGeometry_base.GetEnvelope()
        pBound = (left, bottom, right, top)
        index_base.insert(lID, pBound)  #
        dData_base = pFeature_base.GetField(sAttribute_name_base)
        aData_base.append(dData_base)
        lID = lID + 1

    # read the new file
    pDataset_new = pDriver_geojson.Open(sFilename_new, 0)
    pLayer_new = pDataset_new.GetLayer()

    # create a union geometery for the all the intersection polygon
    pGeometry_union = ogr.Geometry(ogr.wkbPolygon)
    lID_polygon = 1
    nFeature_new = pLayer_new.GetFeatureCount()
    # for pFeature_new in pLayer_new:
    for j in range(nFeature_new):
        pFeature_new = pLayer_new.GetFeature(j)
        iFlag_channel_new = aChannel_new[j]
        if iFlag_channel_new == 1:  # check whethere it is a main channel
            pGeometry_new = pFeature_new.GetGeometryRef()
            left, right, bottom, top = pGeometry_new.GetEnvelope()
            pBound = (left, bottom, right, top)
            aIntersect = list(index_base.search(pBound))
            dData_new = pFeature_new.GetField(sAttribute_name_new)
            nIntersect = len(aIntersect)
            for k in aIntersect:
                pFeature_base = pLayer_base.GetFeature(k)
                iFlag_channel_base = aChannel_base[k]
                if iFlag_channel_base == 1:
                    pGeometry_base = pFeature_base.GetGeometryRef()
                    dData_base = aData_base[k]
                    dData_diff = dData_new - dData_base

                    if dData_new == -9999 or dData_base == -9999:
                        dData_percent = 0.0
                    else:
                        dData_percent = (
                            (dData_diff) / np.max([dData_new, dData_base]) * 100
                        )
                    if dData_percent >= 100:
                        dData_percent = 100
                    if dData_percent < -100:
                        dData_percent = -100

                    iFlag_intersect = pGeometry_new.Intersects(pGeometry_base)
                    if iFlag_intersect == True:
                        pGeometry_intersect = pGeometry_new.Intersection(pGeometry_base)
                        pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                        iCount = pGeometry_intersect.GetGeometryCount()
                        if pGeometrytype_intersect == "POLYGON":
                            aCoords_gcs = get_geometry_coordinates(pGeometry_intersect)
                            dArea = calculate_polygon_area(
                                aCoords_gcs[:, 0], aCoords_gcs[:, 1]
                            )
                            if dArea > dArea_threshold_in:
                                pFeatureOut.SetGeometry(pGeometry_intersect)
                                pFeatureOut.SetField("polygonid", lID_polygon)
                                pFeatureOut.SetField("diff", dData_diff)
                                pFeatureOut.SetField("perc", dData_percent)
                                pFeatureOut.SetField("area", dArea)
                                pLayerOut.CreateFeature(pFeatureOut)
                                lID_polygon = lID_polygon + 1

                        else:
                            print("pGeometrytype_intersect:", pGeometrytype_intersect)
                            if pGeometrytype_intersect == "MULTIPOLYGON":
                                for i in range(iCount):
                                    pPolygon = pGeometry_intersect.GetGeometryRef(i)
                                    aCoords_gcs = get_geometry_coordinates(pPolygon)
                                    dArea = calculate_polygon_area(
                                        aCoords_gcs[:, 0], aCoords_gcs[:, 1]
                                    )
                                    if dArea > dArea_threshold_in:
                                        pFeatureOut.SetGeometry(pPolygon)
                                        pFeatureOut.SetField("polygonid", lID_polygon)
                                        pFeatureOut.SetField("diff", dData_diff)
                                        pFeatureOut.SetField("perc", dData_percent)
                                        pFeatureOut.SetField("area", dArea)
                                        pLayerOut.CreateFeature(pFeatureOut)
                                        lID_polygon = lID_polygon + 1

    # close files
    pDataset_base = None
    pDataset_new = None
    pSpatial_reference_gcs = None
    end_time = datetime.now()
    # get difference
    delta = end_time - start_time
    sec = delta.total_seconds()
    print("difference in seconds:", sec)
    min = sec / 60
    print("difference in minutes:", min)

    return None
