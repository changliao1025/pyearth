import os
import numpy as np
from osgeo import ogr, gdal, osr
from pyearth.system.define_global_variables import *
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_intersect_on_great_circle import find_great_circle_intersection

def split_polygon_cross_idl(aCoord_gcs):
    #find out the two index where the edge crosses the IDL
    nPoint = len(aCoord_gcs)
    aIndex = []
    #this algo require the coordinates to be in CCW order

    for i in range(nPoint - 1):
        dLongitude = aCoord_gcs[i,0]
        dLongitude_next = aCoord_gcs[i + 1,0]
        if dLongitude > 0 and  dLongitude < 180.0 and dLongitude_next < 0:
            aIndex.append(i)
            continue
        if dLongitude < 0 and dLongitude_next > 0:
            aIndex.append(i)
            continue

    #check the closing edge
    dLongitude = aCoord_gcs[nPoint - 1,0]
    dLongitude_next = aCoord_gcs[0,0]
    if dLongitude > 0 and  dLongitude < 180.0 and dLongitude_next < 0:
        aIndex.append(nPoint - 1)
    if dLongitude < 0 and dLongitude_next > 0:
        aIndex.append(nPoint - 1)

    if len(aIndex) != 2:
        print('Warning: no intersection found')
        print(aCoord_gcs)
        return

    #get the two intersection points
    lon1= aCoord_gcs[aIndex[0],0]
    lat1= aCoord_gcs[aIndex[0],1]
    lon2= aCoord_gcs[aIndex[0] + 1,0]
    lat2= aCoord_gcs[aIndex[0] + 1,1]
    target_lon = 180.0
    d, dLat0 = find_great_circle_intersection(lon1, lat1, lon2, lat2, target_lon)

    lon1= aCoord_gcs[aIndex[1],0]
    lat1= aCoord_gcs[aIndex[1],1]
    if aIndex[1] == nPoint - 1:
        lon2= aCoord_gcs[0,0]
        lat2= aCoord_gcs[0,1]
    else:
        lon2= aCoord_gcs[aIndex[1] + 1,0]
        lat2= aCoord_gcs[aIndex[1] + 1,1]
    target_lon = 180.0
    d, dLat1 = find_great_circle_intersection(lon1, lat1, lon2, lat2, target_lon)

    #compare the two intersection points which is top and bottom
    if dLat0 > dLat1:
        dLat_top = dLat0
        dLat_bottom = dLat1
    else:
        dLat_top = dLat1
        dLat_bottom = dLat0

    aCoord_gcs_left = list()
    aCoord_gcs_right = list()

    #left part dLongitude > 0
    aIndex_dummy = np.array(aIndex)
    iFlag_added_right = 0
    iFlag_added_left = 0
    for i in range(nPoint):
        dLongitude = aCoord_gcs[i,0]
        if dLongitude > 0:
            if i <= np.min(aIndex_dummy):
                aCoord_gcs_left.append(aCoord_gcs[i])
                if i in aIndex:
                    if iFlag_added_left == 0:
                        iFlag_added_left = 1
                        aCoord_gcs_left.append([180-1.0E-10, dLat_bottom])
                        aCoord_gcs_left.append([180-1.0E-10, dLat_top])
            else:
                if i in aIndex:
                    aCoord_gcs_left.append(aCoord_gcs[i])
                    if iFlag_added_left == 0:
                        iFlag_added_left = 1
                        aCoord_gcs_left.append([180-1.0E-10, dLat_bottom])
                        aCoord_gcs_left.append([180-1.0E-10, dLat_top])
                else:
                    aCoord_gcs_left.append(aCoord_gcs[i])

        if dLongitude < 0:
            if i <= np.min(aIndex_dummy):
                aCoord_gcs_right.append(aCoord_gcs[i])
                if i in aIndex:
                    if iFlag_added_right == 0:
                        iFlag_added_right = 1
                        aCoord_gcs_right.append([-180+1.0E-10, dLat_top])
                        aCoord_gcs_right.append([-180+1.0E-10, dLat_bottom])
                    else:
                        pass
            else:
                if i in aIndex:
                    aCoord_gcs_right.append(aCoord_gcs[i])
                    if iFlag_added_right == 0:
                        aCoord_gcs_right.append([-180+1.0E-10, dLat_top])
                        aCoord_gcs_right.append([-180+1.0E-10, dLat_bottom])
                else:
                    aCoord_gcs_right.append(aCoord_gcs[i])

    return [aCoord_gcs_left, aCoord_gcs_right]