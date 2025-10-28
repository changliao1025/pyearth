

import math
from typing import List, Tuple

from osgeo import ogr

# Define the coordinate type alias for clarity (Lon, Lat)
Coord = Tuple[float, float]



def reorder_idl_polygon(vertices: List[Coord]) -> List[Coord]:

    # Handle empty input
    if not vertices:
        return []

    # Step 1: Strip the closing vertex for processing
    # The closing vertex (duplicate of first vertex) will be re-added at the end
    # to ensure the reordered polygon is properly closed
    # Ensure polygon is closed; append first vertex if not closed.
    if vertices and (vertices[0] != vertices[-1]):
        vertices = vertices + [vertices[0]]

    unique_vertices = vertices[:-1]
    nPoint = len(unique_vertices)
    aVertex_current = unique_vertices.copy()

    #use a rotating method to find the polygon that is valid


    iFlag_valid = False
    iCount_rotate = 0
    while iFlag_valid == False:
        iCount_rotate += 1
        if iCount_rotate > nPoint:
            break
        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pLinearRing = ogr.Geometry(ogr.wkbLinearRing)
        for lon, lat in aVertex_current:
            pLinearRing.AddPoint(lon, lat)
        pLinearRing.CloseRings()
        pPolygon.AddGeometry(pLinearRing)
        if pPolygon.IsValid():
            iFlag_valid = True
            break
        else:
            #rotate the first point to the end
            aVertex_current = aVertex_current[1:] + [aVertex_current[0]]

    if iFlag_valid == False:
        print('Error: cannot form a valid polygon by rotating vertices')
        return vertices
    else:
        new_vertices = aVertex_current + [aVertex_current[0]]

    return new_vertices