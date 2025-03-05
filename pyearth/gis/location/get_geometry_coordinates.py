import numpy as np
from osgeo import ogr
def get_geometry_coordinates(geometry):

    sGeometry_type = geometry.GetGeometryName()
    if sGeometry_type =='POINT':
        return get_point_coords(geometry)
    elif sGeometry_type == 'LINESTRING':
        return get_linestring_coords(geometry)
    elif sGeometry_type =='POLYGON':
        return get_polygon_exterior_coords(geometry)
    elif sGeometry_type =='LINEARRING':
        return get_linearring_coords(geometry)
    elif sGeometry_type =='MULTIPOLYGON':
        return get_multipolygon_exterior_coords(geometry)
    else:
        print(sGeometry_type)
        raise ValueError("Unsupported geometry type.")

def get_polygon_exterior_coords(polygon_geometry):
    exterior_coords = []
    ring = polygon_geometry.GetGeometryRef(0)  # Get the exterior ring
    npoints = ring.GetPointCount()
    for i in range(npoints):
        point = ring.GetPoint(i)
        exterior_coords.append((point[0], point[1]))
    return np.array(exterior_coords)

def get_multipolygon_exterior_coords(polygon_geometry):
    aExterior_coords = list()
    sGeometry_type = polygon_geometry.GetGeometryName()
    #print(polygon_geometry.GetGeometryType())
    if sGeometry_type == 'MULTIPOLYGON':
        nPart = polygon_geometry.GetGeometryCount()
        for i in range(nPart):
            exterior_coords = []
            polygon = polygon_geometry.GetGeometryRef(i)
            ring = polygon.GetGeometryRef(0)  # Get the exterior ring of each polygon
            npoints = ring.GetPointCount()
            for j in range(npoints):
                point = ring.GetPoint(j)
                exterior_coords.append((point[0], point[1]))

            aExterior_coords.append(exterior_coords)
        return aExterior_coords

def get_linestring_coords(linestring_geometry):
    coords = []
    for i in range(linestring_geometry.GetPointCount()):
        point = linestring_geometry.GetPoint(i)
        coords.append((point[0], point[1]))
    return np.array(coords)

def get_point_coords(point_geometry):
    point = point_geometry.GetPoint()
    return np.array([(point[0], point[1])])

def get_linearring_coords(linearring_geometry):
    coords = []
    for i in range(linearring_geometry.GetPointCount()):
        point = linearring_geometry.GetPoint(i)
        coords.append((point[0], point[1]))
    return np.array(coords)
