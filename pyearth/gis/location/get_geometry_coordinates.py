import numpy as np
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
    else:
        raise ValueError("Unsupported geometry type.")

def get_polygon_exterior_coords(polygon_geometry):
    exterior_coords = []
    ring = polygon_geometry.GetGeometryRef(0)  # Get the exterior ring
    for i in range(ring.GetPointCount()):
        point = ring.GetPoint(i)
        exterior_coords.append((point[0], point[1]))
    return np.array(exterior_coords)

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
