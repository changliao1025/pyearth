import numpy as np
from osgeo import ogr
def check_ccw(coords):
    """
    Determines if a polygon's coordinates are in counter-clockwise (CCW) order.

    Args:
        coords: A numpy array of shape (n, 2) representing the polygon's coordinates
                in (longitude, latitude) format.

    Returns:
        bool: True if the coordinates are CCW, False if they are CW.
    """
    # Shoelace formula to calculate the signed area
    x = coords[:, 0]  # Longitude
    y = coords[:, 1]  # Latitude
    signed_area = np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]) + (x[-1] * y[0] - x[0] * y[-1])

    # If the signed area is positive, the coordinates are CCW
    return signed_area > 0

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
    # Convert to numpy array
    coords_array = np.array(exterior_coords)
    # Check if the coordinates are CCW
    if not check_ccw(coords_array):
        # Reverse the order if not CCW
        coords_array = coords_array[::-1]

    return coords_array

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
