import math
import numpy as np
from scipy.optimize import fsolve
from pyearth.gis.location.convert_longitude_latitude_to_sphere_3d import convert_longitude_latitude_to_sphere_3d

def deg_to_rad(deg):
    """Converts degrees to radians."""
    return deg * np.pi / 180

def rad_to_deg(rad):
    """Converts radians to degrees."""
    return rad * 180 / np.pi

def latlon_to_cartesian(lon, lat ):
    """Converts latitude and longitude to Cartesian coordinates on the unit sphere."""
    lat, lon = deg_to_rad(lat), deg_to_rad(lon)
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.array([x, y, z])

def convert_sphere_3d_to_longitude_latitude(x, y, z):
    longitude = np.arctan2(y, x)
    latitude = np.arcsin(z / np.sqrt(x**2 + y**2 + z**2))
    return np.degrees(longitude), np.degrees(latitude)

def great_circle_param(t, p1, p2):
    """Parameterizes the great circle."""
    theta = (1 - t) * np.arctan2(p1[1], p1[0]) + t * np.arctan2(p2[1], p2[0])
    phi = np.arcsin((1 - t) * p1[2] + t * p2[2])
    return np.array([np.cos(theta) * np.cos(phi), np.sin(theta) * np.cos(phi), np.sin(phi)])

def project_point_onto_plane(p, normal):
    """
    Project a point onto a plane defined by a normal vector.

    Args:
        p (np.array): Point to be projected.
        normal (np.array): Normal vector of the plane.

    Returns:
        np.array: Projected point on the plane.
    """
    # Project point onto the plane
    projection = p - np.dot(p, normal) * normal
    return projection

def great_circle_intersection(p1, p2, p3, v):
    """
    Finds the intersection point(s) of a great circle and a line in 3D.

    Args:
        p1: A 3D point on the great circle.
        p2: Another 3D point on the great circle.
        p3: A point on the line.
        v: The direction vector of the line.

    Returns:
        A list of intersection points, or an empty list if no intersection.
    """
    def line_param(t, p3, v):
        return p3 + t * v

    def find_intersection(t, p1, p2, p3, v):
        p_gc = great_circle_param(t, p1, p2)
        p_line = line_param(t, p3, v)
        return np.linalg.norm(p_gc - p_line)

    # Convert points to unit vectors
    p1 = p1 / np.linalg.norm(p1)
    p2 = p2 / np.linalg.norm(p2)
    p3 = p3 / np.linalg.norm(p3)
    v = v / np.linalg.norm(v)

    # Calculate normal vectors to the planes defined by the great circles
    normal1 = np.cross(p1, p2)
    normal1 = normal1 / np.linalg.norm(normal1)

    # Project p3 onto the plane defined by normal1
    p3_proj = project_point_onto_plane(p3, normal1)
    p3_proj = p3_proj / np.linalg.norm(p3_proj)

    #t_intersect = fsolve(find_intersection, 0.5, args=(p1, p2, p3, v))
    #intersection_point = line_param(t_intersect, p3, v)
    return p3_proj

def calculate_intersect_on_great_circle(dLongitude1_in, dLatitude1_in,
                                        dLongitude2_in, dLatitude2_in,
                                        dLongitude3_in, dLatitude3_in,
                                        iFlag_radian=None):

    if iFlag_radian is None:
        dLongitude1_radian_in, dLatitude1_radian_in = np.radians([dLongitude1_in, dLatitude1_in])
        dLongitude2_radian_in, dLatitude2_radian_in = np.radians([dLongitude2_in, dLatitude2_in])
        dLongitude3_radian_in, dLatitude3_radian_in = np.radians([dLongitude3_in, dLatitude3_in])
    else:
        dLongitude1_radian_in, dLatitude1_radian_in = dLongitude1_in, dLatitude1_in
        dLongitude2_radian_in, dLatitude2_radian_in = dLongitude2_in, dLatitude2_in
        dLongitude3_radian_in, dLatitude3_radian_in = dLongitude3_in, dLatitude3_radian_in

    pA = latlon_to_cartesian( dLongitude1_in,dLatitude1_in)
    p2_proj = latlon_to_cartesian( dLongitude2_in,dLatitude2_in)
    pC = latlon_to_cartesian( dLongitude3_in, dLatitude3_in)
    p0 = np.array([0, 0, 0])
    v1 = pA - p0
    v2 = pC - p0
    n = np.cross(v1, v2)
    n /= np.linalg.norm(n)
    d = np.dot(p2_proj, n)
    p2_proj -= d * n
    AC = pC - pA
    direction_vector = np.cross(AC, n)
    direction_vector /= np.linalg.norm(direction_vector)
    intersection_point = great_circle_intersection(pA, pC, p2_proj, direction_vector)
    #intersection_point /= np.linalg.norm(intersection_point)
    longitude_intersect, latitude_intersect = convert_sphere_3d_to_longitude_latitude(*intersection_point)
    return float(longitude_intersect), float(latitude_intersect)

def find_great_circle_intersection(lon1, lat1, lon2, lat2, target_lon):
    """
    Find the location on the great circle that has the specified longitude.

    Args:
        lon1, lat1: Longitude and latitude of the first point (in degrees).
        lon2, lat2: Longitude and latitude of the second point (in degrees).
        target_lon: The target longitude (in degrees).

    Returns:
        (target_lon, target_lat): The longitude and latitude of the intersection point (in degrees).
    """
    # Convert coordinates to radians
    lon1, lat1, lon2, lat2, target_lon = map(deg_to_rad, [lon1, lat1, lon2, lat2, target_lon])

    # Calculate the difference in longitudes
    d_lon = lon2 - lon1

    # Calculate the latitude of the intersection point using spherical interpolation
    Bx = math.cos(lat2) * math.cos(d_lon)
    By = math.cos(lat2) * math.sin(d_lon)
    lat_intersection = math.atan2(math.sin(lat1) + math.sin(lat2), math.sqrt((math.cos(lat1) + Bx) ** 2 + By ** 2))

    # Calculate the longitude of the intersection point
    lon_intersection = target_lon

    # Convert the intersection point back to degrees
    lon_intersection, lat_intersection = map(rad_to_deg, [lon_intersection, lat_intersection])

    return lon_intersection, lat_intersection

if __name__ == '__main__':
    dLongitude1_in = -78.56458333333426
    dLatitude1_in = 53.73124999999924
    dLongitude2_in = -78.568386
    dLatitude2_in = 53.734444
    dLongitude3_in = -78.56875000000096
    dLatitude3_in = 53.735416666665905
    iFlag_radian = None
    dLongitude_out, dLatitude_out = calculate_intersect_on_great_circle(dLongitude1_in, dLatitude1_in,
                                                                        dLongitude2_in, dLatitude2_in,
                                                                        dLongitude3_in, dLatitude3_in,
                                                                        iFlag_radian)
    print(dLongitude_out, dLatitude_out)
