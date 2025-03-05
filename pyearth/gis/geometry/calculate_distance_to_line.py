import numpy as np
from pyearth.gis.location.convert_longitude_latitude_to_sphere_3d import convert_longitude_latitude_to_sphere_3d

def calculate_distance_to_line(dLongitude1_in, dLatitude1_in,
                                dLongitude2_in, dLatitude2_in,
                                dLongitude3_in, dLatitude3_in,
                                iFlag_radian = None):
    if iFlag_radian is None:
        dLongitude1_radian_in, dLatitude1_radian_in = np.radians(np.array((dLongitude1_in, dLatitude1_in) ))
        dLongitude2_radian_in, dLatitude2_radian_in = np.radians(np.array((dLongitude2_in, dLatitude2_in) )) #this is the middle one
        dLongitude3_radian_in, dLatitude3_radian_in = np.radians(np.array((dLongitude3_in, dLatitude3_in) ))
        pass
    else:
        dLongitude1_radian_in, dLatitude1_radian_in = dLongitude1_in , dLatitude1_in
        dLongitude2_radian_in, dLatitude2_radian_in = dLongitude2_in , dLatitude2_in #this is the middle one
        dLongitude3_radian_in, dLatitude3_radian_in = dLongitude3_in , dLatitude3_in
        pass

    # The points in 3D space
    x1,y1,z1 = convert_longitude_latitude_to_sphere_3d(dLongitude1_radian_in, dLatitude1_radian_in)
    x2,y2,z2 = convert_longitude_latitude_to_sphere_3d(dLongitude2_radian_in, dLatitude2_radian_in)
    x3,y3,z3 = convert_longitude_latitude_to_sphere_3d(dLongitude3_radian_in, dLatitude3_radian_in)

    #Calculate the perpendicular distance from a point to a line segment in 3d
    #https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    # Vectors from line_start to line_end and from line_start to the point
    #x1 and x3 are the two points on the line
    line_vec = np.array([x3 - x1, y3 - y1, z3 - z1])
    point_vec = np.array([x2 - x1, y2 - y1, z2 - z1])

    # Project point_vec onto line_vec to find the projection point on the line
    line_len = np.linalg.norm(line_vec)
    line_unitvec = line_vec / line_len
    projection_length = np.dot(point_vec, line_unitvec)

    # Clamp the projection length to the length of the line segment
    projection_length = max(0, min(line_len, projection_length))

    # Find the projection point on the line segment
    projection_point = np.array([x1, y1, z1]) + projection_length * line_unitvec

    # Calculate the distance from the point to the projection point
    distance = np.linalg.norm(np.array([x2, y2, z2]) - projection_point)

    #consider the radius of the earth in m
    distance = distance * 6371000.0

    return float(distance)