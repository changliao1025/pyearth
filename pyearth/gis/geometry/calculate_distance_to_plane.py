import math
import numpy as np
from pyearth.gis.location.convert_longitude_latitude_to_sphere_3d import convert_longitude_latitude_to_sphere_3d

def calculate_distance_to_plane_new(dLongitude1_in, dLatitude1_in,
                                dLongitude2_in, dLatitude2_in,
                                dLongitude3_in, dLatitude3_in,
                                iFlag_radian = None):
    """Calculate the distance of a point to a plane defined by three points in 3D space."""
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
    # Convert the three points to 3D coordinates
    # The points in 3D space
    x1,y1,z1 = convert_longitude_latitude_to_sphere_3d(dLongitude1_radian_in, dLatitude1_radian_in)
    x2,y2,z2 = convert_longitude_latitude_to_sphere_3d(dLongitude2_radian_in, dLatitude2_radian_in)
    x3,y3,z3 = convert_longitude_latitude_to_sphere_3d(dLongitude3_radian_in, dLatitude3_radian_in)

    # Calculate two vectors on the plane
    v1_x, v1_y, v1_z = x2 - x1, y2 - y1, z2 - z1
    v2_x, v2_y, v2_z = x3 - x1, y3 - y1, z3 - z1

    # Compute the normal vector using the cross product
    normal_x = v1_y * v2_z - v1_z * v2_y
    normal_y = v1_z * v2_x - v1_x * v2_z
    normal_z = v1_x * v2_y - v1_y * v2_x

    # Check if the normal vector is zero (points are collinear)
    if abs(normal_x) < 1e-10 and abs(normal_y) < 1e-10 and abs(normal_z) < 1e-10:
        distance = 0.0
        return distance

    # Calculate the plane equation coefficients (A, B, C, D)
    A, B, C = normal_x, normal_y, normal_z
    D = -(A * x1 + B * y1 + C * z1)

    # Calculate the distance of the second point to the plane
    distance = abs(A * x2 + B * y2 + C * z2 + D) / math.sqrt(A**2 + B**2 + C**2)

    return distance


#https://stackoverflow.com/questions/8204998/how-to-check-if-a-pointlonc-latc-lie-on-a-great-circle-running-from-lona-lata

def calculate_distance_to_plane_old(dLongitude1_in, dLatitude1_in,
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
    #The formula is x+b*y+c*z=0

    # Check for zero denominators
    denominator_c = z1 * y3 - z3 * y1
    denominator_b = y1 * z3 - y3 * z1

    if denominator_c == 0 or denominator_b == 0:
        print(z1, y3, z3, y1)
        print(y1, z3, y3, z1)
        raise ValueError("Division by zero encountered in the calculation of coefficients 'b' and 'c'. Check input points.")

    # Calculate coefficients
    c = (-x1 * y3 + x3 * y1) / denominator_c
    b = (-x1 * z3 + x3 * z1) / denominator_b
    distance = abs(  x2 + b * y2 + c * z2 )
    return float(distance)

