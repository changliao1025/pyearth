
cimport cython
from libc.math cimport M_PI, sin, cos, asin, acos, sqrt, abs
cimport numpy as cnp
import numpy as np
from pyearth.gis.location.kernel cimport convert_longitude_latitude_to_sphere_3d

""" Low-level function for pyflowline
"""
# Authors: Chang Liao

#constant

cdef double dRadius = 6378137.0

@cython.boundscheck(False)
cpdef calculate_distance_based_on_longitude_latitude(
    dLongitude1_in, dLatitude1_in,
    dLongitude2_in, dLatitude2_in,
    bFlag_radian=False
):
    """
    Calculate the great circle distance between two points on the earth.
    Automatically handles both scalar and array inputs.

    Parameters
    ----------
    dLongitude1_in, dLatitude1_in, dLongitude2_in, dLatitude2_in : double or array-like
        Coordinates of the two points. In degrees by default, or radians if bFlag_radian=True.
        Can be scalars or arrays.
    bFlag_radian : bint, optional
        If True, input coordinates are in radians. If False (default), input is in degrees.

    Returns
    -------
    double or np.ndarray
        Great circle distance in meters. Returns scalar if inputs are scalars, array if inputs are arrays.
    """
    # Check if inputs are arrays
    try:
        # Try to access shape attribute to detect arrays
        if (hasattr(dLongitude1_in, 'shape') or hasattr(dLatitude1_in, 'shape') or
            hasattr(dLongitude2_in, 'shape') or hasattr(dLatitude2_in, 'shape')):
            # Convert to numpy arrays and use the numpy version
            import numpy as np
            lon1_arr = np.asarray(dLongitude1_in, dtype=np.float64)
            lat1_arr = np.asarray(dLatitude1_in, dtype=np.float64)
            lon2_arr = np.asarray(dLongitude2_in, dtype=np.float64)
            lat2_arr = np.asarray(dLatitude2_in, dtype=np.float64)
            return calculate_distance_based_on_longitude_latitude_numpy(
                lon1_arr, lat1_arr, lon2_arr, lat2_arr, bFlag_radian)

        # If we can iterate over the inputs, they might be lists/tuples
        try:
            iter(dLongitude1_in)
            # Convert to numpy arrays and use the numpy version
            import numpy as np
            lon1_arr = np.asarray(dLongitude1_in, dtype=np.float64)
            lat1_arr = np.asarray(dLatitude1_in, dtype=np.float64)
            lon2_arr = np.asarray(dLongitude2_in, dtype=np.float64)
            lat2_arr = np.asarray(dLatitude2_in, dtype=np.float64)
            return calculate_distance_based_on_longitude_latitude_numpy(
                lon1_arr, lat1_arr, lon2_arr, lat2_arr, bFlag_radian)
        except TypeError:
            pass  # Not iterable, proceed with scalar calculation

    except:
        pass  # Any error, proceed with scalar calculation

    # Scalar calculation
    cdef double lon1, lat1, lon2, lat2
    cdef double dLongitude1_scalar = <double>dLongitude1_in
    cdef double dLatitude1_scalar = <double>dLatitude1_in
    cdef double dLongitude2_scalar = <double>dLongitude2_in
    cdef double dLatitude2_scalar = <double>dLatitude2_in

    if not bFlag_radian:
        lon1 = dLongitude1_scalar / 180.0 * M_PI
        lat1 = dLatitude1_scalar / 180.0 * M_PI
        lon2 = dLongitude2_scalar / 180.0 * M_PI
        lat2 = dLatitude2_scalar / 180.0 * M_PI
    else:
        lon1 = dLongitude1_scalar
        lat1 = dLatitude1_scalar
        lon2 = dLongitude2_scalar
        lat2 = dLatitude2_scalar
    cdef double dLongtitude_diff = lon2 - lon1
    cdef double dLatitude_diff = lat2 - lat1
    cdef double a = sin(dLatitude_diff/2)*sin(dLatitude_diff/2) + cos(lat1) * cos(lat2) * sin(dLongtitude_diff/2)*sin(dLongtitude_diff/2)
    cdef double b = 2 * asin(sqrt(a))
    cdef double c = b * dRadius
    return c

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float64_t, ndim=1] calculate_distance_based_on_longitude_latitude_array(
    double[:] dLongitude1_in, double[:] dLatitude1_in,
    double[:] dLongitude2_in, double[:] dLatitude2_in,
    bint bFlag_radian=False
):
    """
    Calculate the great circle distance between arrays of points on the earth.

    Parameters
    ----------
    dLongitude1_in, dLatitude1_in, dLongitude2_in, dLatitude2_in : double[:]
        Arrays of coordinates. In degrees by default, or radians if bFlag_radian=True.
    bFlag_radian : bint, optional
        If True, input coordinates are in radians. If False (default), input is in degrees.

    Returns
    -------
    np.ndarray
        Array of distances in meters.
    """
    cdef int i
    cdef int n = dLongitude1_in.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1] result = np.zeros(n, dtype=np.float64)
    cdef double lon1, lat1, lon2, lat2
    cdef double dLongtitude_diff, dLatitude_diff
    cdef double a, b, c
    for i in range(n):
        if not bFlag_radian:
            lon1 = dLongitude1_in[i] / 180.0 * M_PI
            lat1 = dLatitude1_in[i] / 180.0 * M_PI
            lon2 = dLongitude2_in[i] / 180.0 * M_PI
            lat2 = dLatitude2_in[i] / 180.0 * M_PI
        else:
            lon1 = dLongitude1_in[i]
            lat1 = dLatitude1_in[i]
            lon2 = dLongitude2_in[i]
            lat2 = dLatitude2_in[i]
        dLongtitude_diff = lon2 - lon1
        dLatitude_diff = lat2 - lat1
        a = sin(dLatitude_diff/2)*sin(dLatitude_diff/2) + cos(lat1) * cos(lat2) * sin(dLongtitude_diff/2)*sin(dLongtitude_diff/2)
        b = 2 * asin(sqrt(a))
        c = b * dRadius
        result[i] = c
    return result

# For compatibility, also provide a version that works with NumPy arrays directly
@cython.boundscheck(False)
@cython.wraparound(False)
def calculate_distance_based_on_longitude_latitude_numpy(
    dLongitude1_in, dLatitude1_in, dLongitude2_in, dLatitude2_in, bFlag_radian=False):
    """
    Calculate the great circle distance between arrays of points using NumPy arrays.

    Parameters
    ----------
    dLongitude1_in, dLatitude1_in, dLongitude2_in, dLatitude2_in : np.ndarray
        Arrays of coordinates. In degrees by default, or radians if bFlag_radian=True.
    bFlag_radian : bool, optional
        If True, input coordinates are in radians. If False (default), input is in degrees.

    Returns
    -------
    np.ndarray
        Array of distances in meters.
    """
    # Convert to radians if needed
    if not bFlag_radian:
        dLongitude_radian1_in = dLongitude1_in / 180.0 * M_PI
        dLatitude_radian1_in = dLatitude1_in / 180.0 * M_PI
        dLongitude_radian2_in = dLongitude2_in / 180.0 * M_PI
        dLatitude_radian2_in = dLatitude2_in / 180.0 * M_PI
    else:
        dLongitude_radian1_in = dLongitude1_in
        dLatitude_radian1_in = dLatitude1_in
        dLongitude_radian2_in = dLongitude2_in
        dLatitude_radian2_in = dLatitude2_in
    dLongtitude_diff = dLongitude_radian2_in - dLongitude_radian1_in
    dLatitude_diff = dLatitude_radian2_in - dLatitude_radian1_in
    a = np.sin(dLatitude_diff/2)**2 + np.cos(dLatitude_radian1_in) * np.cos(dLatitude_radian2_in) * np.sin(dLongtitude_diff/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) * dRadius
    return c

@cython.boundscheck(False)  # deactivate bnds checking
cpdef calculate_angle_between_vectors_coordinates(double x1, double y1, double z1, double x2, double y2, double z2):
    """Return the angle between two vectors in any dimension space,
    in degrees.
    """
    cdef double a, b, c, d, e, f
    a = x1*x2 + y1*y2 + z1*z2
    b = sqrt( x1*x1 + y1*y1 + z1*z1   )
    c = sqrt( x2*x2 + y2*y2 + z2*z2   )
    if b == 0 or c == 0:
        raise ValueError("Zero-length vector encountered. Cannot calculate angle.")
    d = a / (b* c)
    # Clamp d to [-1, 1] to avoid floating point errors outside acos domain
    if d > 1:
        d = 1
    if d < -1:
        d = -1
    e = acos(d)
    f = e / M_PI * 180.0
    return f

@cython.boundscheck(False)  # deactivate bnds checking
cpdef calculate_angle_between_point(dLongitude_degree1_in, dLatitude_degree1_in, dLongitude_degree2_in, dLatitude_degree2_in, dLongitude_degree3_in, dLatitude_degree3_in):
    #all in degree

    cdef double angle3deg
    cdef double x1, y1, z1
    cdef double x2, y2, z2
    cdef double x3, y3, z3
    cdef double x4, y4, z4
    cdef double x5, y5, z5

    # The points in 3D space
    x1, y1, z1 = convert_longitude_latitude_to_sphere_3d(dLongitude_degree1_in, dLatitude_degree1_in)
    x2, y2, z2 = convert_longitude_latitude_to_sphere_3d(dLongitude_degree2_in, dLatitude_degree2_in)
    x3, y3, z3 = convert_longitude_latitude_to_sphere_3d(dLongitude_degree3_in, dLatitude_degree3_in)
    # Vectors in 3D space

    x4 = x1 - x2
    y4 = y1 - y2
    z4 = z1 - z2
    #c3vec[i] = aCoordinate3[i] - aCoordinate2[i]
    x5 = x3 - x2
    y5 = y3 - y2
    z5 = z3 - z2

    angle3deg = calculate_angle_between_vectors_coordinates( x4, y4, z4, x5, y5, z5)
    return  angle3deg

@cython.boundscheck(False)  # Disable bounds checking for performance
@cython.wraparound(False)   # Disable negative indexing for performance
def calculate_distance_to_plane(double dLongitude_degree1_in, double dLatitude_degree1_in,
                                    double dLongitude_degree2_in, double dLatitude_degree2_in,
                                    double dLongitude_degree3_in, double dLatitude_degree3_in):
    """
    Calculate the distance of a point to a plane defined by three points in 3D space.
    """

    cdef double x1, y1, z1
    cdef double x2, y2, z2
    cdef double x3, y3, z3
    cdef double v1_x, v1_y, v1_z
    cdef double v2_x, v2_y, v2_z
    cdef double normal_x, normal_y, normal_z
    cdef double A, B, C, D
    cdef double distance

    # Convert the three points to 3D coordinates
    x1, y1, z1 = convert_longitude_latitude_to_sphere_3d(dLongitude_degree1_in, dLatitude_degree1_in)
    x2, y2, z2 = convert_longitude_latitude_to_sphere_3d(dLongitude_degree2_in, dLatitude_degree2_in)
    x3, y3, z3 = convert_longitude_latitude_to_sphere_3d(dLongitude_degree3_in, dLatitude_degree3_in)

    # Calculate two vectors on the plane
    v1_x = x2 - x1
    v1_y = y2 - y1
    v1_z = z2 - z1

    v2_x = x3 - x1
    v2_y = y3 - y1
    v2_z = z3 - z1

    # Compute the normal vector using the cross product
    normal_x = v1_y * v2_z - v1_z * v2_y
    normal_y = v1_z * v2_x - v1_x * v2_z
    normal_z = v1_x * v2_y - v1_y * v2_x

    # Check if the normal vector is zero (points are collinear)
    if abs(normal_x) < 1e-10 and abs(normal_y) < 1e-10 and abs(normal_z) < 1e-10:
        # Collinear points: return zero, but warn user in documentation
        return 0.0

    # Calculate the plane equation coefficients (A, B, C, D)
    A = normal_x
    B = normal_y
    C = normal_z
    D = -(A * x1 + B * y1 + C * z1)

    # Calculate the distance of the second point to the plane
    distance = abs(A * x2 + B * y2 + C * z2 + D) / sqrt(A**2 + B**2 + C**2)

    return distance

@cython.boundscheck(False)  # deactivate bnds checking
cpdef convert_360_to_180(double dLongitude_in):
    """[This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html]

    Args:
        dLongitude_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    cdef int a
    cdef double dLongitude_out
    a = int(dLongitude_in /180.0)
    dLongitude_out = dLongitude_in - a * 360.0
    return dLongitude_out