import numpy as np
from pyearth.system.define_global_variables import *  

from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import calculate_distance_based_on_longitude_latitude



def calculate_spherical_triangle_area(aLongitutde_in, 
                                      aLatitude_in, 
                                      iFlag_radian = None, 
                                      dRadius_in = None):
    """
    Compute the area of a spherical triangle.

    Parameters:
        A, B, C: Tuple or list representing the spherical coordinates (longitude, latitude) of three points.

    Returns:
        The area of the spherical triangle.
    """
    
    if iFlag_radian is None:# convert decimal degrees to radians
        aLongitutde_radian_in = np.radian(aLongitutde_radian_in)
        aLatitude_radian_in = np.radian(aLatitude_in)
    else:
        aLongitutde_radian_in = aLongitutde_in
        aLatitude_radian_in = aLatitude_in
        pass        
    
    dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to = aLongitutde_radian_in[0], aLatitude_radian_in[0], aLongitutde_radian_in[1], aLatitude_radian_in[1]

    a = calculate_distance_based_on_longitude_latitude(dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to, iFlag_radian=1)

    dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to = aLongitutde_radian_in[1], aLatitude_radian_in[1], aLongitutde_radian_in[2], aLatitude_radian_in[2] 

    b = calculate_distance_based_on_longitude_latitude(dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to, iFlag_radian=1)

    dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to = aLongitutde_radian_in[2], aLatitude_radian_in[2], aLongitutde_radian_in[0], aLatitude_radian_in[0] 

    c = calculate_distance_based_on_longitude_latitude(dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to, iFlag_radian=1)

    s = 0.5 * (a + b + c)

    tanqe = np.sqrt(np.tan(0.5 * s) * np.tan(0.5 * (s - a)) * np.tan(0.5 * (s - b)) * np.tan(0.5 * (s - c)))

    e = 4.0 * np.arctan(tanqe)

    if iFlag_radian is not None: #requested in radian
        return e
    else:       
        if dRadius_in is not None:
            f = e * dRadius_in**2
        else:
            f = e * earth_radius**2
        
        return f

    