import numpy as np
from pyearth.system.define_global_variables import *

def calculate_distance_based_on_longitude_latitude(dLongitude_from,
                                                   dLatitude_from,
                                                   dLongitude_to,
                                                   dLatitude_to,
                                                   iFlag_radian = None,
                                                   dRadius_in = None):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    https://en.wikipedia.org/wiki/Great-circle_distance
    The current algorothm is the  haversine formula 

    Args:
        dLongitude_radian_from (float): The longitude of the start point
        dLatitude_radian_from (float):  The latitude of the start point
        dLongitude_radian_to (float):  The longitude of the end point
        dLatitude_radian_to (float):  The latitude of the end point

    Returns:
        float: The great circle distance
    """
    
    if iFlag_radian is None:# convert decimal degrees to radians
        dLongitude_radian_from, dLatitude_radian_from, dLongitude_radian_to, dLatitude_radian_to = np.radians([dLongitude_from,
                                                                                                          dLatitude_from,
                                                                                                          dLongitude_to,
                                                                                                          dLatitude_to])
    else: #already in radian
        dLongitude_radian_from, dLatitude_radian_from, dLongitude_radian_to, dLatitude_radian_to = dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to
        pass

    # haversine formula
    dlon = dLongitude_radian_to - dLongitude_radian_from
    dlat = dLatitude_radian_to - dLatitude_radian_from

    a = np.sin(dlat/2)**2 + np.cos(dLatitude_radian_from) * np.cos(dLatitude_radian_to) * np.sin(dlon/2)**2
    
    
    c = 2 * np.arcsin(np.sqrt(a))

    if iFlag_radian is not None:
        return c
            
    else:  
        if dRadius_in is not None:
            d = c * dRadius_in
        else:
            d = c * earth_radius
        
        return d    

    
