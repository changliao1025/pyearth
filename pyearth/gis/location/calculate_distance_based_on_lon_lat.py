from math import radians, cos, sin, asin, sqrt

def calculate_distance_based_on_lon_lat(dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    https://en.wikipedia.org/wiki/Great-circle_distance
    The current algorothm is the  haversine formula 

    Args:
        dLongitude_from (float): The longitude of the start point
        dLatitude_from (float):  The latitude of the start point
        dLongitude_to (float):  The longitude of the end point
        dLatitude_to (float):  The latitude of the end point

    Returns:
        float: The great circle distance
    """
    # convert decimal degrees to radians 
    dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to = map(radians, [dLongitude_from, dLatitude_from, dLongitude_to, dLatitude_to])

    # haversine formula 
    dlon = dLongitude_to - dLongitude_from 
    dlat = dLatitude_to - dLatitude_from 
    a = sin(dlat/2)**2 + cos(dLatitude_from) * cos(dLatitude_to) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6378137.0# Radius of earth in kilometers. Use 3956 for miles
    return c * r