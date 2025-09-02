import numpy as np
from pyearth.system.define_global_variables import *

def calculate_distance_based_on_longitude_latitude(aLongitude_from,
                                                       aLatitude_from,
                                                       aLongitude_to,
                                                       aLatitude_to,
                                                       iFlag_radian=None,
                                                       dRadius_in=None):
    """
    Calculate the great circle distances between arrays of points
    on the earth (specified in decimal degrees) using the haversine formula
    https://en.wikipedia.org/wiki/Great-circle_distance

    This function is a vectorized version of calculate_distance_based_on_longitude_latitude
    that can operate on arrays of coordinates efficiently.

    Args:
        aLongitude_from (array-like): The longitudes of the start points
        aLatitude_from (array-like): The latitudes of the start points
        aLongitude_to (array-like): The longitudes of the end points
        aLatitude_to (array-like): The latitudes of the end points
        iFlag_radian (int, optional): Flag to indicate if input is in radians. If None, input is assumed to be in degrees.
        dRadius_in (float, optional): Custom earth radius in meters. If None, default earth_radius is used.

    Returns:
        ndarray: The great circle distances in meters (or in radians if iFlag_radian is not None)
    """
    # Convert inputs to numpy arrays to ensure vectorized operations
    aLongitude_from = np.asarray(aLongitude_from)
    aLatitude_from = np.asarray(aLatitude_from)
    aLongitude_to = np.asarray(aLongitude_to)
    aLatitude_to = np.asarray(aLatitude_to)

    if iFlag_radian is None:  # convert decimal degrees to radians
        dLongitude_radian_from = np.radians(aLongitude_from)
        dLatitude_radian_from = np.radians(aLatitude_from)
        dLongitude_radian_to = np.radians(aLongitude_to)
        dLatitude_radian_to = np.radians(aLatitude_to)
    else:  # already in radian
        dLongitude_radian_from = aLongitude_from
        dLatitude_radian_from = aLatitude_from
        dLongitude_radian_to = aLongitude_to
        dLatitude_radian_to = aLatitude_to

    # haversine formula
    dlon = dLongitude_radian_to - dLongitude_radian_from
    dlat = dLatitude_radian_to - dLatitude_radian_from

    a = np.sin(dlat/2)**2 + np.cos(dLatitude_radian_from) * np.cos(dLatitude_radian_to) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    if iFlag_radian is not None:
        return c
    else:
        radius = dRadius_in if dRadius_in is not None else earth_radius
        d = c * radius
        return d