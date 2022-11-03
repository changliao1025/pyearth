import numpy as np
def degree_to_meter(dLatitude_mean_in, dResolution_degree_in):
    """[summary]

    Args:
        dLatitude_in (float): The latitide of interest
        dResolution_degree_in (float): The desired resolution in degree

    Returns:
        float: The resolution in meter at this latitude
    """
    dRadius = 6378137.0 
    dRadius2 = dRadius * np.cos( dLatitude_mean_in / 180.0 * np.pi)
    dResolution_meter = dResolution_degree_in / 360.0 * (2*np.pi * dRadius2)

    return dResolution_meter


def meter_to_degree(dLatitude_mean_in, dResolution_meter_in ):
    """[summary]

    Args:
        dResolution_meter_in (float): [description]
        dLatitude_mean_in (float): [description]

    Returns:
        float: The desired resoluton in degree
    """
    dLatitude_mean_in = abs(dLatitude_mean_in)

    dRadius = 6378137.0
    dRadius2 = dRadius * np.cos( dLatitude_mean_in / 180.0 * np.pi)

 
    dResolution_degree= dResolution_meter_in/(2*np.pi * dRadius2) * 360.0

    return dResolution_degree