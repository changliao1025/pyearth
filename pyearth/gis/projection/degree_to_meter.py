import numpy as np
def degree_to_meter(dLatitude, dResolution_degree):
    """[summary]

    Args:
        dLatitude ([type]): [description]
        dResolution_degree ([type]): [description]

    Returns:
        [type]: [description]
    """
    dRadius = 6378137.0 
    dRadius2 = dRadius * np.cos( dLatitude / 180.0 * np.pi)
    dResolution_meter = dResolution_degree / 360.0 * (2*np.pi * dRadius2)

    return dResolution_meter