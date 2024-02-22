import numpy as np
def convert_longitude_latitude_to_sphere_3d(dLongitude_in, dLatitude_in, iFlag_radian= None):
    """Convert a point given latitude and longitude in radians to
    3-dimensional space, assuming a sphere radius of one."""

    if iFlag_radian is None:
        dLongitude_radian, dLatitude_radian = np.radians([dLongitude_in, dLatitude_in])
    else:
        dLongitude_radian, dLatitude_radian = dLongitude_in, dLatitude_in
        pass

    a = np.cos(dLatitude_radian) * np.cos(dLongitude_radian)
    b = np.cos(dLatitude_radian) * np.sin(dLongitude_radian)
    c = np.sin(dLatitude_radian)
    d = np.array((a,b,c))
    
    return d