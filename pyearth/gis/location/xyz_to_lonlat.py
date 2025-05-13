import numpy as np
def xyz_to_lonlat(x, y, z):
    # Normalize the coordinates
    norm = np.sqrt(x**2 + y**2 + z**2)
    x /= norm
    y /= norm
    z /= norm

    # Convert to spherical coordinates
    lon = np.arctan2(y, x)
    lat = np.arcsin(z)

    # Convert radians to degrees
    lon = np.degrees(lon)
    lat = np.degrees(lat)

    return float(lon), float(lat)