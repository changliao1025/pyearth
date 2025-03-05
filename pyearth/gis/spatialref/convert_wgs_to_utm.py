# convert_wgs_to_utm function, see https://stackoverflow.com/a/40140326/4556479
import math


def get_utm_zone(lon: float, lat: float) -> int:
    """Return UTM zone based on longitude and latitude"""
    return int((lon + 180) / 6) + 1

def get_utm_epsg_code(lon: float, lat: float):
    """Based on lat and lng, return best utm epsg-code"""
    iUTM_zone = get_utm_zone(lon, lat)
    utm_band = '{:02d}'.format(iUTM_zone)
    if lat >= 0:
        epsg_code = '326' + utm_band
    else:
        epsg_code = '327' + utm_band
    #convert back to int
    epsg_code = int(epsg_code)
    return epsg_code
