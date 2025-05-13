from osgeo import  osr
from pyearth.gis.spatialref.convert_wgs_to_utm import get_utm_epsg_code
def get_utm_spatial_reference_wkt(dLongitude_in: float, dLatitude_in: float) -> str:
    dLongitude_in = float(dLongitude_in)
    dLatitude_in = float(dLatitude_in)
    if -180 <= dLongitude_in <= 180:
        epsg_code = get_utm_epsg_code(dLongitude_in, dLatitude_in)
        utm_sr = osr.SpatialReference()
        utm_sr.ImportFromEPSG(epsg_code)
        utm_projection = utm_sr.ExportToWkt()
        return utm_projection
    else:
        raise ValueError("Longitude must be in the range [-180, 180].")