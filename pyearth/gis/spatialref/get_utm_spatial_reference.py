from osgeo import ogr, osr, gdal
def get_utm_spatial_reference(dLongitude_in):
    if -180 <= dLongitude_in <= 180:
        zone = int((dLongitude_in + 180) / 6) + 1
        hemisphere = 'N' if dLongitude_in >= 0 else 'S'
        epsg_code = 32600 + zone if hemisphere == 'N' else 32700 + zone
        utm_sr = osr.SpatialReference()
        utm_sr.ImportFromEPSG(epsg_code)
        return utm_sr
    else:
        raise ValueError("Longitude must be in the range [-180, 180].")