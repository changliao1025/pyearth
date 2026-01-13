from osgeo import osr


def is_wgs84_projection(projection_wkt: str) -> bool:
    _SPATIAL_REF_WGS84 = osr.SpatialReference()
    _SPATIAL_REF_WGS84.ImportFromEPSG(4326)
    """Return True if the provided projection WKT matches EPSG:4326."""
    if not projection_wkt:
        return False

    spatial_ref = osr.SpatialReference()
    if spatial_ref.ImportFromWkt(projection_wkt) != 0:
        return False

    return spatial_ref.IsSame(_SPATIAL_REF_WGS84) == 1
