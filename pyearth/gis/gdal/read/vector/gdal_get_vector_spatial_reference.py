from osgeo import ogr

def gdal_get_vector_spatial_reference_wkt(sFilename_in):
    dataset = ogr.Open(sFilename_in)
    layer = dataset.GetLayer(0)
    spatial_ref = layer.GetSpatialRef()
    pProjection_reference = spatial_ref.ExportToWkt()
    return pProjection_reference

