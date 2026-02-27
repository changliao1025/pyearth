
import os
from osgeo import ogr, osr
from pyearth.gis.geometry.calculate_polygon_area  import calculate_polygon_area
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename

def save_points_as_polygon(aPoint_in, sFilename_out):
    pDriver = get_vector_driver_from_filename(sFilename_out)
    if os.path.exists(sFilename_out):
        os.remove(sFilename_out)

    pDataset = pDriver.CreateDataSource(sFilename_out)
    pSrcSpatialRef = osr.SpatialReference()
    pSrcSpatialRef.ImportFromEPSG(4326)
    pLayer = pDataset.CreateLayer('buffer_polygon', pSrcSpatialRef, geom_type=ogr.wkbPolygon)
    aLon = list()
    aLat = list()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for pVertex in aPoint_in:
        ring.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)
        aLon.append(pVertex.dLongitude_degree)
        aLat.append(pVertex.dLatitude_degree)

    ring.CloseRings()
    pArea_field = ogr.FieldDefn('area', ogr.OFTReal)
    pArea_field.SetWidth(20)
    pArea_field.SetPrecision(2)
    pLayer.CreateField(pArea_field)
    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)
    # close the ring geometry before creating the polygon
    dArea = calculate_polygon_area(aLon, aLat)
    pPolygon = ogr.Geometry(ogr.wkbPolygon)
    pPolygon.AddGeometry(ring)
    pFeature.SetGeometry(pPolygon)
    pFeature.SetField("area", dArea )
    pLayer.CreateFeature(pFeature)
    #flush the cache and close the file
    pDataset.FlushCache()
    pDataset.Destroy()
    pSrcSpatialRef = None

    return