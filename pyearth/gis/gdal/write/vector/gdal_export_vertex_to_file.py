import os
from osgeo import ogr, osr
from pyearth.gis.geometry.calculate_polygon_area  import calculate_polygon_area
from pyearth.gis.gdal.gdal_vector_format_support import get_format_from_extension

def export_vertex_to_vector(aVertex_in,
        sFilename_vector_in,
        iFlag_projected_in=None,
        pSpatial_reference_in=None,
        aAttribute_data=None):
    """
    Export vertices to any supported vector format.

    Args:
        aVertex_in (list): List of vertex objects with coordinate information
        sFilename_vector_in (str): Output filename with extension
        iFlag_projected_in (int, optional): Flag for coordinate system (0=geographic, 1=projected). Defaults to 0.
        pSpatial_reference_in (osr.SpatialReference, optional): Spatial reference system. Defaults to WGS84.
        aAttribute_data (list, optional): Optional attribute data for each vertex. Defaults to None.
    """

    if os.path.exists(sFilename_vector_in):
        os.remove(sFilename_vector_in)

    iFlag_projected_in = 0 if iFlag_projected_in is None else 1

    if  pSpatial_reference_in is None:
        pSpatial_reference_in = osr.SpatialReference()
        pSpatial_reference_in.ImportFromEPSG(4326)    # WGS84 lat/lon

    # Determine format from file extension
    sFormat = get_format_from_extension(sFilename_vector_in)

    iFlag_attribute = 1 if aAttribute_data is not None else 0
    aAttribute = aAttribute_data if aAttribute_data is not None else []

    #nVertex = len(aVertex_in)
    pDriver = ogr.GetDriverByName(sFormat)
    pDataset = pDriver.CreateDataSource(sFilename_vector_in)
    pLayer = pDataset.CreateLayer('vertex', pSpatial_reference_in, ogr.wkbPoint)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('pointid', ogr.OFTInteger64)) #long type for high resolution
    if iFlag_attribute ==1:
        pLayer.CreateField(ogr.FieldDefn('connectivity', ogr.OFTInteger64)) #long type for high resolution
        pass

    pLayerDefn = pLayer.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayerDefn)


    for lID, pVertex in enumerate(aVertex_in):
        pPoint = ogr.Geometry(ogr.wkbPoint)
        if iFlag_projected_in == 1:
            pPoint.AddPoint(pVertex.dx, pVertex.dy)
        else:
            pPoint.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)

        pGeometry_out = ogr.CreateGeometryFromWkb(pPoint.ExportToWkb())
        pFeature_out.SetGeometry(pGeometry_out)
        pFeature_out.SetField("pointid", lID + 1)

        if iFlag_attribute == 1:
            pFeature_out.SetField("connectivity", int(aAttribute[lID]))

        pLayer.CreateFeature(pFeature_out)

    pDataset.FlushCache()
    pDataset = pLayer = pFeature_out  = None

    return


def export_vertex_as_polygon(aVertex_in, sFilename_out, pSpatial_reference_in=None):
    """
    Export vertices as a polygon to any supported vector format.

    Args:
        aVertex_in (list): List of vertex objects with coordinate information
        sFilename_out (str): Output filename with extension
        pSpatial_reference_in (osr.SpatialReference, optional): Spatial reference system. Defaults to WGS84.
    """
    if os.path.exists(sFilename_out):
        os.remove(sFilename_out)

    # Determine format from file extension
    sFormat = get_format_from_extension(sFilename_out)

    pDriver = ogr.GetDriverByName(sFormat)
    pDataset = pDriver.CreateDataSource(sFilename_out)

    if pSpatial_reference_in is None:
        pSrcSpatialRef = osr.SpatialReference()
        pSrcSpatialRef.ImportFromEPSG(4326)
    else:
        pSrcSpatialRef = pSpatial_reference_in

    pLayer = pDataset.CreateLayer('polygon', pSrcSpatialRef, geom_type=ogr.wkbPolygon)
    aLon = list()
    aLat = list()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for pVertex in aVertex_in:
        ring.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)
        aLon.append(pVertex.dLongitude_degree)
        aLat.append(pVertex.dLatitude_degree)

    ring.AddPoint(aVertex_in[0].dLongitude_degree, aVertex_in[0].dLatitude_degree)
    pArea_field = ogr.FieldDefn('area', ogr.OFTReal)
    pArea_field.SetWidth(20)
    pArea_field.SetPrecision(2)
    pLayer.CreateField(pArea_field)
    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)
    #add the first point to close the ring
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



