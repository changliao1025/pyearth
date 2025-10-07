import os
from osgeo import ogr
def gdal_write_wkt_to_geojson(wkt, sFilename_out):
    pDriver = ogr.GetDriverByName("GeoJSON")
    """
    Save a WKT geometry to a GeoJSON file.

    Args:
        wkt (str): Well-Known Text representation of the geometry
        sFilename_out (str): Output file path
        sBasin (str): Basin identifier for logging
        pDriver: OGR driver for GeoJSON

    Returns:
        bool: True if successful, False otherwise
    """
    if wkt is None:
        return False

    # Remove existing file if it exists
    if os.path.exists(sFilename_out):
        os.remove(sFilename_out)

    # Create output dataset
    pDataset_out = pDriver.CreateDataSource(sFilename_out)
    if pDataset_out is None:
        print(f"Could not create output dataset: {sFilename_out}")
        return False

    # Create layer
    pLayer_out = pDataset_out.CreateLayer("watershed_boundary", geom_type=ogr.wkbPolygon)
    if pLayer_out is None:
        print(f"Could not create layer in output dataset: {sFilename_out}")
        pDataset_out = None
        return False

    # Create geometry from WKT
    pGeometry_out = ogr.CreateGeometryFromWkt(wkt)
    if pGeometry_out is None:
        print(f"Could not create geometry from WKT: {wkt}")
        pDataset_out = None
        return False

    # Create and save feature
    pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
    pFeature_out.SetGeometry(pGeometry_out)

    if pLayer_out.CreateFeature(pFeature_out) != ogr.OGRERR_NONE:
        print(f"Could not create feature in output layer: {sFilename_out}")
        pFeature_out = None
        pDataset_out = None
        return False

    # Clean up
    pFeature_out = None
    pDataset_out = None


    return True