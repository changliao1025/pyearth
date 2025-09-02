import os
from osgeo import ogr
def copy_geometery_without_attributes(sFilename_geojson_in, sFilename_geojson_out):
    """
    Copy a GeoJSON file keeping only the geometry, removing all attributes.

    Args:
        sFilename_geojson_in (str): Input GeoJSON file path
        sFilename_geojson_out (str, optional): Output GeoJSON file path.
                                              If None, creates a new file with '_geometry_only' suffix

    Returns:
        str: Path to the output GeoJSON file
    """

    # Check if input file exists
    if not os.path.exists(sFilename_geojson_in):
        raise FileNotFoundError(f"Input file not found: {sFilename_geojson_in}")


    if os.path.exists(sFilename_geojson_out):
        os.remove(sFilename_geojson_out)

    # Open input dataset
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset_in = pDriver.Open(sFilename_geojson_in, 0)  # 0 = read-only

    if pDataset_in is None:
        raise ValueError(f"Could not open input file: {sFilename_geojson_in}")

    pLayer_in = pDataset_in.GetLayer(0)

    # Get spatial reference from input layer
    pSRS = pLayer_in.GetSpatialRef()

    # Create output dataset
    if os.path.exists(sFilename_geojson_out):
        os.remove(sFilename_geojson_out)

    pDataset_out = pDriver.CreateDataSource(sFilename_geojson_out)

    if pDataset_out is None:
        raise ValueError(f"Could not create output file: {sFilename_geojson_out}")

    # Get geometry type from first feature
    pLayer_in.ResetReading()
    pFeature_first = pLayer_in.GetNextFeature()
    if pFeature_first is None:
        raise ValueError("Input file contains no features")

    geometry_type = pFeature_first.GetGeometryRef().GetGeometryType()

    # Create output layer with only geometry (no attributes)
    pLayer_out = pDataset_out.CreateLayer("geometry_only", pSRS, geometry_type)

    if pLayer_out is None:
        raise ValueError("Could not create output layer")

    # Get layer definition for creating new features
    pLayerDefn_out = pLayer_out.GetLayerDefn()

    # Reset reading and copy geometries
    pLayer_in.ResetReading()

    for pFeature_in in pLayer_in:
        # Get geometry from input feature
        pGeometry = pFeature_in.GetGeometryRef()

        if pGeometry is not None:
            # Create new feature with only geometry
            pFeature_out = ogr.Feature(pLayerDefn_out)
            pFeature_out.SetGeometry(pGeometry.Clone())

            # Add feature to output layer
            pLayer_out.CreateFeature(pFeature_out)

            # Clean up
            pFeature_out = None

    # Clean up
    pDataset_in = None
    pDataset_out = None

    print(f"Geometry-only GeoJSON created: {sFilename_geojson_out}")
    return sFilename_geojson_out