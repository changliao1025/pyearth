import os
from typing import Tuple
from osgeo import ogr, gdal
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename


def gdal_get_vector_boundary(
    sFilename_boundary_in: str,
) -> Tuple[str, Tuple[float, float, float, float]]:
    """Extract the boundary geometry and extent from a vector file.

    This function reads a vector file and extracts the boundary geometry
    (polygon or multipolygon) along with its spatial extent.

    Parameters
    ----------
    sFilename_boundary_in : str
        Path to the vector file containing boundary geometry.
        Supports GeoJSON, Shapefile, GeoPackage, KML, and other OGR formats.

    Returns
    -------
    tuple
        A tuple containing:
        - boundary_wkt (str): Boundary geometry in Well Known Text (WKT) format
        - extent (tuple): Spatial extent as (min_x, max_x, min_y, max_y)

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If the file extension is not supported.
    RuntimeError
        If OGR cannot open the file.
    ValueError
        If the layer cannot be accessed, has no features, or contains
        unsupported geometry types.

    Notes
    -----
    - Only processes POLYGON and MULTIPOLYGON geometries
    - Returns the geometry from the first feature in the layer
    - For MULTIPOLYGON, creates a new MultiPolygon containing all sub-geometries
    - Automatically determines the appropriate OGR driver from file extension
    """

    if not os.path.isfile(sFilename_boundary_in):
        raise FileNotFoundError(f"File does not exist: {sFilename_boundary_in}")

    dataset = None
    try:
        # Get the appropriate driver based on file extension
        driver = get_vector_driver_from_filename(sFilename_boundary_in)

        dataset = driver.Open(sFilename_boundary_in, gdal.GA_ReadOnly)
        if dataset is None:
            raise RuntimeError(f"Unable to open vector file {sFilename_boundary_in}.")

        layer = dataset.GetLayer(0)
        if layer is None:
            raise ValueError(f"Unable to access layer from {sFilename_boundary_in}.")

        feature_count = layer.GetFeatureCount()
        if feature_count == 0:
            raise ValueError(f"No features found in {sFilename_boundary_in}.")

        # Process the first feature
        boundary_geom = None
        for feature in layer:
            geometry = feature.GetGeometryRef()
            if geometry is None:
                continue

            geometry_type = geometry.GetGeometryName()

            if geometry_type == "POLYGON":
                boundary_geom = geometry.Clone()
                break

            elif geometry_type == "MULTIPOLYGON":
                # Create a new MultiPolygon with all sub-geometries
                boundary_geom = ogr.Geometry(ogr.wkbMultiPolygon)
                num_geometries = geometry.GetGeometryCount()
                for i in range(num_geometries):
                    sub_geometry = geometry.GetGeometryRef(i)
                    if sub_geometry is not None:
                        boundary_geom.AddGeometry(sub_geometry)
                break

            else:
                raise ValueError(
                    f"Unsupported geometry type '{geometry_type}'. "
                    f"Expected POLYGON or MULTIPOLYGON."
                )

        if boundary_geom is None:
            raise ValueError(
                f"No valid polygon geometry found in {sFilename_boundary_in}."
            )

        boundary_wkt = boundary_geom.ExportToWkt()
        if not boundary_wkt:
            raise ValueError("Failed to export boundary geometry to WKT.")

        extent = boundary_geom.GetEnvelope()
        if extent is None:
            raise ValueError("Failed to compute geometry envelope.")

        return boundary_wkt, extent

    finally:
        dataset = None
