import os
from typing import Optional

from osgeo import ogr
from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_format_from_filename,
    get_vector_driver_from_filename,
)


def gdal_write_wkt_to_vector_file(
    wkt: str, sFilename_out: str, layer_name: str = "geometry", overwrite: bool = True
) -> None:
    """Save a WKT geometry to a vector file.

    Parameters
    ----------
    wkt : str
        Well-Known Text representation of the geometry.
    sFilename_out : str
        Output file path. The file format is determined from the extension.
        Supported formats include GeoJSON, Shapefile, GeoPackage, KML, etc.
    layer_name : str, optional
        Name for the output layer. Default is "geometry".
    overwrite : bool, optional
        If True, overwrite existing file. Default is True.

    Raises
    ------
    ValueError
        If WKT is None or empty, if the file extension is not supported,
        or if geometry creation from WKT fails.
    RuntimeError
        If dataset, layer, or feature creation fails.
    FileExistsError
        If the output file exists and overwrite is False.

    Notes
    -----
    - Automatically determines the output format from file extension
    - Creates a single-feature layer with the provided geometry
    - Geometry type is inferred from the WKT string
    """

    if not wkt:
        raise ValueError("WKT geometry string cannot be None or empty.")

    # Check if file exists and handle overwrite
    if os.path.exists(sFilename_out):
        if not overwrite:
            raise FileExistsError(f"Output file already exists: {sFilename_out}")
        os.remove(sFilename_out)

    # Get the appropriate driver based on file extension
    driver = get_vector_driver_from_filename(sFilename_out)

    # Create geometry from WKT first to infer geometry type
    geometry = ogr.CreateGeometryFromWkt(wkt)
    if geometry is None:
        raise ValueError(f"Could not create geometry from WKT. Invalid WKT string.")

    geom_type = geometry.GetGeometryType()

    dataset = None
    layer = None
    feature = None

    try:
        # Create output dataset
        dataset = driver.CreateDataSource(sFilename_out)
        if dataset is None:
            raise RuntimeError(f"Could not create output dataset: {sFilename_out}")

        # Create layer with appropriate geometry type
        layer = dataset.CreateLayer(layer_name, geom_type=geom_type)
        if layer is None:
            raise RuntimeError(
                f"Could not create layer in output dataset: {sFilename_out}"
            )

        # Create and configure feature
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(geometry)

        # Write feature to layer
        result = layer.CreateFeature(feature)
        if result != ogr.OGRERR_NONE:
            raise RuntimeError(
                f"Could not create feature in output layer: {sFilename_out}"
            )

    finally:
        # Clean up in reverse order
        feature = None
        layer = None
        dataset = None
        geometry = None
