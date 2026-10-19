"""Vector boundary geometry extraction utilities using OGR."""

import os
from typing import Optional, Tuple, Union

from osgeo import ogr

PathLike = Union[str, os.PathLike]

__all__ = ["gdal_get_vector_boundary"]


def gdal_get_vector_boundary(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> Tuple[str, Tuple[float, float, float, float]]:
    """Extract the boundary geometry and envelope from a vector file.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the vector file containing polygon or multipolygon boundary geometry.

    Returns
    -------
    tuple
        (boundary_wkt, (min_x, max_x, min_y, max_y))

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    RuntimeError
        If OGR cannot open the file.
    ValueError
        If no valid polygon geometry is found.
    """
    if filename is None:
        filename = kwargs.get("sFilename_boundary_in")

    if filename is None:
        raise ValueError("Filename must be provided.")

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File does not exist: {filepath}")

    dataset = None
    layer = None
    try:
        dataset = ogr.Open(filepath, 0)
        if dataset is None:
            raise RuntimeError(f"Unable to open vector file {filepath}.")

        layer = dataset.GetLayer(0)
        if layer is None:
            raise ValueError(f"Unable to access layer from {filepath}.")

        feature_count = layer.GetFeatureCount()
        if feature_count == 0:
            raise ValueError(f"No features found in {filepath}.")

        boundary_geom = None
        layer.ResetReading()
        for feature in layer:
            if feature is None:
                continue
            geometry = feature.GetGeometryRef()
            if geometry is None:
                continue

            geom_name = geometry.GetGeometryName().upper()
            if geom_name == "POLYGON":
                boundary_geom = geometry.Clone()
                break
            elif geom_name == "MULTIPOLYGON":
                boundary_geom = geometry.Clone()
                break

        if boundary_geom is None:
            raise ValueError(f"No valid polygon geometry found in {filepath}.")

        boundary_wkt = boundary_geom.ExportToWkt()
        extent = boundary_geom.GetEnvelope()

        return boundary_wkt, extent

    finally:
        layer = None
        dataset = None
