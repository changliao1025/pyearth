"""Vector spatial reference extraction utilities using OGR."""

import os
from typing import Optional, Union

from osgeo import ogr, osr

PathLike = Union[str, os.PathLike]

__all__ = [
    "gdal_get_vector_spatial_reference_wkt",
    "gdal_get_vector_spatial_reference",
]


def gdal_get_vector_spatial_reference_wkt(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> str:
    """Return the spatial reference of a vector file as a WKT string.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the vector file.

    Returns
    -------
    str
        The spatial reference in Well Known Text (WKT) format.

    Raises
    ------
    FileNotFoundError
        If the vector file does not exist.
    RuntimeError
        If OGR cannot open the vector file.
    ValueError
        If the layer cannot be accessed or lacks spatial reference metadata.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")

    if filename is None:
        raise ValueError("Filename must be provided.")

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist.")

    dataset = None
    layer = None
    try:
        dataset = ogr.Open(filepath, 0)  # 0 = read-only
        if dataset is None:
            raise RuntimeError(f"Unable to open vector file {filepath}.")

        layer = dataset.GetLayer(0)
        if layer is None:
            raise ValueError(f"Unable to access layer from {filepath}.")

        spatial_ref = layer.GetSpatialRef()
        if spatial_ref is None:
            raise ValueError(
                f"Vector file {filepath} does not define spatial reference metadata."
            )

        projection_wkt = spatial_ref.ExportToWkt()
        if not projection_wkt:
            raise ValueError(
                f"Vector file {filepath} contains empty spatial reference metadata."
            )

        return projection_wkt

    finally:
        layer = None
        dataset = None


def gdal_get_vector_spatial_reference(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> osr.SpatialReference:
    """Return the spatial reference of a vector file as an osr.SpatialReference object.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the vector file.

    Returns
    -------
    osr.SpatialReference
        Spatial reference object with traditional GIS axis mapping order.
    """
    wkt = gdal_get_vector_spatial_reference_wkt(filename, **kwargs)
    sr = osr.SpatialReference()
    sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    sr.ImportFromWkt(wkt)
    return sr
