"""Raster spatial reference extraction utilities using GDAL."""

import os
from typing import Optional, Union
from osgeo import gdal, osr

PathLike = Union[str, os.PathLike]

__all__ = [
    "gdal_get_raster_spatial_reference_wkt",
    "gdal_get_raster_spatial_reference",
]


def gdal_get_raster_spatial_reference_wkt(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> str:
    """Return the spatial reference of a raster as a WKT string.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the raster file readable by GDAL.

    Returns
    -------
    str
        The spatial reference in Well Known Text (WKT) format.

    Raises
    ------
    FileNotFoundError
        If the raster path does not exist.
    RuntimeError
        If GDAL cannot open the raster.
    ValueError
        If the raster lacks spatial reference information.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")

    if filename is None:
        raise ValueError("Filename must be provided.")

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist.")

    dataset = gdal.Open(filepath, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError(f"Unable to open raster {filepath}.")

    try:
        spatial_ref = getattr(dataset, "GetSpatialRef", lambda: None)()
        if spatial_ref is not None:
            try:
                wkt = spatial_ref.ExportToWkt()
                if wkt:
                    return wkt
            except Exception:
                pass

        wkt_fallback = dataset.GetProjectionRef() or dataset.GetProjection()
        if wkt_fallback:
            return wkt_fallback

        raise ValueError(
            f"Raster {filepath} does not define spatial reference metadata."
        )
    finally:
        dataset = None


def gdal_get_raster_spatial_reference(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> osr.SpatialReference:
    """Return the spatial reference of a raster as an osr.SpatialReference object.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the raster file readable by GDAL.

    Returns
    -------
    osr.SpatialReference
        Spatial reference object with traditional GIS axis mapping order.
    """
    wkt = gdal_get_raster_spatial_reference_wkt(filename, **kwargs)
    sr = osr.SpatialReference()
    sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    sr.ImportFromWkt(wkt)
    return sr
