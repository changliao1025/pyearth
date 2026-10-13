"""Write Well-Known Text (WKT) geometries to vector GIS files."""

import os
from pathlib import Path
from typing import Any, Optional, Union

from osgeo import ogr, osr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename

ogr.UseExceptions()
osr.UseExceptions()

PathLike = Union[str, os.PathLike]

__all__ = ["gdal_write_wkt_to_vector_file"]


def _resolve_spatial_reference(
    srs: Optional[Union[osr.SpatialReference, int, str]]
) -> Optional[osr.SpatialReference]:
    """Resolve an EPSG code, WKT string, or SpatialReference instance."""
    if srs is None:
        return None
    if isinstance(srs, osr.SpatialReference):
        return srs
    spatial_ref = osr.SpatialReference()
    if isinstance(srs, int):
        spatial_ref.ImportFromEPSG(srs)
        return spatial_ref
    if isinstance(srs, str):
        try:
            spatial_ref.ImportFromWkt(srs)
            return spatial_ref
        except Exception:
            spatial_ref.SetFromUserInput(srs)
            return spatial_ref
    return None


def gdal_write_wkt_to_vector_file(
    wkt: str,
    sFilename_out: PathLike,
    layer_name: str = "geometry",
    spatial_reference: Optional[Union[osr.SpatialReference, int, str]] = None,
    overwrite: bool = True,
) -> None:
    """Save a WKT geometry to a vector file.

    Parameters
    ----------
    wkt : str
        Well-Known Text representation of the geometry.
    sFilename_out : PathLike
        Output file path. Format is determined from the extension (GeoJSON, Shapefile, etc.).
    layer_name : str, optional
        Name for the output layer. Default is "geometry".
    spatial_reference : osr.SpatialReference, int, or str, optional
        Spatial reference system. Can be an osr.SpatialReference object, EPSG code (int),
        or WKT / user input string.
    overwrite : bool, optional
        If True, overwrite existing file. Default is True.

    Raises
    ------
    ValueError
        If WKT is None or empty, or if geometry creation fails.
    FileExistsError
        If the output file exists and overwrite is False.
    RuntimeError
        If dataset, layer, or feature creation fails.
    """
    if not wkt or not wkt.strip():
        raise ValueError("WKT geometry string cannot be None or empty.")

    filename_str = os.fspath(sFilename_out)
    parent_dir = os.path.dirname(os.path.abspath(filename_str))
    if parent_dir:
        os.makedirs(parent_dir, exist_ok=True)

    driver = get_vector_driver_from_filename(filename_str)

    if os.path.exists(filename_str):
        if not overwrite:
            raise FileExistsError(f"Output file already exists: {filename_str}")
        try:
            driver.DeleteDataSource(filename_str)
        except Exception:
            if os.path.isfile(filename_str):
                os.remove(filename_str)

    # Create geometry from WKT
    geometry = ogr.CreateGeometryFromWkt(wkt)
    if geometry is None:
        raise ValueError(f"Could not create geometry from WKT: {wkt[:100]}...")

    geom_type = geometry.GetGeometryType()
    srs = _resolve_spatial_reference(spatial_reference)

    dataset = None
    layer = None
    feature = None

    try:
        dataset = driver.CreateDataSource(filename_str)
        if dataset is None:
            raise RuntimeError(f"Could not create output dataset: {filename_str}")

        layer = dataset.CreateLayer(layer_name, srs=srs, geom_type=geom_type)
        if layer is None:
            raise RuntimeError(
                f"Could not create layer '{layer_name}' in output dataset: {filename_str}"
            )

        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(geometry)

        result = layer.CreateFeature(feature)
        if result != ogr.OGRERR_NONE:
            raise RuntimeError(
                f"Could not create feature in output layer: {filename_str}"
            )

        dataset.FlushCache()
    finally:
        feature = None
        layer = None
        dataset = None
        geometry = None
