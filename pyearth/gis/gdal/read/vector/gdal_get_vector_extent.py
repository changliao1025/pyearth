"""Vector dataset extent calculation utilities using OGR."""

import os
from typing import Optional, Tuple, Union
from osgeo import ogr

PathLike = Union[str, os.PathLike]

__all__ = ["gdal_get_vector_extent"]


def gdal_get_vector_extent(
    filename: Optional[PathLike] = None,
    return_union_geometry: bool = False,
    **kwargs,
) -> Union[
    Tuple[float, float, float, float],
    Tuple[Tuple[float, float, float, float], Optional[ogr.Geometry]],
]:
    """Get the spatial extent of a vector file, optionally with a combined geometry.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the input vector file.
    return_union_geometry : bool, default=False
        If True, return a tuple containing both the extent and a combined geometry
        of all features.

    Returns
    -------
    tuple
        If return_union_geometry is False:
            (min_x, max_x, min_y, max_y)
        If return_union_geometry is True:
            ((min_x, max_x, min_y, max_y), union_geometry)

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    RuntimeError
        If the vector file cannot be opened.
    ValueError
        If the layer cannot be accessed.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")
    if "iFlag_return_union_geometry" in kwargs:
        return_union_geometry = bool(kwargs["iFlag_return_union_geometry"])

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

        extent = layer.GetExtent()
        if extent is None:
            raise ValueError(f"Unable to retrieve extent from {filepath}.")

        if not return_union_geometry:
            return extent

        layer_geom_type = layer.GetGeomType()
        base_geom_type = layer_geom_type & 0xFF

        if base_geom_type == ogr.wkbPoint:
            union_geom = ogr.Geometry(ogr.wkbMultiPoint)
        elif base_geom_type == ogr.wkbLineString:
            union_geom = ogr.Geometry(ogr.wkbMultiLineString)
        elif base_geom_type == ogr.wkbPolygon:
            union_geom = ogr.Geometry(ogr.wkbMultiPolygon)
        else:
            return extent, None

        layer.ResetReading()
        for feature in layer:
            if feature is None:
                continue

            geometry = feature.GetGeometryRef()
            if geometry is None:
                continue

            geom_clone = geometry.Clone()
            geom_clone.FlattenTo2D()
            flat_geom_type = geom_clone.GetGeometryType()

            if flat_geom_type in (
                ogr.wkbMultiPoint,
                ogr.wkbMultiLineString,
                ogr.wkbMultiPolygon,
            ):
                for i in range(geom_clone.GetGeometryCount()):
                    part = geom_clone.GetGeometryRef(i)
                    if part is not None:
                        union_geom.AddGeometry(part)
            elif flat_geom_type in (ogr.wkbPoint, ogr.wkbLineString, ogr.wkbPolygon):
                union_geom.AddGeometry(geom_clone)

        return extent, union_geom

    finally:
        layer = None
        dataset = None
