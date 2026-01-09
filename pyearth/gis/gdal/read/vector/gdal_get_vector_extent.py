import os
from typing import Optional, Tuple, Union

from osgeo import gdal, osr, ogr


def gdal_get_vector_extent(
    sFilename_in: str, iFlag_return_union_geometry: bool = False
) -> Union[
    Tuple[float, float, float, float],
    Tuple[Tuple[float, float, float, float], Optional[ogr.Geometry]],
    None,
]:
    """Get the spatial extent of a vector file, optionally with a union geometry.

    This function can also return a single geometry that is a collection of all
    feature geometries in the source file.

    Parameters
    ----------
    sFilename_in : str
        Path to the input vector file.
    iFlag_return_union_geometry : bool, optional
        If True, return a tuple containing both the extent and a geometry collection
        of all features. Default is False.

    Returns
    -------
    tuple or None
        If iFlag_return_union_geometry is False:
            Returns the extent as (min_x, max_x, min_y, max_y), or None on error.
        If iFlag_return_union_geometry is True:
            Returns (extent, union_geometry) where union_geometry is a single OGR
            Geometry object (e.g., MultiPolygon) containing all geometries from
            the layer, or (None, None) on error.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    RuntimeError
        If the vector file cannot be opened.
    ValueError
        If the layer cannot be accessed or has no features.

    Notes
    -----
    When creating a union geometry, the function:
    - Determines the appropriate multi-geometry type based on the layer's geometry type
    - Flattens all geometries to 2D
    - Handles multi-geometries by adding their individual parts
    - Returns None for the union if the geometry type is not supported (Point, LineString, or Polygon)
    """

    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"File {sFilename_in} does not exist.")

    dataset = None
    try:
        dataset = ogr.Open(sFilename_in, 0)  # 0 = read-only
        if dataset is None:
            raise RuntimeError(f"Unable to open vector file {sFilename_in}.")

        layer = dataset.GetLayer(0)
        if layer is None:
            raise ValueError(f"Unable to access layer from {sFilename_in}.")

        extent = layer.GetExtent()
        if extent is None:
            raise ValueError(f"Unable to retrieve extent from {sFilename_in}.")

        if not iFlag_return_union_geometry:
            return extent

        # Create a collection of all geometries
        layer_geom_type = layer.GetGeomType()

        # Determine the target collection type using the flattened geometry type
        base_geom_type = layer_geom_type & 0xFF  # Flatten to 2D base type

        if base_geom_type == ogr.wkbPoint:
            union_geom = ogr.Geometry(ogr.wkbMultiPoint)
        elif base_geom_type == ogr.wkbLineString:
            union_geom = ogr.Geometry(ogr.wkbMultiLineString)
        elif base_geom_type == ogr.wkbPolygon:
            union_geom = ogr.Geometry(ogr.wkbMultiPolygon)
        else:
            geom_type_name = ogr.GeometryTypeToName(layer_geom_type)
            # For unsupported types, return extent without union geometry
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

            # Handle multi-geometries by adding their individual parts
            if flat_geom_type in [
                ogr.wkbMultiPoint,
                ogr.wkbMultiLineString,
                ogr.wkbMultiPolygon,
            ]:
                for i in range(geom_clone.GetGeometryCount()):
                    part = geom_clone.GetGeometryRef(i)
                    if part is not None:
                        union_geom.AddGeometry(part)
            elif flat_geom_type in [ogr.wkbPoint, ogr.wkbLineString, ogr.wkbPolygon]:
                union_geom.AddGeometry(geom_clone)

        return extent, union_geom

    finally:
        dataset = None
