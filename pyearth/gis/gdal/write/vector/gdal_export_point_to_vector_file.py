"""Export points and point collections to vector GIS files."""

import os
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

from osgeo import ogr, osr

from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area

ogr.UseExceptions()
osr.UseExceptions()

PathLike = Union[str, os.PathLike]

__all__ = [
    "export_point_to_vector_file",
    "export_point_as_polygon_file",
]


def _resolve_spatial_reference(
    srs: Optional[Union[osr.SpatialReference, int, str]], default_epsg: int = 4326
) -> osr.SpatialReference:
    """Resolve an EPSG code, WKT string, or SpatialReference instance."""
    if srs is None:
        spatial_ref = osr.SpatialReference()
        spatial_ref.ImportFromEPSG(default_epsg)
        return spatial_ref
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
    spatial_ref.ImportFromEPSG(default_epsg)
    return spatial_ref


def _extract_point_coords(point: Any, projected: bool) -> Tuple[float, float]:
    """Extract (x, y) or (lon, lat) from points (objects, dicts, tuples)."""
    if isinstance(point, (tuple, list)):
        if len(point) < 2:
            raise ValueError(
                f"Point tuple or list must contain at least 2 elements: {point}"
            )
        return float(point[0]), float(point[1])

    if isinstance(point, dict):
        if projected:
            x = point.get("dx", point.get("x"))
            y = point.get("dy", point.get("y"))
        else:
            x = point.get(
                "dLongitude_degree",
                point.get("longitude", point.get("lon", point.get("x"))),
            )
            y = point.get(
                "dLatitude_degree",
                point.get("latitude", point.get("lat", point.get("y"))),
            )
        if x is not None and y is not None:
            return float(x), float(y)
        raise ValueError(f"Could not extract coordinates from point dict: {point}")

    if projected:
        if hasattr(point, "dx") and hasattr(point, "dy"):
            return float(point.dx), float(point.dy)
        if hasattr(point, "x") and hasattr(point, "y"):
            return float(point.x), float(point.y)
    else:
        if hasattr(point, "dLongitude_degree") and hasattr(point, "dLatitude_degree"):
            return float(point.dLongitude_degree), float(point.dLatitude_degree)
        if hasattr(point, "longitude") and hasattr(point, "latitude"):
            return float(point.longitude), float(point.latitude)
        if hasattr(point, "x") and hasattr(point, "y"):
            return float(point.x), float(point.y)

    raise ValueError(f"Point object missing coordinate attributes: {point}")


def export_point_to_vector_file(
    aPoint_in: List[Any],
    sFilename_vector_in: PathLike,
    iFlag_projected_in: Optional[int] = None,
    pSpatial_reference_in: Optional[Union[osr.SpatialReference, int, str]] = None,
    aAttribute_data: Optional[List[Any]] = None,
    overwrite: bool = True,
) -> None:
    """Export a list of points to a vector GIS file.

    Parameters
    ----------
    aPoint_in : list
        List of point objects, dicts, or tuples with coordinates.
    sFilename_vector_in : PathLike
        Output filename. Format determined from extension (GeoJSON, SHP, GPKG, etc.).
    iFlag_projected_in : int, optional
        0 for geographic coordinates (default), 1 for projected coordinates.
    pSpatial_reference_in : osr.SpatialReference, int, or str, optional
        Spatial reference system. Defaults to EPSG:4326 (WGS84).
    aAttribute_data : list, optional
        Optional attribute data (connectivity, ID, or value) for each point.
    overwrite : bool, optional
        If True, overwrite existing file. Default is True.

    Raises
    ------
    ValueError
        If point list is empty or coordinates cannot be extracted.
    FileExistsError
        If output file exists and overwrite is False.
    RuntimeError
        If dataset, layer, or feature creation fails.
    """
    if not aPoint_in:
        raise ValueError("Point list cannot be empty.")

    filename_str = os.fspath(sFilename_vector_in)
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

    is_projected = bool(iFlag_projected_in == 1)
    spatial_ref = _resolve_spatial_reference(pSpatial_reference_in)

    has_attributes = aAttribute_data is not None
    attribute_data: List[Any] = aAttribute_data if aAttribute_data is not None else []

    dataset = None
    layer = None
    feature = None

    try:
        dataset = driver.CreateDataSource(filename_str)
        if dataset is None:
            raise RuntimeError(f"Could not create output dataset: {filename_str}")

        layer = dataset.CreateLayer("point", spatial_ref, ogr.wkbPoint)
        if layer is None:
            raise RuntimeError(
                f"Could not create layer in output dataset: {filename_str}"
            )

        layer.CreateField(ogr.FieldDefn("pointid", ogr.OFTInteger64))
        if has_attributes:
            layer.CreateField(ogr.FieldDefn("connectivity", ogr.OFTInteger64))

        layer_defn = layer.GetLayerDefn()

        for point_id, pPoint in enumerate(aPoint_in):
            x, y = _extract_point_coords(pPoint, is_projected)

            point_geom = ogr.Geometry(ogr.wkbPoint)
            point_geom.AddPoint(x, y)

            feat = ogr.Feature(layer_defn)
            feat.SetGeometry(point_geom)
            feat.SetField("pointid", point_id + 1)

            if has_attributes and point_id < len(attribute_data):
                feat.SetField("connectivity", int(attribute_data[point_id]))

            result = layer.CreateFeature(feat)
            if result != ogr.OGRERR_NONE:
                raise RuntimeError(
                    f"Could not create feature for point index {point_id} in {filename_str}."
                )

        dataset.FlushCache()
    finally:
        feature = None
        layer = None
        dataset = None


def export_point_as_polygon_file(
    aPoint_in: List[Any],
    sFilename_out: PathLike,
    pSpatial_reference_in: Optional[Union[osr.SpatialReference, int, str]] = None,
    overwrite: bool = True,
) -> None:
    """Export points as a closed polygon to a vector GIS file.

    Parameters
    ----------
    aPoint_in : list
        List of point objects, dicts, or tuples with longitude/latitude coordinates.
    sFilename_out : PathLike
        Output filename. Format determined from extension.
    pSpatial_reference_in : osr.SpatialReference, int, or str, optional
        Spatial reference system. Defaults to EPSG:4326.
    overwrite : bool, optional
        If True, overwrite existing file. Default is True.

    Raises
    ------
    ValueError
        If point list contains fewer than 3 points.
    FileExistsError
        If output file exists and overwrite is False.
    RuntimeError
        If dataset, layer, or feature creation fails.
    """
    if not aPoint_in or len(aPoint_in) < 3:
        raise ValueError("Polygon requires at least 3 points.")

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

    spatial_ref = _resolve_spatial_reference(pSpatial_reference_in)

    dataset = None
    layer = None
    feature = None

    try:
        dataset = driver.CreateDataSource(filename_str)
        if dataset is None:
            raise RuntimeError(f"Could not create output dataset: {filename_str}")

        layer = dataset.CreateLayer("polygon", spatial_ref, geom_type=ogr.wkbPolygon)
        if layer is None:
            raise RuntimeError(
                f"Could not create layer in output dataset: {filename_str}"
            )

        area_field = ogr.FieldDefn("area", ogr.OFTReal)
        area_field.SetWidth(20)
        area_field.SetPrecision(2)
        layer.CreateField(area_field)

        ring = ogr.Geometry(ogr.wkbLinearRing)
        lon_list = []
        lat_list = []

        for pPoint in aPoint_in:
            lon, lat = _extract_point_coords(pPoint, projected=False)
            ring.AddPoint(lon, lat)
            lon_list.append(lon)
            lat_list.append(lat)

        # Close ring if first and last point are not identical
        if lon_list[0] != lon_list[-1] or lat_list[0] != lat_list[-1]:
            ring.AddPoint(lon_list[0], lat_list[0])

        area = calculate_polygon_area(lon_list, lat_list)

        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)

        layer_defn = layer.GetLayerDefn()
        feature = ogr.Feature(layer_defn)
        feature.SetGeometry(polygon)
        feature.SetField("area", float(area))

        result = layer.CreateFeature(feature)
        if result != ogr.OGRERR_NONE:
            raise RuntimeError(f"Could not create polygon feature in {filename_str}.")

        dataset.FlushCache()
    finally:
        feature = None
        layer = None
        dataset = None
