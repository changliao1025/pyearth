"""Vector polygon file validation utilities using OGR."""

import os
from typing import Optional, Union
from osgeo import ogr

PathLike = Union[str, os.PathLike]

__all__ = ["gdal_validate_polygon_file"]


def gdal_validate_polygon_file(
    filename: Optional[PathLike] = None,
    print_err: bool = False,
    **kwargs,
) -> bool:
    """Validate that a vector file exists, contains layers, and all geometries are valid polygons.

    Accepts both standard Polygon and MultiPolygon geometries.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the vector polygon file (e.g. .geojson, .shp, .gpkg).
    print_err : bool, default=False
        If True, print error descriptions if validation fails.

    Returns
    -------
    bool
        True if the file can be opened and all features contain valid Polygon
        or MultiPolygon geometries, False otherwise.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")

    if filename is None:
        return False

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        if print_err:
            print(f"File does not exist: {filepath}")
        return False

    dataset = None
    layer = None
    try:
        dataset = ogr.Open(filepath, 0)
        if dataset is None:
            if print_err:
                print(f"Unable to open vector dataset: {filepath}")
            return False

        layer = dataset.GetLayer(0)
        if layer is None:
            if print_err:
                print(f"No layers found in dataset: {filepath}")
            return False

        feature_count = 0
        for feature in layer:
            feature_count += 1
            geom = feature.GetGeometryRef()
            if geom is None:
                if print_err:
                    print(f"Feature {feature.GetFID()} has null geometry")
                return False

            geom_name = geom.GetGeometryName().upper()
            if geom_name not in ("POLYGON", "MULTIPOLYGON"):
                if print_err:
                    print(
                        f"Feature {feature.GetFID()} has non-polygon geometry: {geom_name}"
                    )
                return False

            if not geom.IsValid():
                if print_err:
                    print(f"Feature {feature.GetFID()} has invalid geometry topology")
                return False

        return feature_count > 0

    except Exception as e:
        if print_err:
            print(f"Error validating polygon file: {e}")
        return False

    finally:
        layer = None
        dataset = None

