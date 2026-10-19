"""Vector dataset reading utilities using OGR."""

import os
from typing import Any, Dict, List, Optional, Tuple, Union

from osgeo import ogr, osr

ogr.UseExceptions()
osr.UseExceptions()

PathLike = Union[str, os.PathLike]

__all__ = [
    "VectorDataset",
    "gdal_read_vector_file",
]


class VectorDataset:
    """Container for vector dataset features, metadata, and attributes.

    Supports indexing and unpacking as a 7-tuple:
    (geometries, spatial_reference, feature_count, min_x, max_x, min_y, max_y)
    as well as dictionary-style and attribute access.
    """

    def __init__(
        self,
        geometries: List[ogr.Geometry],
        spatial_reference: Optional[osr.SpatialReference],
        feature_count: int,
        min_x: float,
        max_x: float,
        min_y: float,
        max_y: float,
        attributes: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        self.geometries = geometries
        self.spatial_reference = spatial_reference
        self.feature_count = feature_count
        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.attributes = attributes if attributes is not None else []

    def __iter__(self):
        return iter(
            (
                self.geometries,
                self.spatial_reference,
                self.feature_count,
                self.min_x,
                self.max_x,
                self.min_y,
                self.max_y,
            )
        )

    def __getitem__(self, item: Union[int, str]) -> Any:
        if isinstance(item, int):
            return (
                self.geometries,
                self.spatial_reference,
                self.feature_count,
                self.min_x,
                self.max_x,
                self.min_y,
                self.max_y,
            )[item]
        if isinstance(item, str) and hasattr(self, item):
            return getattr(self, item)
        raise KeyError(item)

    def __len__(self) -> int:
        return 7

    def __repr__(self) -> str:
        return (
            f"VectorDataset(features={self.feature_count}, "
            f"bounds=({self.min_x}, {self.min_y}, {self.max_x}, {self.max_y}))"
        )


def gdal_read_vector_file(
    filename: Optional[PathLike] = None,
    metadata_only: bool = False,
    **kwargs,
) -> Union[
    Tuple[Optional[osr.SpatialReference], int, float, float, float, float],
    VectorDataset,
]:
    """Read a vector dataset (Shapefile, GeoJSON, GPKG, etc.) using OGR.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the vector file.
    metadata_only : bool, default=False
        If True, return only metadata (spatial reference, feature count, extent)
        without extracting feature geometries.

    Returns
    -------
    VectorDataset or tuple
        If metadata_only=True:
            (spatial_reference, feature_count, min_x, max_x, min_y, max_y)
        If metadata_only=False:
            VectorDataset instance (unpackable as 7-tuple:
            (geometries, spatial_reference, feature_count, min_x, max_x, min_y, max_y))

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    RuntimeError
        If OGR cannot open the file.
    ValueError
        If the dataset has no layers.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")
    if "iFlag_metadata_only" in kwargs:
        metadata_only = bool(kwargs["iFlag_metadata_only"])

    if filename is None:
        raise ValueError("Filename must be provided.")

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist.")

    dataset = ogr.Open(filepath, 0)  # 0 = read-only
    if dataset is None:
        raise RuntimeError(f"Unable to open vector dataset {filepath}.")

    layer = None
    try:
        layer = dataset.GetLayer(0)
        if layer is None:
            raise ValueError(f"No layers found in vector dataset {filepath}.")

        feature_count = layer.GetFeatureCount()
        spatial_ref = layer.GetSpatialRef()

        extent = layer.GetExtent()
        if extent is not None:
            min_x, max_x, min_y, max_y = extent
        else:
            min_x = max_x = min_y = max_y = 0.0

        if metadata_only:
            return spatial_ref, feature_count, min_x, max_x, min_y, max_y

        geometries: List[ogr.Geometry] = []
        attributes: List[Dict[str, Any]] = []

        layer.ResetReading()
        for feature in layer:
            if feature is None:
                continue

            geom = feature.GetGeometryRef()
            if geom is not None:
                geometries.append(geom.Clone())

            # Extract feature fields
            field_count = feature.GetFieldCount()
            if field_count > 0:
                record = {}
                for i in range(field_count):
                    field_defn = feature.GetFieldDefnRef(i)
                    if field_defn is not None:
                        field_name = field_defn.GetName()
                        record[field_name] = feature.GetField(i)
                attributes.append(record)

        return VectorDataset(
            geometries=geometries,
            spatial_reference=spatial_ref,
            feature_count=feature_count,
            min_x=min_x,
            max_x=max_x,
            min_y=min_y,
            max_y=max_y,
            attributes=attributes,
        )
    finally:
        layer = None
        dataset = None
