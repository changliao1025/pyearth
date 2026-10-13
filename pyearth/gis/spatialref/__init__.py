"""Spatial reference, coordinate system utilities, and reprojections."""

from pyearth.gis.spatialref.convert_between_degree_and_meter import (
    degree_to_meter,
    meter_to_degree,
)
from pyearth.gis.spatialref.is_wgs84_projection import is_wgs84_projection
from pyearth.gis.spatialref.reproject_coordinates import (
    reproject_coordinates,
    reproject_coordinates_batch,
)
from pyearth.gis.spatialref.utm_utility import (
    get_utm_epsg_code,
    get_utm_spatial_reference_wkt,
    get_utm_zone,
)

__all__ = [
    # UTM utilities
    "get_utm_zone",
    "get_utm_epsg_code",
    "get_utm_spatial_reference_wkt",
    # Degree-meter conversions
    "degree_to_meter",
    "meter_to_degree",
    # Projection validation
    "is_wgs84_projection",
    # Reprojections
    "reproject_coordinates",
    "reproject_coordinates_batch",
]
