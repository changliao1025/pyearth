"""GIS location, coordinate conversion, and grid indexing utilities."""

from pyearth.gis.location.Google_MetersPerPixel import (
    google_meters_per_pixel,
)
from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
    convert_longitude_latitude_to_sphere_3d,
    convert_sphere_3d_to_longitude_latitude,
)
from pyearth.gis.location.find_index_by_longitude_latitude import (
    find_index_by_longitude_latitude,
)
from pyearth.gis.location.get_geometry_coordinates import (
    get_geometry_coordinates,
    get_linearring_coords,
    get_linestring_coords,
    get_multipolygon_exterior_coords,
    get_point_coords,
    get_polygon_exterior_coords,
)
from pyearth.gis.location.get_hydrosheds_continent_from_extent import (
    get_hydrosheds_continent_from_extent,
)
from pyearth.gis.location.world2Pixel import (
    pixel_to_world,
    world_to_pixel,
)

__all__ = [
    # Google tile resolutions
    "google_meters_per_pixel",
    # World and pixel conversions
    "world_to_pixel",
    "pixel_to_world",
    # 3D sphere conversions
    "convert_longitude_latitude_to_sphere_3d",
    "convert_sphere_3d_to_longitude_latitude",
    # Grid indexing
    "find_index_by_longitude_latitude",
    # HydroSHEDS extent mapping
    "get_hydrosheds_continent_from_extent",
    # Geometry coordinates extraction
    "get_geometry_coordinates",
    "get_polygon_exterior_coords",
    "get_multipolygon_exterior_coords",
    "get_linestring_coords",
    "get_point_coords",
    "get_linearring_coords",
]
