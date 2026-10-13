"""Spherical and planar geometric operations, polygons, angles, and simplification."""

from pyearth.gis.geometry.calculate_angle_between_point import (
    calculate_angle_between_point,
)
from pyearth.gis.geometry.calculate_angle_between_point_normal import (
    calculate_angle_between_point_normal,
)
from pyearth.gis.geometry.calculate_angle_between_vectors_degrees import (
    calculate_angle_between_vectors_degrees,
)
from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
    calculate_distance_based_on_longitude_latitude,
)
from pyearth.gis.geometry.calculate_distance_to_line import calculate_distance_to_line
from pyearth.gis.geometry.calculate_distance_to_plane import calculate_distance_to_plane
from pyearth.gis.geometry.calculate_intersect_on_great_circle import (
    calculate_intersect_on_great_circle,
    find_great_circle_intersection_with_meridian,
    project_point_onto_plane,
)
from pyearth.gis.geometry.calculate_polygon_area import (
    calculate_polygon_area,
    calculate_polygon_file_area,
    haversine,
    spherical_polygon_area,
)
from pyearth.gis.geometry.calculate_spherical_triangle_area import (
    calculate_spherical_triangle_area,
)
from pyearth.gis.geometry.check_convex_polygon import (
    cross_product_2d,
    is_convex_polygon,
)
from pyearth.gis.geometry.check_counter_clockwise import (
    calculate_signed_area_shoelace,
    calculate_signed_area_spherical_polar,
    check_counter_clockwise,
    get_polygon_orientation,
    unwrap_longitudes,
)
from pyearth.gis.geometry.clean_geometry import clean_geometry
from pyearth.gis.geometry.convert_longitude_range import (
    convert_180_to_360,
    convert_360_to_180,
)
from pyearth.gis.geometry.create_box_from_longitude_latitude import (
    create_box_from_longitude_latitude,
    find_nearest_resolution,
)
from pyearth.gis.geometry.douglas_peucker_geodetic import douglas_peucker_geodetic
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import (
    extract_unique_vertices_and_connectivity,
)
from pyearth.gis.geometry.get_output_geometry_type import (
    GEOMETRY_TYPE_MAPPING,
    get_output_geometry_type,
)
from pyearth.gis.geometry.remove_duplicate_closure import remove_duplicate_closure
from pyearth.gis.geometry.visvalingam_whyatt_geodetic import (
    VWSimplifier,
    visvalingam_whyatt_geodetic,
)

__all__ = [
    # Angles
    "calculate_angle_between_point",
    "calculate_angle_between_point_normal",
    "calculate_angle_between_vectors_degrees",
    # Distances
    "calculate_distance_based_on_longitude_latitude",
    "calculate_distance_to_line",
    "calculate_distance_to_plane",
    # Great circle
    "calculate_intersect_on_great_circle",
    "find_great_circle_intersection_with_meridian",
    "project_point_onto_plane",
    # Areas
    "calculate_polygon_area",
    "calculate_polygon_file_area",
    "spherical_polygon_area",
    "calculate_spherical_triangle_area",
    "haversine",
    # Polygon checks & cleaning
    "is_convex_polygon",
    "cross_product_2d",
    "check_counter_clockwise",
    "get_polygon_orientation",
    "calculate_signed_area_shoelace",
    "calculate_signed_area_spherical_polar",
    "unwrap_longitudes",
    "remove_duplicate_closure",
    "clean_geometry",
    # Coordinate range & boxes
    "convert_360_to_180",
    "convert_180_to_360",
    "create_box_from_longitude_latitude",
    "find_nearest_resolution",
    # Simplification
    "douglas_peucker_geodetic",
    "visvalingam_whyatt_geodetic",
    "VWSimplifier",
    # Mesh vertices & topology
    "extract_unique_vertices_and_connectivity",
    "get_output_geometry_type",
    "GEOMETRY_TYPE_MAPPING",
]
