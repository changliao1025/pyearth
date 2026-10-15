"""Mesh algorithms: coordinate conversions, nvector adapter, polygon operations, slerp."""

from pyearth.toolbox.mesh.algorithm.convert_coordinates import (
    convert_gcs_coordinates_to_meshcell,
    convert_gcs_coordinates_to_polyline,
    convert_pcs_coordinates_to_polyline,
)
from pyearth.toolbox.mesh.algorithm.convert_nvector_point import (
    nvector_to_point,
    point_to_nvector,
)
from pyearth.toolbox.mesh.algorithm.find_minimal_enclosing_polygon import (
    find_minimal_enclosing_polygon,
)
from pyearth.toolbox.mesh.algorithm.save_points_as_polygon import (
    save_points_as_polygon,
)
from pyearth.toolbox.mesh.algorithm.split_by_length import (
    slerp,
    split_line_by_length,
)

__all__ = [
    "convert_gcs_coordinates_to_meshcell",
    "convert_gcs_coordinates_to_polyline",
    "convert_pcs_coordinates_to_polyline",
    "point_to_nvector",
    "nvector_to_point",
    "find_minimal_enclosing_polygon",
    "save_points_as_polygon",
    "slerp",
    "split_line_by_length",
]
