"""Backward-compatibility alias for pyearth.gis.idl_handler."""

from pyearth.gis.idl_handler import (  # noqa: F401
    DEFAULT_CUTOFF,
    IDL_OFFSET,
    IDL_TOLERANCE,
    IdlHandler,
    IdlResult,
    IdlStrategy,
    calculate_signed_area_shoelace,
    check_cross_international_date_line_geometry,
    check_cross_international_date_line_polygon,
    convert_international_date_line_polygon_to_unwrapped_polygon,
    split_international_date_line_polygon_coordinates,
    unwrap_longitudes,
)
