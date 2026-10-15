"""GIS core package for PyEarth: spatial reference, GDAL/ENVI I/O, location, geometry, IDL and polar handlers."""

from pyearth.gis.idl_handler import (
    DEFAULT_CUTOFF,
    IDL_OFFSET,
    IDL_TOLERANCE,
    IdlHandler,
    IdlResult,
    IdlStrategy,
)
from pyearth.gis.polar_handler import (
    LATITUDE_BOUNDARY,
    POLAR_TOLERANCE,
    PolarHandler,
    PolarResult,
    PoleType,
    polygon_includes_pole,
)

__all__ = [
    # IDL Handler
    "IdlHandler",
    "IdlStrategy",
    "IdlResult",
    "IDL_TOLERANCE",
    "IDL_OFFSET",
    "DEFAULT_CUTOFF",
    # Polar Handler
    "PolarHandler",
    "PoleType",
    "PolarResult",
    "polygon_includes_pole",
    "POLAR_TOLERANCE",
    "LATITUDE_BOUNDARY",
]
