"""Vector reading and inspection utilities using OGR."""

from pyearth.gis.gdal.read.vector.gdal_get_vector_boundary import (
    gdal_get_vector_boundary,
)
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import (
    gdal_get_vector_extent,
)
from pyearth.gis.gdal.read.vector.gdal_get_vector_spatial_reference import (
    gdal_get_vector_spatial_reference,
    gdal_get_vector_spatial_reference_wkt,
)
from pyearth.gis.gdal.read.vector.gdal_read_vector_file import (
    gdal_read_vector_file,
)

__all__ = [
    "gdal_read_vector_file",
    "gdal_get_vector_extent",
    "gdal_get_vector_boundary",
    "gdal_get_vector_spatial_reference",
    "gdal_get_vector_spatial_reference_wkt",
]
