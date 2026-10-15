"""GDAL raster and vector reading utilities."""

from pyearth.gis.gdal.read.raster import (
    gdal_get_raster_extent,
    gdal_get_raster_spatial_reference,
    gdal_get_raster_spatial_reference_wkt,
    gdal_read_ascii_file,
    gdal_read_envi_file,
    gdal_read_envi_file_multiple_band,
    gdal_read_geotiff_file,
    gdal_read_geotiff_file_multiple_band,
)
from pyearth.gis.gdal.read.vector import (
    gdal_get_vector_boundary,
    gdal_get_vector_extent,
    gdal_get_vector_spatial_reference,
    gdal_get_vector_spatial_reference_wkt,
    gdal_read_vector_file,
)

__all__ = [
    # Raster readers
    "gdal_read_geotiff_file",
    "gdal_read_geotiff_file_multiple_band",
    "gdal_read_envi_file",
    "gdal_read_envi_file_multiple_band",
    "gdal_read_ascii_file",
    "gdal_get_raster_extent",
    "gdal_get_raster_spatial_reference",
    "gdal_get_raster_spatial_reference_wkt",
    # Vector readers
    "gdal_read_vector_file",
    "gdal_get_vector_extent",
    "gdal_get_vector_boundary",
    "gdal_get_vector_spatial_reference",
    "gdal_get_vector_spatial_reference_wkt",
]
