"""GDAL/OGR format support, validation, type conversions, and I/O."""

from pyearth.gis.gdal.gdal_check_file_type import gdal_check_file_type
from pyearth.gis.gdal.gdal_check_raster_valid import gdal_check_raster_valid
from pyearth.gis.gdal.gdal_raster_format_support import (
    SUPPORTED_RASTER_FORMATS,
    gdal_raster_format_support,
    get_available_raster_formats,
    get_extension_from_raster_format,
    get_raster_driver_from_filename,
    get_raster_driver_from_format,
    get_raster_format_from_extension,
    get_raster_format_from_filename,
    print_supported_raster_formats,
)
from pyearth.gis.gdal.gdal_to_numpy_datatype import (
    gdal_to_numpy_datatype,
    numpy_dtype_to_gdal,
    numpy_to_gdal_type,
)
from pyearth.gis.gdal.gdal_validate_polygon_file import gdal_validate_polygon_file
from pyearth.gis.gdal.gdal_vector_format_support import (
    SUPPORTED_VECTOR_FORMATS,
    check_parquet_support,
    gdal_vector_format_support,
    get_available_vector_formats,
    get_extension_from_vector_format,
    get_vector_driver_from_filename,
    get_vector_driver_from_format,
    get_vector_format_from_extension,
    get_vector_format_from_filename,
    has_parquet_support,
    print_supported_vector_formats,
)
from pyearth.gis.gdal.read import (
    gdal_get_raster_extent,
    gdal_get_raster_spatial_reference,
    gdal_get_raster_spatial_reference_wkt,
    gdal_get_vector_boundary,
    gdal_get_vector_extent,
    gdal_get_vector_spatial_reference,
    gdal_get_vector_spatial_reference_wkt,
    gdal_read_ascii_file,
    gdal_read_envi_file,
    gdal_read_envi_file_multiple_band,
    gdal_read_geotiff_file,
    gdal_read_geotiff_file_multiple_band,
    gdal_read_vector_file,
)
from pyearth.gis.gdal.write import (
    export_point_as_polygon_file,
    export_point_to_vector_file,
    gdal_write_envi_file,
    gdal_write_envi_file_multiple_band,
    gdal_write_geotiff_file,
    gdal_write_geotiff_file_multiple_band,
    gdal_write_wkt_to_vector_file,
)

__all__ = [
    # File type inspection and validation
    "gdal_check_file_type",
    "gdal_check_raster_valid",
    "gdal_validate_polygon_file",
    # Raster formats
    "SUPPORTED_RASTER_FORMATS",
    "gdal_raster_format_support",
    "get_available_raster_formats",
    "get_raster_format_from_extension",
    "get_raster_format_from_filename",
    "get_extension_from_raster_format",
    "get_raster_driver_from_format",
    "get_raster_driver_from_filename",
    "print_supported_raster_formats",
    # Vector formats
    "SUPPORTED_VECTOR_FORMATS",
    "gdal_vector_format_support",
    "get_available_vector_formats",
    "get_vector_format_from_extension",
    "get_vector_format_from_filename",
    "get_extension_from_vector_format",
    "get_vector_driver_from_format",
    "get_vector_driver_from_filename",
    "print_supported_vector_formats",
    "check_parquet_support",
    "has_parquet_support",
    # Type conversion
    "gdal_to_numpy_datatype",
    "numpy_dtype_to_gdal",
    "numpy_to_gdal_type",
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
    # Raster writers
    "gdal_write_geotiff_file",
    "gdal_write_geotiff_file_multiple_band",
    "gdal_write_envi_file",
    "gdal_write_envi_file_multiple_band",
    # Vector writers
    "export_point_to_vector_file",
    "export_point_as_polygon_file",
    "gdal_write_wkt_to_vector_file",
]
