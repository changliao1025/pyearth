#############
API Reference
#############

This page provides the API reference for PyEarth. For topic-specific guides and narrative details, see:

* :doc:`gis/gis` - GIS and geospatial operations
* :doc:`system/system` - System utilities and environment settings
* :doc:`toolbox/toolbox` - Geoprocessing toolbox utilities

================
Module Structure
================

PyEarth is organized into three primary modules:

.. code-block:: text

   pyearth/
   ├── gis/                    # GIS and geospatial operations
   │   ├── envi/              # ENVI format support
   │   ├── gdal/              # GDAL wrappers for raster and vector I/O
   │   │   ├── read/          # Raster and vector readers
   │   │   └── write/         # Raster and vector writers
   │   ├── geometry/          # Spherical geometry, polygons, simplification
   │   ├── location/          # Coordinate conversion and spatial indexing
   │   ├── spatialref/        # Spatial reference and projection operations
   │   ├── idl_handler.py     # Unified Antimeridian (IDL) handler
   │   └── polar_handler.py   # Unified Polar region enclosure detector
   ├── system/                 # System utilities
   │   ├── python/            # Python and Conda environment detection
   │   ├── create_symlink.py  # Cross-platform symbolic link utility
   │   ├── define_global_variables.py # Constants and platform metadata
   │   ├── filename.py        # Path and filename parsers
   │   └── set_linker_directory.py    # Dynamic library linker config
   └── toolbox/                # Geoprocessing and scientific toolbox
       ├── analysis/          # Clip, filter, intersect, difference, union
       ├── conversion/        # Format conversion, rasterization, vectorization
       ├── data/              # Outlier removal, time series aggregation, netCDF
       ├── date/              # Julian date, leap year, day-of-year, timer
       ├── geometry/          # Regular geometric calculations
       ├── management/        # Raster and vector file management
       ├── math/              # Statistical remap, KDE, gap filling
       ├── mesh/              # Point, line, polygon, and mesh algorithms
       └── reader/            # XML parsing, line count, config reading

========================
GIS Module (pyearth.gis)
========================

ENVI Operations
---------------

.. code-block:: python

   from pyearth.gis.envi.envi_write_header import envi_write_header

GDAL Operations
---------------

Raster I/O
~~~~~~~~~~

.. code-block:: python

   # Reading rasters
   from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import (
       gdal_read_geotiff_file,
       gdal_read_geotiff_file_multiple_band,
   )
   from pyearth.gis.gdal.read.raster.gdal_read_envi_file import (
       gdal_read_envi_file,
       gdal_read_envi_file_multiple_band,
   )
   from pyearth.gis.gdal.read.raster.gdal_read_ascii_file import gdal_read_ascii_file
   from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
   from pyearth.gis.gdal.read.raster.gdal_get_raster_spatial_reference import (
       gdal_get_raster_spatial_reference_wkt,
   )

   # Writing rasters
   from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import (
       gdal_write_geotiff_file,
       gdal_write_geotiff_file_multiple_band,
   )
   from pyearth.gis.gdal.write.raster.gdal_write_envi_file import (
       gdal_write_envi_file,
       gdal_write_envi_file_multiple_band,
   )

Vector I/O
~~~~~~~~~~

.. code-block:: python

   # Reading vectors
   from pyearth.gis.gdal.read.vector.gdal_read_vector_file import (
       gdal_read_shapefile,
       gdal_read_vector_file,
   )
   from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
   from pyearth.gis.gdal.read.vector.gdal_get_vector_boundary import gdal_get_vector_boundary
   from pyearth.gis.gdal.read.vector.gdal_get_vector_spatial_reference import (
       gdal_get_vector_spatial_reference_wkt,
   )

   # Writing vectors
   from pyearth.gis.gdal.write.vector.gdal_write_wkt_to_vector_file import (
       gdal_write_wkt_to_vector_file,
   )
   from pyearth.gis.gdal.write.vector.gdal_export_point_to_vector_file import (
       export_point_to_vector_file,
       export_point_as_polygon_file,
   )

Format Support & Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pyearth.gis.gdal.gdal_raster_format_support import (
       gdal_raster_format_support,
       get_raster_format_from_filename,
       get_raster_driver_from_filename,
   )
   from pyearth.gis.gdal.gdal_vector_format_support import (
       gdal_vector_format_support,
       get_vector_format_from_filename,
       get_vector_driver_from_filename,
       has_parquet_support,
   )
   from pyearth.gis.gdal.gdal_check_file_type import gdal_check_file_type
   from pyearth.gis.gdal.gdal_check_raster_valid import gdal_check_raster_valid
   from pyearth.gis.gdal.gdal_validate_polygon_file import gdal_validate_polygon_file
   from pyearth.gis.gdal.gdal_to_numpy_datatype import (
       gdal_to_numpy_datatype,
       numpy_dtype_to_gdal,
   )

Geometry Operations
-------------------

Distance, Area, and Angles
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
       calculate_distance_based_on_longitude_latitude,
   )
   from pyearth.gis.geometry.calculate_distance_to_line import calculate_distance_to_line
   from pyearth.gis.geometry.calculate_distance_to_plane import calculate_distance_to_plane
   from pyearth.gis.geometry.calculate_polygon_area import (
       calculate_polygon_area,
       calculate_polygon_file_area,
   )
   from pyearth.gis.geometry.calculate_spherical_triangle_area import (
       calculate_spherical_triangle_area,
   )
   from pyearth.gis.geometry.calculate_angle_between_point import (
       calculate_angle_between_point,
   )
   from pyearth.gis.geometry.calculate_angle_between_vectors_degrees import (
       calculate_angle_between_vectors_degrees,
   )

Polygon Operations & Simplification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pyearth.gis.geometry.check_convex_polygon import check_convex_polygon
   from pyearth.gis.geometry.check_counter_clockwise import (
       check_counter_clockwise,
       get_polygon_orientation,
   )
   from pyearth.gis.geometry.clean_geometry import clean_geometry
   from pyearth.gis.geometry.remove_duplicate_closure import remove_duplicate_closure
   from pyearth.gis.geometry.get_output_geometry_type import get_output_geometry_type
   from pyearth.gis.geometry.douglas_peucker_geodetic import douglas_peucker_geodetic
   from pyearth.gis.geometry.visvalingam_whyatt_geodetic import (
       visvalingam_whyatt_geodetic,
   )

International Date Line Handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pyearth.gis.idl_handler import (
       IdlHandler,
       IdlStrategy,
       IdlResult,
       IDL_TOLERANCE,
   )

Polar Region Handling
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pyearth.gis.polar_handler import (
       PolarHandler,
       PoleType,
       PolarResult,
       polygon_includes_pole,
   )

Location Utilities
------------------

.. code-block:: python

   from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
       convert_longitude_latitude_to_sphere_3d,
       convert_sphere_3d_to_longitude_latitude,
   )
   from pyearth.gis.location.get_geometry_coordinates import (
       get_geometry_coordinates,
       get_polygon_exterior_coords,
   )
   from pyearth.gis.location.find_index_by_longitude_latitude import (
       find_index_by_longitude_latitude,
   )
   from pyearth.gis.location.world2Pixel import world2Pixel
   from pyearth.gis.location.get_hydrosheds_continent_from_extent import (
       get_hydrosheds_continent_from_extent,
   )

Spatial Reference
-----------------

.. code-block:: python

   from pyearth.gis.spatialref.reproject_coordinates import (
       reproject_coordinates,
       reproject_coordinates_batch,
   )
   from pyearth.gis.spatialref.utm_utility import (
       get_utm_zone,
       get_utm_epsg_code,
       get_utm_spatial_reference_wkt,
   )
   from pyearth.gis.spatialref.convert_between_degree_and_meter import (
       degree_to_meter,
       meter_to_degree,
   )
   from pyearth.gis.spatialref.is_wgs84_projection import is_wgs84_projection

==============================
System Module (pyearth.system)
==============================

Configuration & Global Variables
--------------------------------

.. code-block:: python

   from pyearth.system.define_global_variables import (
       earth_radius,
       missing_value,
       cms2cmd,
       mms2mmd,
       feet2meter,
       inch2mm,
       nmonth,
       iMonth_start,
       iMonth_end,
       sPlatform_os,
       sMachine,
       sUsername,
       sWorkspace_home,
       sWorkspace_scratch,
       slash,
       sPath_proj_lib,
       print_environment_info,
   )

File and Path Utilities
-----------------------

.. code-block:: python

   from pyearth.system.create_symlink import create_symlink
   from pyearth.system.filename import (
       get_extension_from_filename,
       get_extension_from_path,
       get_filename_without_extension,
       get_filename_from_path_without_extension,
       get_folder_path,
   )

Linker & Environment Introspection
----------------------------------

.. code-block:: python

   from pyearth.system.set_linker_directory import configure_linker_path
   from pyearth.system.python.get_python_environment import (
       get_python_environment,
       get_conda_environment,
   )

==================================
Toolbox Module (pyearth.toolbox)
==================================

Analysis Operations
-------------------

.. code-block:: python

   # Extract and clipping
   from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import (
       clip_raster_by_polygon_file,
   )
   from pyearth.toolbox.analysis.extract.clip_vector_by_polygon_file import (
       clip_vector_by_polygon_file,
   )
   from pyearth.toolbox.analysis.extract.clip_vector_by_bounding_box import (
       clip_vector_by_bounding_box,
   )
   from pyearth.toolbox.analysis.extract.filter_vector_by_attribute import (
       filter_vector_by_attribute,
   )
   from pyearth.toolbox.analysis.extract.filter_vector_by_polygon import (
       filter_vector_by_polygon,
   )

   # Intersect, Difference, Union
   from pyearth.toolbox.analysis.intersect.intersect_polygon_with_polygon_file import (
       intersect_polygon_with_polygon_file,
   )
   from pyearth.toolbox.analysis.difference.difference_polygon_with_polygon_file import (
       difference_polygon_with_polygon_file,
   )
   from pyearth.toolbox.analysis.union.union_polygon_file import union_polygon_file

Conversion Operations
---------------------

.. code-block:: python

   from pyearth.toolbox.conversion.convert_vector_format import convert_vector_format
   from pyearth.toolbox.conversion.rasterize_vector import rasterize_vector
   from pyearth.toolbox.conversion.vectorize_raster import vectorize_raster
   from pyearth.toolbox.conversion.convert_polygon_to_polyline_file import (
       convert_polygon_to_polyline_file,
   )
   from pyearth.toolbox.conversion.convert_vector_to_global_raster import (
       convert_vector_to_global_raster,
   )

Data Processing
---------------

.. code-block:: python

   from pyearth.toolbox.data.remove_outliers import remove_outliers
   from pyearth.toolbox.data.cgpercentiles import cgpercentiles
   from pyearth.toolbox.data.convert_time_series_daily_to_monthly import (
       convert_time_series_daily_to_monthly,
   )
   from pyearth.toolbox.data.list.check_if_duplicates import check_if_duplicates
   from pyearth.toolbox.data.list.list_alg import list_max, list_min
   from pyearth.toolbox.data.vector.copy_geometry_without_attributes import (
       copy_geometry_without_attributes,
   )

Date/Time Utilities
-------------------

.. code-block:: python

   from pyearth.toolbox.date.day_of_year import day_of_year
   from pyearth.toolbox.date.day_in_month import day_in_month
   from pyearth.toolbox.date.leap_year import leap_year
   from pyearth.toolbox.date.julian import to_jd, from_jd
   from pyearth.toolbox.date.dt2cal import dt2cal
   from pyearth.toolbox.date.timer import timer

Geometry Calculations
---------------------

.. code-block:: python

   from pyearth.toolbox.geometry.calculate_hexagon_area import calculate_hexagon_area

Management Operations
---------------------

.. code-block:: python

   # Raster management
   from pyearth.toolbox.management.raster.merge_raster_files import merge_raster_files
   from pyearth.toolbox.management.raster.reproject import reproject_raster
   from pyearth.toolbox.management.raster.resample import resample_raster

   # Vector management
   from pyearth.toolbox.management.vector.merge_vector_files import merge_vector_files
   from pyearth.toolbox.management.vector.merge_features import merge_features
   from pyearth.toolbox.management.vector.reproject import reproject_vector
   from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
   from pyearth.toolbox.management.vector.fields import (
       get_field_value,
       get_field_and_value,
       add_field_to_vector_file,
   )
   from pyearth.toolbox.management.vector.polygon_calculator import (
       polygon_difference,
       polygon_difference_rtree,
   )

Math Operations
---------------

.. code-block:: python

   from pyearth.toolbox.math.stat.remap import remap
   from pyearth.toolbox.math.stat.scipy_bivariate_kde import scipy_bivariate_kde
   from pyearth.toolbox.math.gap_fill_by_window import gap_fill_by_window

Mesh Classes and Algorithms
---------------------------

.. code-block:: python

   # Mesh classes
   from pyearth.toolbox.mesh.point import pypoint
   from pyearth.toolbox.mesh.line import pyline
   from pyearth.toolbox.mesh.polyline import pypolyline
   from pyearth.toolbox.mesh.polygon import pypolygon
   from pyearth.toolbox.mesh.circle import pycircle
   from pyearth.toolbox.mesh.nvector import pynvector

   # Mesh algorithms
   from pyearth.toolbox.mesh.algorithm.convert_coordinates import (
       convert_gcs_coordinates_to_meshcell,
       convert_gcs_coordinates_to_polyline,
       convert_pcs_coordinates_to_polyline,
   )
   from pyearth.toolbox.mesh.algorithm.convert_nvector_point import (
       point_to_nvector,
       nvector_to_point,
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

Reader Utilities
----------------

.. code-block:: python

   from pyearth.toolbox.reader.line_count import line_count
   from pyearth.toolbox.reader.parse_xml_file import parse_xml_file
   from pyearth.toolbox.reader.text_reader_string import text_reader_string
   from pyearth.toolbox.reader.read_configuration_file import read_configuration_file
   from pyearth.toolbox.reader.split_string_into_chunk import split_string_into_chunk
   from pyearth.toolbox.reader.beta.parse_xml_file_e3sm import parse_xml_file_e3sm

=======================
Quick Reference by Task
=======================

Working with Rasters
--------------------

.. code-block:: python

   # Read GeoTIFF
   from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
   data, geotransform, projection, nrows, ncols, nodata, extent, wkt = \
       gdal_read_geotiff_file('input.tif')

   # Write GeoTIFF
   from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import gdal_write_geotiff_file
   gdal_write_geotiff_file('output.tif', data, geotransform, projection)

   # Clip raster
   from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import (
       clip_raster_by_polygon_file,
   )
   clip_raster_by_polygon_file('input.tif', 'boundary.geojson', 'clipped.tif')

Working with Vectors
--------------------

.. code-block:: python

   # Read shapefile
   from pyearth.gis.gdal.read.vector.gdal_read_vector_file import gdal_read_shapefile
   geometries = gdal_read_shapefile('input.shp')

   # Convert format
   from pyearth.toolbox.conversion.convert_vector_format import convert_vector_format
   convert_vector_format('input.shp', 'output.geojson')

   # Merge files
   from pyearth.toolbox.management.vector.merge_vector_files import merge_vector_files
   merge_vector_files(['file1.shp', 'file2.shp'], 'merged.shp')

Coordinate Operations
---------------------

.. code-block:: python

   # Calculate great-circle distance
   from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
       calculate_distance_based_on_longitude_latitude,
   )
   distance = calculate_distance_based_on_longitude_latitude(-122.0, 47.0, -73.0, 40.0)

   # Reproject coordinates
   from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates
   x, y = reproject_coordinates(-122.0, 47.0, 'EPSG:4326', 'EPSG:32610')

   # Convert to 3D Cartesian coordinates
   from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
       convert_longitude_latitude_to_sphere_3d,
   )
   x, y, z = convert_longitude_latitude_to_sphere_3d(-122.0, 47.0)

See Also
--------

* :doc:`gis/gis` - Detailed GIS module documentation
* :doc:`system/system` - System utilities documentation
* :doc:`toolbox/toolbox` - Detailed Toolbox module documentation
* :doc:`../getting-started` - Getting started guide with examples
* :doc:`../faq` - Frequently asked questions
