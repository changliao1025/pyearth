#############
API Reference
#############

This page provides the complete API reference for PyEarth. For more detailed information about specific modules, see:

* :doc:`gis/gis` - GIS and geospatial operations
* :doc:`system/system` - System utilities
* :doc:`toolbox/toolbox` - Toolbox utilities

Module Structure
================

PyEarth is organized into three main modules:

.. code-block:: text

   pyearth/
   ├── gis/                    # GIS and geospatial operations
   │   ├── envi/              # ENVI format support
   │   ├── gdal/              # GDAL wrappers for raster and vector I/O
   │   ├── geometry/          # Geometric calculations on sphere
   │   ├── location/          # Coordinate conversion and spatial indexing
   │   └── spatialref/        # Spatial reference and projection operations
   ├── system/                 # System utilities
   │   └── python/            # Python environment detection
   └── toolbox/                # General-purpose toolbox
       ├── analysis/          # Spatial analysis operations
       ├── conversion/        # Format conversion utilities
       ├── data/              # Data processing and cleaning
       ├── date/              # Date and time utilities
       ├── geometry/          # Geometric calculations
       ├── management/        # Data management operations
       ├── math/              # Mathematical operations
       ├── mesh/              # Mesh generation and manipulation
       └── reader/            # File reading utilities

*************************
GIS Module (pyearth.gis)
*************************

ENVI Operations
===============

.. code-block:: python

   from pyearth.gis.envi.envi_write_header import envi_write_header

GDAL Operations
===============

Raster I/O
----------

.. code-block:: python

   # Reading rasters
   from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import (
       gdal_read_geotiff_file,
       gdal_read_geotiff_file_multiple_band
   )
   from pyearth.gis.gdal.read.raster.gdal_read_envi_file import (
       gdal_read_envi_file,
       gdal_read_envi_file_multiple_band
   )
   from pyearth.gis.gdal.read.raster.gdal_read_ascii_file import gdal_read_ascii_file

   # Writing rasters
   from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import (
       gdal_write_geotiff_file,
       gdal_write_geotiff_file_multiple_band
   )
   from pyearth.gis.gdal.write.raster.gdal_write_envi_file import (
       gdal_write_envi_file,
       gdal_write_envi_file_multiple_band
   )

Vector I/O
----------

.. code-block:: python

   # Reading vectors
   from pyearth.gis.gdal.read.vector.gdal_read_vector_file import gdal_read_shapefile

   # Writing vectors
   from pyearth.gis.gdal.write.vector.gdal_write_wkt_to_vector_file import (
       gdal_write_wkt_to_vector_file
   )
   from pyearth.gis.gdal.write.vector.gdal_export_point_to_vector_file import (
       export_point_to_vector_file,
       export_point_as_polygon_file
   )

Format Support
--------------

.. code-block:: python

   from pyearth.gis.gdal.gdal_raster_format_support import (
       gdal_raster_format_support,
       get_raster_format_from_filename,
       get_raster_driver_from_filename
   )
   from pyearth.gis.gdal.gdal_vector_format_support import (
       gdal_vector_format_support,
       get_vector_format_from_filename,
       get_vector_driver_from_filename
   )

Geometry Operations
===================

Distance and Area
-----------------

.. code-block:: python

   from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
       calculate_distance_based_on_longitude_latitude
   )
   from pyearth.gis.geometry.calculate_polygon_area import (
       calculate_polygon_area,
       calculate_polygon_file_area
   )
   from pyearth.gis.geometry.calculate_spherical_triangle_area import (
       calculate_spherical_triangle_area
   )

Polygon Operations
------------------

.. code-block:: python

   from pyearth.gis.geometry.check_convex_polygon import is_convex_polygon
   from pyearth.gis.geometry.check_counter_clockwise import (
       check_counter_clockwise,
       get_polygon_orientation
   )
   from pyearth.gis.geometry.clean_geometry import clean_geometry

International Date Line
-----------------------

.. code-block:: python

   from pyearth.gis.geometry.international_date_line_utility import (
       check_cross_international_date_line_geometry,
       split_international_date_line_polygon_coordinates,
       unwrap_longitudes
   )

Simplification
--------------

.. code-block:: python

   from pyearth.gis.geometry.douglas_peucker_geodetic import douglas_peucker_geodetic
   from pyearth.gis.geometry.visvalingam_whyatt_geodetic import visvalingam_whyatt_geodetic

Location Utilities
==================

.. code-block:: python

   from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import (
       convert_longitude_latitude_to_sphere_3d,
       convert_sphere_3d_to_longitude_latitude
   )
   from pyearth.gis.location.get_geometry_coordinates import (
       get_geometry_coordinates,
       get_polygon_exterior_coords
   )
   from pyearth.gis.location.find_index_by_longitude_latitude import (
       find_index_by_longitude_latitude
   )

Spatial Reference
=================

.. code-block:: python

   from pyearth.gis.spatialref.reproject_coordinates import (
       reproject_coordinates,
       reproject_coordinates_batch
   )
   from pyearth.gis.spatialref.utm_utility import (
       get_utm_zone,
       get_utm_epsg_code,
       get_utm_spatial_reference_wkt
   )
   from pyearth.gis.spatialref.convert_between_degree_and_meter import (
       degree_to_meter,
       meter_to_degree
   )

*******************************
System Module (pyearth.system)
*******************************

.. code-block:: python

   from pyearth.system.define_global_variables import (
       earth_radius,
       print_environment_info
   )
   from pyearth.system.create_symlink import create_symlink
   from pyearth.system.filename import (
       get_extension_from_filename,
       get_filename_without_extension,
       get_folder_path
   )
   from pyearth.system.python.get_python_environment import (
       get_python_environment,
       get_conda_environment
   )

***********************************
Toolbox Module (pyearth.toolbox)
***********************************

Analysis
========

.. code-block:: python

   # Extract operations
   from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import (
       clip_raster_by_polygon_file
   )
   from pyearth.toolbox.analysis.extract.clip_vector_by_polygon_file import (
       clip_vector_by_polygon_file
   )

   # Intersection operations
   from pyearth.toolbox.analysis.intersect.intersect_polygon_with_polygon_file import (
       intersect_polygon_with_polygon_file
   )

   # Difference operations
   from pyearth.toolbox.analysis.difference.difference_polygon_with_polygon_file import (
       difference_polygon_with_polygon_file
   )

Conversion
==========

.. code-block:: python

   from pyearth.toolbox.conversion.convert_vector_format import convert_vector_format
   from pyearth.toolbox.conversion.rasterize_vector import rasterize_vector
   from pyearth.toolbox.conversion.vectorize_raster import vectorize_raster

Data Processing
===============

.. code-block:: python

   from pyearth.toolbox.data.remove_outliers import remove_outliers
   from pyearth.toolbox.data.cgpercentiles import cgpercentiles
   from pyearth.toolbox.data.convert_time_series_daily_to_monthly import (
       convert_time_series_daily_to_monthly
   )

Date/Time
=========

.. code-block:: python

   from pyearth.toolbox.date.day_of_year import day_of_year
   from pyearth.toolbox.date.day_in_month import day_in_month
   from pyearth.toolbox.date.leap_year import leap_year
   from pyearth.toolbox.date.julian import to_jd, from_jd
   from pyearth.toolbox.date.timer import timer

Geometry
========

.. code-block:: python

   from pyearth.toolbox.geometry.calculate_hexagon_area import calculate_hexagon_area
   from pyearth.toolbox.geometry.create_gcs_buffer_zone import (
       create_point_buffer_zone,
       create_polyline_buffer_zone,
       create_buffer_zone_polygon_file
   )

Management
==========

.. code-block:: python

   # Raster management
   from pyearth.toolbox.management.raster.merge_rasters import merge_rasters
   from pyearth.toolbox.management.raster.reproject import reproject_raster
   from pyearth.toolbox.management.raster.resample import resample_raster

   # Vector management
   from pyearth.toolbox.management.vector.merge_files import merge_files
   from pyearth.toolbox.management.vector.merge_features import merge_features
   from pyearth.toolbox.management.vector.reproject import reproject_vector
   from pyearth.toolbox.management.vector.fields import (
       get_field_value,
       add_field_to_vector_file
   )

Math
====

.. code-block:: python

   from pyearth.toolbox.math.stat.remap import remap
   from pyearth.toolbox.math.stat.scipy_bivariate_kde import scipy_bivariate_kde
   from pyearth.toolbox.math.gap_fill_by_window import gap_fill_by_window

Mesh
====

.. code-block:: python

   # Mesh classes
   from pyearth.toolbox.mesh.point import pypoint
   from pyearth.toolbox.mesh.line import pyline
   from pyearth.toolbox.mesh.polyline import pypolyline
   from pyearth.toolbox.mesh.polygon import pypolygon
   from pyearth.toolbox.mesh.circle import pycircle
   from pyearth.toolbox.mesh.nvector import pynvector

   # Mesh generation
   from pyearth.toolbox.mesh.latlon.create_latlon_mesh import create_latlon_mesh
   from pyearth.toolbox.mesh.square.create_square_mesh import create_square_mesh
   from pyearth.toolbox.mesh.hexaon.create_hexagon_mesh import create_hexagon_mesh

   # Mesh algorithms
   from pyearth.toolbox.mesh.algorithm.find_minimal_enclosing_polygon import (
       find_minimal_enclosing_polygon
   )
   from pyearth.toolbox.mesh.algorithm.save_points_as_polygon import (
       save_points_as_polygon
   )

Reader
======

.. code-block:: python

   from pyearth.toolbox.reader.line_count import line_count
   from pyearth.toolbox.reader.parse_xml_file import parse_xml_file

************************
Quick Reference by Task
************************

Working with Rasters
====================

.. code-block:: python

   # Read GeoTIFF
   from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
   data, geotransform, projection, nrows, ncols, nodata, extent, wkt = \
       gdal_read_geotiff_file('input.tif')

   # Write GeoTIFF
   from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import gdal_write_geotiff_file
   gdal_write_geotiff_file('output.tif', data, geotransform, projection)

   # Clip raster
   from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import \
       clip_raster_by_polygon_file
   clip_raster_by_polygon_file('input.tif', 'boundary.geojson', 'clipped.tif')

Working with Vectors
====================

.. code-block:: python

   # Read shapefile
   from pyearth.gis.gdal.read.vector.gdal_read_vector_file import gdal_read_shapefile
   geometries = gdal_read_shapefile('input.shp')

   # Convert format
   from pyearth.toolbox.conversion.convert_vector_format import convert_vector_format
   convert_vector_format('input.shp', 'output.geojson')

   # Merge files
   from pyearth.toolbox.management.vector.merge_files import merge_files
   merge_files(['file1.shp', 'file2.shp'], 'merged.shp')

Coordinate Operations
=====================

.. code-block:: python

   # Calculate distance
   from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import \
       calculate_distance_based_on_longitude_latitude
   distance = calculate_distance_based_on_longitude_latitude(-122.0, 47.0, -73.0, 40.0)

   # Reproject coordinates
   from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates
   x, y = reproject_coordinates(-122.0, 47.0, 'EPSG:4326', 'EPSG:32610')

   # Convert to 3D
   from pyearth.gis.location.convert_between_longitude_latitude_and_sphere_3d import \
       convert_longitude_latitude_to_sphere_3d
   x, y, z = convert_longitude_latitude_to_sphere_3d(-122.0, 47.0)

Date/Time Operations
====================

.. code-block:: python

   # Day of year
   from pyearth.toolbox.date.day_of_year import day_of_year
   doy = day_of_year(2024, 3, 15)

   # Check leap year
   from pyearth.toolbox.date.leap_year import leap_year
   is_leap = leap_year(2024)

   # Julian date conversion
   from pyearth.toolbox.date.julian import to_jd, from_jd
   from datetime import datetime
   jd = to_jd(datetime(2024, 3, 15))
   dt = from_jd(jd)

See Also
========

* :doc:`gis/gis` - Detailed GIS module documentation
* :doc:`toolbox/toolbox` - Detailed Toolbox module documentation
* :doc:`system/system` - System utilities documentation
* :doc:`../getting-started` - Getting started guide with examples
* :doc:`../faq` - Frequently asked questions
