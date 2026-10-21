################
Toolbox Module
################

The Toolbox module provides a comprehensive collection of utilities for common Earth science tasks, including spatial analysis, format conversion, date/time operations, mesh algorithms, data management, and file reading.

*************************
Analysis Operations
*************************

Extract Operations
==================

Clipping and Filtering
-----------------------

* **clip_raster_by_polygon_file** - Clip raster to polygon boundary
* **clip_vector_by_polygon_file** - Clip vector features by polygon boundary
* **clip_vector_by_bounding_box** - Clip vector features by rectangular bounding box
* **exclude_vector_by_polygon_file** - Exclude vector features within polygon
* **exclude_vector_by_polygon_files** - Exclude features by multiple polygons
* **filter_vector_by_attribute** - Filter vector features by attribute values
* **filter_vector_by_polygon** - Filter vector features spatially by polygon

.. code-block:: python

   from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import (
       clip_raster_by_polygon_file,
   )

   # Clip raster to watershed boundary
   clip_raster_by_polygon_file(
       'dem.tif',
       'watershed.geojson',
       'dem_clipped.tif'
   )

Intersection Operations
========================

* **intersect_polygon_with_polygon_file** - Intersect two polygon layers
* **intersect_polyline_with_polygon_files** - Intersect polylines with polygon files

Difference Operations
======================

* **difference_polygon_with_polygon_file** - Calculate polygon difference
* **difference_polyline_with_polygon_file** - Calculate polyline difference with polygon

Union Operations
================

* **union_polygon_file** - Calculate geometric union of features in a polygon file

Image Processing
================

* **raster_process** - Image operations and filtering on raster data

*************************
Conversion Operations
*************************

Format Conversion
=================

* **convert_vector_format** - Convert between vector formats (Shapefile, GeoJSON, etc.)
* **convert_polygon_to_polyline_file** - Convert polygon to polyline features

Raster/Vector Conversion
=========================

* **rasterize_vector** - Convert vector features to raster format
* **vectorize_raster** - Convert raster to vector polygons
* **convert_vector_to_global_raster** - Convert vector to global raster grid

.. code-block:: python

   from pyearth.toolbox.conversion.rasterize_vector import rasterize_vector

   # Rasterize polygon shapefile
   rasterize_vector(
       'polygons.shp',
       'output_raster.tif',
       dResolution_x_in=0.01,
       dResolution_y_in=0.01,
       sAttribute_name='value'
   )

*************************
Data Processing
*************************

Statistical Operations
======================

* **cgpercentiles** - Calculate percentiles from data array
* **remove_outliers** - Remove outliers using various methods (IQR, Z-score, percentile)

Time Series Operations
======================

* **convert_time_series_daily_to_monthly** - Aggregate daily data to monthly values

Data Cleaning & List Utilities
==============================

* **gap_fill_by_window** - Fill missing data using window-based interpolation
* **check_if_duplicates** - Check for duplicate values in sequences
* **list_max**, **list_min** - Extract extremal values from sequences

Ocean/Land Utilities
====================

* **define_land_ocean_mask** - Create land/ocean mask from Natural Earth data
* **create_land_ocean_vector_mask_naturalearth** - Generate vector land/ocean mask

Vector Data Operations
======================

* **copy_geometry_without_attributes** - Copy geometries without attribute table

NetCDF Utilities (Beta)
=======================

* **add_variable_to_netcdf** - Add a new data variable into an existing NetCDF file
* **replace_variable_in_netcdf** - Update/replace a variable within a NetCDF file

.. code-block:: python

   from pyearth.toolbox.data.remove_outliers import remove_outliers
   import numpy as np

   # Remove outliers from data
   data = np.array([1, 2, 3, 4, 5, 100, 6, 7, 8])
   clean_data = remove_outliers(data, method='iqr', threshold=1.5)

*************************
Date/Time Utilities
*************************

Date Calculations
=================

* **day_in_month** - Get number of days in a month
* **day_of_year** - Calculate day of year (1-365/366)
* **leap_year** - Check if year is a leap year

Julian Date Conversion
======================

* **to_jd** - Convert datetime to Julian date
* **from_jd** - Convert Julian date to datetime

Calendar Conversion
===================

* **dt2cal** - Convert numpy datetime64 array to calendar dates

Timer Utilities
===============

* **timer** - Context manager and decorator for timing code execution

.. code-block:: python

   from pyearth.toolbox.date.day_of_year import day_of_year
   from pyearth.toolbox.date.leap_year import leap_year
   from pyearth.toolbox.date.timer import timer

   # Calculate day of year
   doy = day_of_year(2026, 3, 15)  # Returns 74

   # Check leap year
   is_leap = leap_year(2026)  # Returns False

   # Time code execution
   with timer("Processing"):
       pass

*************************
Geometry Operations
*************************

Geometric Calculations
======================

* **calculate_hexagon_area** - Calculate area of regular hexagon based on edge length

.. code-block:: python

   from pyearth.toolbox.geometry.calculate_hexagon_area import calculate_hexagon_area

   # Calculate hexagon area (in square meters)
   area = calculate_hexagon_area(dLength_edge_in=1000.0)

*************************
Management Operations
*************************

Raster Management
=================

Transformation Operations
-------------------------

* **reproject_raster** - Reproject raster to different coordinate system
* **resample_raster** - Resample raster to different resolution
* **merge_raster_files** - Mosaic multiple raster files into a single composite

.. code-block:: python

   from pyearth.toolbox.management.raster.merge_raster_files import merge_raster_files

   # Merge multiple raster files
   merge_raster_files(
       ['tile1.tif', 'tile2.tif', 'tile3.tif'],
       'merged_output.tif'
   )

Vector Management
=================

File Operations
---------------

* **merge_vector_files** - Merge multiple vector files into one
* **merge_features** - Merge/dissolve features within a file
* **remove_small_polygon** - Remove polygons below size threshold
* **reproject_vector** - Reproject vector to different coordinate system

Field Operations
----------------

* **get_field_value** - Extract field values from vector file
* **get_field_and_value** - Get field names and values
* **add_field_to_vector_file** - Add new field to attribute table

Polygon Calculations
--------------------

* **polygon_difference** - Calculate difference between polygon layers
* **polygon_difference_rtree** - Optimized difference using R-tree index
* **polygon_calculator** - Various polygon calculation operations

.. code-block:: python

   from pyearth.toolbox.management.vector.merge_vector_files import merge_vector_files

   # Merge multiple shapefiles
   merge_vector_files(
       ['region1.shp', 'region2.shp', 'region3.shp'],
       'merged_regions.shp'
   )

*************************
Math Operations
*************************

Statistical Functions
=====================

* **remap** - Remap value from one numeric range to another
* **scipy_bivariate_kde** - Bivariate kernel density estimation

Array Operations
================

* **gap_fill_by_window** - Fill gaps using moving window interpolation

.. code-block:: python

   from pyearth.toolbox.math.stat.remap import remap

   # Remap value from [0, 100] to [0, 1]
   normalized = remap(50, old_min=0, old_max=100, new_min=0, new_max=1)

*************************
Mesh Operations
*************************

Mesh Classes
============

Point and Line Objects
----------------------

* **pypoint** - Point class with longitude/latitude coordinates
* **pynvector** - Normal vector class for sphere calculations
* **pyline** - Great circle line segment class
* **pypolyline** - Polyline class containing multiple segments

Polygon Objects
---------------

* **pypolygon** - Spherical polygon class with ordered vertices
* **pycircle** - Circle class on a sphere

.. code-block:: python

   from pyearth.toolbox.mesh.point import pypoint
   from pyearth.toolbox.mesh.polygon import pypolygon

   # Create points
   p1 = pypoint({'dLongitude_degree': -122.0, 'dLatitude_degree': 47.0})
   p2 = pypoint({'dLongitude_degree': -121.0, 'dLatitude_degree': 47.0})
   p3 = pypoint({'dLongitude_degree': -121.0, 'dLatitude_degree': 48.0})

   # Create polygon
   polygon = pypolygon(aPoint=[p1, p2, p3])

Mesh Algorithms
===============

* **convert_gcs_coordinates_to_meshcell** - Convert geographic coordinates to mesh cell
* **convert_gcs_coordinates_to_polyline** - Convert geographic coordinates to polyline
* **convert_pcs_coordinates_to_polyline** - Convert projected coordinates to polyline
* **point_to_nvector**, **nvector_to_point** - Adapter between point and n-vector representations
* **find_minimal_enclosing_polygon** - Compute minimal enclosing polygon (convex hull) of points
* **save_points_as_polygon** - Save point array as polygon geometry file
* **slerp** - Spherical linear interpolation between points
* **split_line_by_length** - Split line segments by length threshold

.. note::
   Advanced grid generation (lat-lon, square, and hexagonal meshes, MPAS mesh generation) is provided by the companion package **pyearthmesh**.

*************************
Reader Utilities
*************************

File Reading
============

* **line_count** - Count lines in text file efficiently
* **parse_xml_file** - Parse XML configuration files
* **text_reader_string** - Read text file and return as string
* **read_configuration_file** - Read and parse key-value configuration files
* **split_string_into_chunk** - Split string into fixed-size chunks

Model XML Utilities (Beta)
==========================

* **parse_xml_file_e3sm** - Parse E3SM model XML files (atmosphere and land configuration)

.. code-block:: python

   from pyearth.toolbox.reader.line_count import line_count
   from pyearth.toolbox.reader.parse_xml_file import parse_xml_file

   # Count lines in large file
   num_lines = line_count('large_file.txt')

   # Parse XML configuration
   config = parse_xml_file('config.xml')

*************************
Usage Examples
*************************

Example 1: Clip and Rasterize Workflow
=======================================

.. code-block:: python

   from pyearth.toolbox.analysis.extract.clip_vector_by_polygon_file import (
       clip_vector_by_polygon_file,
   )
   from pyearth.toolbox.conversion.rasterize_vector import rasterize_vector

   # Clip stream network to watershed
   clip_vector_by_polygon_file(
       'streams.shp',
       'watershed.geojson',
       'streams_clipped.shp'
   )

   # Rasterize clipped streams
   rasterize_vector(
       'streams_clipped.shp',
       'streams_raster.tif',
       dResolution_x_in=30.0,
       dResolution_y_in=30.0,
       sAttribute_name='order'
   )

Example 2: Time Series Processing
==================================

.. code-block:: python

   from pyearth.toolbox.data.convert_time_series_daily_to_monthly import (
       convert_time_series_daily_to_monthly,
   )
   from pyearth.toolbox.data.remove_outliers import remove_outliers
   import numpy as np

   # Daily precipitation data (365 days)
   daily_precip = np.random.gamma(2, 2, 365)

   # Remove outliers
   clean_precip = remove_outliers(daily_precip, method='iqr')

   # Convert to monthly
   monthly_precip = convert_time_series_daily_to_monthly(
       clean_precip,
       iYear_start=2026,
       iMonth_start=1,
       iDay_start=1
   )

Example 3: Mesh Algorithms and Export
=====================================

.. code-block:: python

   from pyearth.toolbox.mesh.point import pypoint
   from pyearth.toolbox.mesh.algorithm.find_minimal_enclosing_polygon import (
       find_minimal_enclosing_polygon,
   )
   from pyearth.toolbox.mesh.algorithm.save_points_as_polygon import (
       save_points_as_polygon,
   )

   points = [
       pypoint({'dLongitude_degree': -122.0, 'dLatitude_degree': 47.0}),
       pypoint({'dLongitude_degree': -121.0, 'dLatitude_degree': 47.0}),
       pypoint({'dLongitude_degree': -121.5, 'dLatitude_degree': 48.0}),
   ]

   hull = find_minimal_enclosing_polygon(points)
   save_points_as_polygon(points, 'points_polygon.geojson')

Example 4: Merge and Reproject Vectors
=======================================

.. code-block:: python

   from pyearth.toolbox.management.vector.merge_vector_files import merge_vector_files
   from pyearth.toolbox.management.vector.reproject import reproject_vector

   # Merge multiple region files
   merge_vector_files(
       ['region_north.shp', 'region_south.shp', 'region_east.shp'],
       'all_regions.shp'
   )

   # Reproject to UTM
   reproject_vector(
       'all_regions.shp',
       'all_regions_utm.shp',
       'EPSG:32610'  # UTM Zone 10N
   )
