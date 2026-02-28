################
Toolbox Module
################

The Toolbox module provides a comprehensive collection of utilities for common Earth science tasks, including data analysis, format conversion, date/time operations, geometry calculations, and file management.

*************************
Analysis Operations
*************************

Extract Operations
==================

Clipping and Filtering
-----------------------

* **clip_raster_by_polygon_file** - Clip raster to polygon boundary
* **clip_vector_by_polygon_file** - Clip vector features by polygon
* **clip_vector_by_polygon_files** - Clip vector by multiple polygon files
* **exclude_vector_by_polygon_file** - Exclude vector features within polygon
* **exclude_vector_by_polygon_files** - Exclude features by multiple polygons
* **filter_vector_by_attribute** - Filter vector features by attribute values
* **filter_vector_by_polygon** - Filter vector features spatially by polygon

.. code-block:: python

   from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import \
       clip_raster_by_polygon_file

   # Clip raster to watershed boundary
   clip_raster_by_polygon_file(
       'dem.tif',
       'watershed.geojson',
       'dem_clipped.tif'
   )

Intersection Operations
========================

* **intersect_polygon_with_polygon_file** - Intersect two polygon layers
* **intersect_polyline_with_polygon_file** - Intersect polylines with polygons
* **intersect_polyline_with_polygon_files** - Intersect with multiple polygon files

Difference Operations
======================

* **difference_polygon_with_polygon_file** - Calculate polygon difference
* **difference_polyline_with_polygon_file** - Calculate polyline difference with polygon

*************************
Conversion Operations
*************************

Format Conversion
=================

* **convert_vector_format** - Convert between vector formats (shapefile, GeoJSON, etc.)
* **convert_polygon_to_polyline_file** - Convert polygon to polyline features

Raster/Vector Conversion
=========================

* **rasterize_vector** - Convert vector features to raster
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

* **convert_time_series_daily_to_monthly** - Aggregate daily data to monthly

Data Cleaning
=============

* **gap_fill_by_window** - Fill missing data using window-based interpolation

Ocean/Land Utilities
====================

* **define_land_ocean_mask** - Create land/ocean mask from Natural Earth data
* **create_land_ocean_vector_mask_naturalearth** - Generate vector land/ocean mask

Vector Data Operations
======================

* **copy_geometry_without_attributes** - Copy geometries without attribute table

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
   doy = day_of_year(2024, 3, 15)  # Returns 75

   # Check leap year
   is_leap = leap_year(2024)  # Returns True

   # Time code execution
   with timer("Processing"):
       # Your code here
       pass

*************************
Geometry Operations
*************************

Geometric Calculations
======================

* **calculate_hexagon_area** - Calculate area of regular hexagon
* **create_gcs_buffer_zone** - Create buffer zone in geographic coordinates

Buffer Operations
=================

* **create_point_buffer_zone** - Create buffer around point
* **create_polyline_buffer_zone** - Create buffer around polyline
* **create_buffer_zone_polygon_file** - Create buffer for polygon file

.. code-block:: python

   from pyearth.toolbox.geometry.calculate_hexagon_area import calculate_hexagon_area

   # Calculate hexagon area
   area = calculate_hexagon_area(dLength_edge_in=1000.0)  # meters

*************************
Management Operations
*************************

Raster Management
=================

Transformation Operations
-------------------------

* **reproject_raster** - Reproject raster to different coordinate system
* **resample_raster** - Resample raster to different resolution
* **merge_rasters** - Mosaic multiple raster files

.. code-block:: python

   from pyearth.toolbox.management.raster.merge_rasters import merge_rasters

   # Merge multiple raster files
   merge_rasters(
       ['tile1.tif', 'tile2.tif', 'tile3.tif'],
       'merged_output.tif'
   )

Vector Management
=================

File Operations
---------------

* **merge_files** - Merge multiple vector files into one
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

   from pyearth.toolbox.management.vector.merge_files import merge_files

   # Merge multiple shapefiles
   merge_files(
       ['region1.shp', 'region2.shp', 'region3.shp'],
       'merged_regions.shp'
   )

*************************
Math Operations
*************************

Statistical Functions
=====================

* **remap** - Remap value from one range to another
* **scipy_bivariate_kde** - Bivariate kernel density estimation
* **kde_support** - Support functions for KDE calculations

Array Operations
================

* **gap_fill_by_window** - Fill gaps using moving window

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

* **pypoint** - Point class with longitude/latitude
* **pynvector** - Normal vector class for sphere calculations
* **pyline** - Line segment class
* **pypolyline** - Polyline class with multiple segments

Polygon Objects
---------------

* **pypolygon** - Polygon class with vertices
* **pycircle** - Circle class on sphere

.. code-block:: python

   from pyearth.toolbox.mesh.point import pypoint
   from pyearth.toolbox.mesh.polygon import pypolygon

   # Create points
   p1 = pypoint({'dLongitude_degree': -122.0, 'dLatitude_degree': 47.0})
   p2 = pypoint({'dLongitude_degree': -121.0, 'dLatitude_degree': 47.0})
   p3 = pypoint({'dLongitude_degree': -121.0, 'dLatitude_degree': 48.0})

   # Create polygon
   polygon = pypolygon(aPoint=[p1, p2, p3])

Mesh Generation
===============

Grid Generation
---------------

* **create_latlon_mesh** - Create latitude/longitude grid mesh
* **create_square_mesh** - Create square grid mesh
* **create_hexagon_mesh** - Create hexagonal grid mesh

Mesh Algorithms
---------------

* **convert_gcs_coordinates_to_meshcell** - Convert coordinates to mesh cell
* **convert_gcs_coordinates_to_polyline** - Convert coordinates to polyline
* **convert_pcs_coordinates_to_polyline** - Convert projected coordinates to polyline
* **find_minimal_enclosing_polygon** - Find convex hull of points
* **save_points_as_polygon** - Save point array as polygon file
* **split_line_by_length** - Split line into segments of specified length
* **split_polyline_by_length** - Split polyline by length

.. code-block:: python

   from pyearth.toolbox.mesh.latlon.create_latlon_mesh import create_latlon_mesh

   # Create global lat/lon mesh
   cells = create_latlon_mesh(
       dLongitude_left_in=-180.0,
       dLongitude_right_in=180.0,
       dLatitude_bot_in=-90.0,
       dLatitude_top_in=90.0,
       dResolution_in=1.0
   )

*************************
Reader Utilities
*************************

File Reading
============

* **line_count** - Count lines in text file efficiently
* **parse_xml_file** - Parse XML configuration files
* **text_reader_string** - Read text file and return as string

String Operations
=================

* **split_string_into_chunk** - Split string into fixed-size chunks

Configuration Parsing
=====================

* **read_configuration_file** - Read and parse configuration files
* **parse_xml_file_e3sm** - Parse E3SM model XML files (atmosphere and land)

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

   from pyearth.toolbox.analysis.extract.clip_vector_by_polygon_file import \
       clip_vector_by_polygon_file
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

   from pyearth.toolbox.data.convert_time_series_daily_to_monthly import \
       convert_time_series_daily_to_monthly
   from pyearth.toolbox.data.remove_outliers import remove_outliers
   import numpy as np

   # Daily precipitation data (365 days)
   daily_precip = np.random.gamma(2, 2, 365)

   # Remove outliers
   clean_precip = remove_outliers(daily_precip, method='iqr')

   # Convert to monthly
   monthly_precip = convert_time_series_daily_to_monthly(
       clean_precip,
       iYear_start=2024,
       iMonth_start=1,
       iDay_start=1
   )

Example 3: Mesh Generation and Export
======================================

.. code-block:: python

   from pyearth.toolbox.mesh.hexaon.create_hexagon_mesh import create_hexagon_mesh
   from pyearth.toolbox.mesh.algorithm.save_points_as_polygon import \
       save_points_as_polygon

   # Create hexagonal mesh
   cells = create_hexagon_mesh(
       iFlag_rotation_in=0,
       dX_left_in=0.0,
       dX_right_in=100.0,
       dY_bot_in=0.0,
       dY_top_in=100.0,
       dResolution_in=10.0
   )

   # Extract center points and save as polygon
   center_points = [cell.pVertex_center for cell in cells]
   save_points_as_polygon(center_points, 'hexagon_centers.geojson')

Example 4: Merge and Reproject Vectors
=======================================

.. code-block:: python

   from pyearth.toolbox.management.vector.merge_files import merge_files
   from pyearth.toolbox.management.vector.reproject import reproject_vector

   # Merge multiple region files
   merge_files(
       ['region_north.shp', 'region_south.shp', 'region_east.shp'],
       'all_regions.shp'
   )

   # Reproject to UTM
   reproject_vector(
       'all_regions.shp',
       'all_regions_utm.shp',
       'EPSG:32610'  # UTM Zone 10N
   )

Example 5: Data Quality Control
================================

.. code-block:: python

   from pyearth.toolbox.data.remove_outliers import remove_outliers
   from pyearth.toolbox.math.gap_fill_by_window import gap_fill_by_window
   import numpy as np

   # Sensor data with outliers and gaps
   data = np.array([1.0, 2.0, 100.0, 4.0, np.nan, 6.0, 7.0, -50.0, 9.0])

   # Remove outliers
   data_clean = remove_outliers(data, method='zscore', threshold=2.0)

   # Fill gaps
   data_filled = gap_fill_by_window(data_clean, iWindow_size_in=3)
