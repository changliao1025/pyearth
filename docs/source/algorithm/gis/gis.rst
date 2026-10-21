############
GIS Module
############

The GIS module provides comprehensive spatial data operations for working with raster and vector data, geometry calculations, coordinate transformations, and spatial reference systems.

*************************
ENVI Format Support
*************************

ENVI Header File Operations
============================

* **envi_write_header** - Generate an ENVI header file from a numpy array

.. code-block:: python

   from pyearth.gis.envi.envi_write_header import envi_write_header

   # Write ENVI header file
   envi_write_header('output.hdr', header_dict)

*************************
GDAL Operations
*************************

Format Support
==============

Raster Formats
--------------

* **gdal_raster_format_support** - Get dictionary of supported raster formats
* **get_raster_format_from_extension** - Determine format from file extension
* **get_raster_format_from_filename** - Get format from complete filename
* **get_raster_driver_from_format** - Get GDAL driver for specific format
* **get_raster_driver_from_filename** - Get GDAL driver from filename

Vector Formats
--------------

* **gdal_vector_format_support** - Get dictionary of supported vector formats
* **get_vector_format_from_extension** - Determine format from file extension
* **get_vector_format_from_filename** - Get format from complete filename
* **get_vector_driver_from_format** - Get GDAL driver for specific format
* **get_vector_driver_from_filename** - Get GDAL driver from filename
* **has_parquet_support** - Check if Parquet/Arrow support is available

File Type Detection
-------------------

* **gdal_check_file_type** - Determine if file is raster, vector, or unknown
* **gdal_check_raster_valid** - Check if raster file is valid
* **gdal_validate_polygon_file** - Validate polygon vector file

Data Type Conversion
--------------------

* **gdal_to_numpy_datatype** - Convert GDAL data type to NumPy dtype
* **numpy_dtype_to_gdal** - Convert NumPy dtype to GDAL data type
* **numpy_to_gdal_type** - Convert NumPy value to appropriate GDAL type

Reading Raster Data
====================

Read Functions
--------------

* **gdal_read_geotiff_file** - Read single-band GeoTIFF raster file
* **gdal_read_geotiff_file_multiple_band** - Read multi-band GeoTIFF file
* **gdal_read_envi_file** - Read single-band ENVI format raster
* **gdal_read_envi_file_multiple_band** - Read multi-band ENVI raster
* **gdal_read_ascii_file** - Read ArcGIS ASCII/AAIGrid raster file

.. code-block:: python

   from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file

   # Read GeoTIFF file
   data, geotransform, projection, nrows, ncols, nodata, extent, wkt = \
       gdal_read_geotiff_file('input.tif')

Metadata Extraction
-------------------

* **gdal_get_raster_extent** - Get bounding box extent of raster
* **gdal_get_raster_spatial_reference_wkt** - Get spatial reference as WKT string

Reading Vector Data
====================

Read Functions
--------------

* **gdal_read_shapefile** - Read shapefile and return geometries
* **gdal_read_vector_file** - General vector file reader

Metadata Extraction
-------------------

* **gdal_get_vector_extent** - Get bounding box extent of vector file
* **gdal_get_vector_boundary** - Get boundary geometry of vector data
* **gdal_get_vector_spatial_reference_wkt** - Get spatial reference as WKT string

Writing Raster Data
====================

Write Functions
---------------

* **gdal_write_geotiff_file** - Write single-band GeoTIFF file
* **gdal_write_geotiff_file_multiple_band** - Write multi-band GeoTIFF file
* **gdal_write_envi_file** - Write single-band ENVI format file
* **gdal_write_envi_file_multiple_band** - Write multi-band ENVI file

.. code-block:: python

   from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import gdal_write_geotiff_file

   # Write GeoTIFF file
   gdal_write_geotiff_file(
       'output.tif',
       data_array,
       geotransform,
       projection_wkt,
       data_type=gdal.GDT_Float32
   )

Writing Vector Data
====================

Write Functions
---------------

* **gdal_write_wkt_to_vector_file** - Write WKT geometry to vector file
* **export_point_to_vector_file** - Export point objects to vector file
* **export_point_as_polygon_file** - Export points as polygon features

*************************
Geometry Operations
*************************

Distance Calculations
=====================

* **calculate_distance_based_on_longitude_latitude** - Great circle distance using haversine formula
* **calculate_distance_to_line** - Calculate distance from point to line on sphere
* **calculate_distance_to_plane** - Calculate distance from point to plane

.. code-block:: python

   from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import \
       calculate_distance_based_on_longitude_latitude

   # Calculate great circle distance
   distance_m = calculate_distance_based_on_longitude_latitude(
       aLongitude_from=-122.4194,
       aLatitude_from=47.6062,
       aLongitude_to=-73.935242,
       aLatitude_to=40.730610
   )

Area Calculations
=================

* **calculate_polygon_area** - Calculate spherical area of polygon
* **calculate_polygon_file_area** - Calculate total area from polygon file
* **calculate_spherical_triangle_area** - Calculate area of spherical triangle
* **spherical_polygon_area** - Alternative spherical polygon area calculation

Angle Calculations
==================

* **calculate_angle_between_point** - Calculate angle between two points
* **calculate_angle_between_point_normal** - Calculate angle with normal vector
* **calculate_angle_between_vectors_degrees** - Angle between 3D vectors in degrees

Polygon Operations
==================

* **check_convex_polygon** - Check if polygon is convex
* **check_counter_clockwise** - Check polygon vertex ordering
* **get_polygon_orientation** - Get orientation as string (CW or CCW)
* **clean_geometry** - Clean and validate geometry objects
* **remove_duplicate_closure** - Remove duplicate closing vertex from polygon coordinate rings

Simplification
--------------

* **douglas_peucker_geodetic** - Simplify polyline using Douglas-Peucker algorithm
* **visvalingam_whyatt_geodetic** - Simplify polyline using Visvalingam-Whyatt algorithm

Geometry Creation
=================

* **create_box_from_longitude_latitude** - Create bounding box geometry
* **extract_unique_vertices_and_connectivity** - Extract unique vertices from mesh

Great Circle Operations
=======================

* **calculate_intersect_on_great_circle** - Find intersection of great circles
* **find_great_circle_intersection_with_meridian** - Find meridian crossing point

International Date Line Handling
=================================

* **IdlHandler** - Centralized IDL processor with configurable strategies
* **IdlStrategy** - Strategy enumeration (CHECK_ONLY, SPLIT, WRAP_TO_360, UNWRAP, DROP, REORDER)
* **IdlResult** - Dataclass encapsulating detection and transformation outputs

Polar Region Handling
=====================

* **PolarHandler** - Centralized detector for polygons enclosing North or South poles via stereographic projection
* **PoleType** - Pole enumeration (NORTH, SOUTH)
* **PolarResult** - Dataclass encapsulating polar enclosure check results and transformed vertices
* **polygon_includes_pole** - Module convenience function to check if coordinates enclose a pole
* **POLAR_TOLERANCE** - Floating point tolerance for point-in-polygon tests (1e-12)
* **LATITUDE_BOUNDARY** - Absolute latitude boundary for poles (90.0 degrees)

.. code-block:: python

   from pyearth.gis.polar_handler import PolarHandler, PoleType, polygon_includes_pole
   import numpy as np

   # Check if polygon encloses South Pole
   coords = np.array([
       [-180.0, -80.0],
       [-90.0, -80.0],
       [0.0, -80.0],
       [90.0, -80.0],
       [180.0, -80.0]
   ])

   handler = PolarHandler(pole=PoleType.SOUTH)
   result = handler.detect(coords)
   print(f"Encloses South Pole: {result.includes_pole}")

   # Or using convenience function
   encloses_south = polygon_includes_pole(coords, pole='south')

Geometry Type Utilities
========================

* **get_output_geometry_type** - Determine appropriate output geometry type for transformed features

Coordinate Range Conversion
============================

* **convert_360_to_180** - Convert longitude from [0, 360] to [-180, 180]
* **convert_180_to_360** - Convert longitude from [-180, 180] to [0, 360]

*************************
Location Utilities
*************************

Coordinate Conversion
=====================

* **convert_longitude_latitude_to_sphere_3d** - Convert lat/lon to 3D Cartesian coordinates
* **convert_sphere_3d_to_longitude_latitude** - Convert 3D Cartesian to lat/lon

Coordinate Extraction
=====================

* **get_geometry_coordinates** - Extract coordinates from OGR geometry
* **get_polygon_exterior_coords** - Get exterior ring coordinates
* **get_multipolygon_exterior_coords** - Get coordinates from multipolygon
* **get_linestring_coords** - Extract linestring coordinates
* **get_point_coords** - Extract point coordinates
* **get_linearring_coords** - Extract linear ring coordinates

Spatial Indexing
================

* **find_index_by_longitude_latitude** - Find array index by coordinates
* **world2Pixel** - Convert world coordinates to pixel coordinates

Regional Classification
=======================

* **get_hydrosheds_continent_from_extent** - Determine HydroSHEDS continent from bounding box

Tile/Pixel Utilities
====================

* **Google_MetersPerPixel** - Calculate meters per pixel at zoom level

*************************
Spatial Reference
*************************

Coordinate Reprojection
=======================

* **reproject_coordinates** - Reproject single point coordinates
* **reproject_coordinates_batch** - Reproject multiple points efficiently

.. code-block:: python

   from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates

   # Reproject from lat/lon to UTM
   x_utm, y_utm = reproject_coordinates(
       dx_in=-122.4194,
       dy_in=47.6062,
       pProjection_source_in='EPSG:4326',
       pProjection_target_in='EPSG:32610'  # UTM Zone 10N
   )

Unit Conversion
===============

* **degree_to_meter** - Convert degree resolution to meters at latitude
* **meter_to_degree** - Convert meter resolution to degrees at latitude

UTM Utilities
=============

* **get_utm_zone** - Calculate UTM zone number from coordinates
* **get_utm_epsg_code** - Get EPSG code for UTM zone
* **get_utm_spatial_reference_wkt** - Create UTM spatial reference WKT

Projection Detection
====================

* **is_wgs84_projection** - Check if projection is WGS84 (EPSG:4326)

*************************
Usage Examples
*************************

Example 1: Read and Reproject Raster
=====================================

.. code-block:: python

   from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
   from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import gdal_write_geotiff_file
   from pyearth.toolbox.management.raster.reproject import reproject_raster

   # Read raster
   data, geotransform, projection, nrows, ncols, nodata, extent, wkt = \
       gdal_read_geotiff_file('input.tif')

   # Reproject to different coordinate system
   reproject_raster('input.tif', 'output_reprojected.tif', 'EPSG:32610')

Example 2: Calculate Polygon Area
==================================

.. code-block:: python

   from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
   import numpy as np

   # Define polygon vertices (longitude, latitude)
   lon = np.array([-122.0, -121.0, -121.0, -122.0, -122.0])
   lat = np.array([47.0, 47.0, 48.0, 48.0, 47.0])

   # Calculate area in square meters
   area_m2 = calculate_polygon_area(lon, lat)

Example 3: Handle International Date Line
==========================================

.. code-block:: python

   from pyearth.gis.idl_handler import IdlHandler, IdlStrategy
   from osgeo import ogr
   import numpy as np

   handler = IdlHandler(strategy=IdlStrategy.SPLIT)

   # Fix an OGR geometry crossing the IDL
   geometry = ogr.CreateGeometryFromWkt('POLYGON((170 10, -170 10, -170 -10, 170 -10, 170 10))')
   fixed_geometries = handler.fix(geometry)

   # Or detect and split polygon coordinates
   coords = np.array([[170, 10], [-170, 10], [-170, -10], [170, -10], [170, 10]])
   result = handler.detect(coords)
   if result.crosses_idl:
       eastern, western = result.sub_polygons

Example 4: Work with Vector Files
==================================

.. code-block:: python

   from pyearth.gis.gdal.read.vector.gdal_read_vector_file import gdal_read_shapefile
   from pyearth.gis.gdal.write.vector.gdal_write_wkt_to_vector_file import \
       gdal_write_wkt_to_vector_file

   # Read shapefile
   geometries = gdal_read_shapefile('input.shp')

   # Write first geometry to new file
   if geometries:
       wkt = geometries[0].ExportToWkt()
       gdal_write_wkt_to_vector_file(wkt, 'output.geojson', layer_name='features')

Example 5: Check Polar Region Enclosure
=======================================

.. code-block:: python

   from pyearth.gis.polar_handler import PolarHandler, PoleType
   import numpy as np

   # Define polygon vertices circling Antarctica
   antarctica_ring = np.array([
       [-180.0, -70.0],
       [-90.0, -70.0],
       [0.0, -70.0],
       [90.0, -70.0],
       [180.0, -70.0]
   ])

   handler = PolarHandler(pole=PoleType.SOUTH)
   result = handler.detect(antarctica_ring)

   if result.includes_pole:
       print("The polygon encloses the South Pole.")
