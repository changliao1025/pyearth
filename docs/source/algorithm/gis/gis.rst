#########
GIS/RS
#########

*************************
envi
*************************

* envi_write_header: generate an ENVI header file from a numpy array

*************************
GDAL
*************************

====
Read
====

* gdal read ArcGIS ASCII raster
* gdal read shapefile
* gdal geotiff raster
* gdal read envi raster

=====
Write
=====

* gdal write envi file
* gdal geotiff file

====================
Other gdal functions
====================

* reproject_coordinates: reproject coordinates from one projection to another
* reproject_coordinates_batch: reproject coordinates from one projection to another for a batch of points
* obtain_raster_metadata_geotiff: obtain metadata from a geotiff raster
* obtain_shapefile_metadata: obtain metadata from a shapefile

==========================
Location related functions
==========================

* calculate_distance_based_on_lon_lat: calculate the great circle distance between two points based on their longitude and latitude
* calculate_polygon_area: calcualte the spherical area of a polygon
* convert_360_to_180: convert longitude from 0-360 to -180 to 180