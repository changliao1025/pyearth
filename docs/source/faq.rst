###########################
Frequently Asked Questions
###########################

Installation Issues
===================

1. Why can't my conda create an environment?
---------------------------------------------

   Turn off your VPN or configure it to bypass conda channels. VPN can interfere with conda's ability to download packages.

2. Why does importing GDAL fail?
---------------------------------

   Consider using the ``conda-forge`` channel:

   .. code-block:: bash

      conda install -c conda-forge gdal

   Make sure GDAL is installed before installing pyearth.

3. PROJ library related issues
-------------------------------

   Reference: https://github.com/OSGeo/gdal/issues/1546

   Make sure you correctly set up the ``PROJ_LIB`` environment variable. The GDAL library depends on the PROJ library, which is often not configured correctly automatically.

   On Linux or macOS, you can set it up in your ``.bash_profile`` or ``.bashrc``:

   **Anaconda:**

   .. code-block:: bash

      export PROJ_LIB=$HOME/.conda/envs/your_env_name/share/proj
      # or
      export PROJ_LIB=$HOME/opt/anaconda3/envs/your_env_name/share/proj

   **Miniconda:**

   .. code-block:: bash

      export PROJ_LIB=/opt/miniconda3/envs/your_env_name/share/proj

Usage Questions
===============

4. What if my model doesn't produce the correct or expected answer?
--------------------------------------------------------------------

   There are several hidden assumptions within the workflow. For example, if you provide DEM and river network data for two different regions, the program won't automatically detect the mismatch.

   **Best practices:**

   * Always perform visual inspection of your input data
   * Verify coordinate systems match across all input files
   * Check that spatial extents overlap appropriately
   * Optionally, enable the ``iFlag_debug`` option in configuration files to output intermediate files for inspection

5. How do I handle International Date Line issues?
---------------------------------------------------

   PyEarth includes utilities for handling geometries that cross the International Date Line:

   .. code-block:: python

      from pyearth.gis.geometry.international_date_line_utility import (
          check_cross_international_date_line_geometry,
          split_international_date_line_polygon_coordinates
      )

6. Which coordinate system should I use?
-----------------------------------------

   * For global analyses, use **WGS84 (EPSG:4326)** with longitude/latitude
   * For regional analyses with metric measurements, use appropriate **UTM zones**
   * PyEarth provides utilities to work with both:

   .. code-block:: python

      from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates
      from pyearth.gis.spatialref.utm_utility import get_utm_zone, get_utm_epsg_code

Package Organization
====================

7. I need visualization functions - where did they go?
-------------------------------------------------------

   Visualization functions have been moved to separate focused packages:

   * **pyearthviz** - For 2D plotting with matplotlib and cartopy
   * **pyearthviz3d** - For 3D globe visualization with GeoVista

   Install them separately:

   .. code-block:: bash

      conda install -c conda-forge pyearthviz pyearthviz3d

8. Where can I find river network analysis functions?
------------------------------------------------------

   River network analysis has been moved to **pyearthriver**:

   .. code-block:: bash

      conda install -c conda-forge pyearthriver

   Or visit: https://github.com/changliao1025/pyearthriver

9. How do I choose between pyearth and pyearthmesh for mesh generation?
------------------------------------------------------------------------

   * **pyearth** - Basic mesh generation (lat/lon, square, hexagon grids)
   * **pyearthmesh** - Advanced mesh generation with more complex geometries and topological operations

   For simple regular grids, use pyearth. For complex mesh workflows, use pyearthmesh.

Performance and Optimization
=============================

10. How can I speed up geometry operations?
--------------------------------------------

    PyEarth includes Cython-optimized functions for performance-critical operations:

    * Geometry calculations use ``pyearth.gis.geometry.kernel``
    * Location operations use ``pyearth.gis.location.kernel``

    These are automatically compiled during installation when Cython is available.

11. How do I handle large raster files?
----------------------------------------

    For large rasters:

    * Use the ``iFlag_metadata_only=1`` parameter to read only metadata without loading data
    * Consider tiling or clipping large files before processing
    * Use GDAL virtual raster (VRT) files for mosaics

    .. code-block:: python

       from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file

       # Read only metadata
       data, geotransform, projection, nrows, ncols, nodata, extent, wkt = \
           gdal_read_geotiff_file('large_file.tif', iFlag_metadata_only=1)

Common Errors
=============

12. "ModuleNotFoundError" when importing pyearth submodules
------------------------------------------------------------

    Make sure you're using the full import path:

    .. code-block:: python

       # Correct
       from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file

       # Incorrect
       from pyearth import gdal_read_geotiff_file

13. Geometry operations return None or unexpected results
----------------------------------------------------------

    Check that:

    * Input geometries are valid (use ``geometry.IsValid()``)
    * Coordinate systems are consistent
    * Longitude values are in expected range (-180 to 180 or 0 to 360)
    * Polygons follow counter-clockwise winding order

Contributing and Support
========================

14. How can I contribute to PyEarth?
-------------------------------------

    Contributions are welcome! Please:

    1. Fork the repository on GitHub
    2. Create a feature branch
    3. Add tests for new functionality
    4. Submit a pull request

    See the :doc:`contribution` page for more details.

15. Where can I get help?
--------------------------

    * Open an issue on GitHub: https://github.com/changliao1025/pyearth/issues
    * Check the documentation: https://pyearth.readthedocs.io
    * Review existing issues for similar problems
