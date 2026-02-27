Getting started
===============


Installation
============

PyEarth depends on several packages including GDAL, which cannot be installed through pip easily.
You are recommended to use conda to install dependencies.

Using Conda (Recommended)
--------------------------

.. code-block:: bash

    # Install from conda-forge
    conda install -c conda-forge pyearth

Using pip
---------

If you already have GDAL installed in your system:

.. code-block:: bash

    pip install pyearth

Building from Source
--------------------

To build from source with Cython extensions:

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/changliao1025/pyearth.git
    cd pyearth

    # Install dependencies
    conda install numpy gdal netcdf4 pandas scipy rtree geographiclib cython

    # Build and install
    pip install -e .

For more detailed installation instructions, see the `README <https://github.com/changliao1025/pyearth#installation>`_.

Usage
=====

Basic Import
------------

To use functions from pyearth, you need to import the package or the functions directly:

.. code-block:: python

    import pyearth

    # Import specific GIS functions
    from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
    from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
        calculate_distance_based_on_longitude_latitude
    )

    # Import toolbox functions
    from pyearth.toolbox.date.day_of_year import day_of_year
    from pyearth.toolbox.mesh.polygon import pypolygon

Quick Examples
--------------

Reading a GeoTIFF File
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file

    # Read raster data
    data, geotransform, projection, nrows, ncols, nodata_value, extent, pSRS_wkt = \
        gdal_read_geotiff_file('path/to/file.tif')

Calculating Great Circle Distance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
        calculate_distance_based_on_longitude_latitude
    )

    # Calculate distance between two points (returns distance in meters)
    distance = calculate_distance_based_on_longitude_latitude(
        aLongitude_from=-122.4194,  # Seattle
        aLatitude_from=47.6062,
        aLongitude_to=-73.935242,    # New York
        aLatitude_to=40.730610
    )

Working with Dates
~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pyearth.toolbox.date.day_of_year import day_of_year
    from pyearth.toolbox.date.leap_year import leap_year

    # Get day of year
    doy = day_of_year(2024, 3, 15)  # Returns 75

    # Check if leap year
    is_leap = leap_year(2024)  # Returns True

Creating Mesh Objects
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pyearth.toolbox.mesh.point import pypoint
    from pyearth.toolbox.mesh.polygon import pypolygon

    # Create points
    point1 = pypoint({'dLongitude_degree': -122.4194, 'dLatitude_degree': 47.6062})

    # Create polygon from points
    points = [
        pypoint({'dLongitude_degree': -122.0, 'dLatitude_degree': 47.0}),
        pypoint({'dLongitude_degree': -121.0, 'dLatitude_degree': 47.0}),
        pypoint({'dLongitude_degree': -121.0, 'dLatitude_degree': 48.0}),
        pypoint({'dLongitude_degree': -122.0, 'dLatitude_degree': 48.0}),
    ]
    polygon = pypolygon(aPoint=points)

Module Organization
===================

PyEarth is organized into three main modules:

GIS Module
----------

The **gis** module provides comprehensive spatial data operations:

* **GDAL wrappers**: Read/write raster and vector formats
* **Geometry operations**: Distance calculations, polygon area, convex hull checks
* **Location utilities**: Coordinate conversions, spatial indexing
* **Spatial reference**: Projection handling, coordinate reprojection

Toolbox Module
--------------

The **toolbox** module contains utilities for common Earth science tasks:

* **Analysis**: Polygon operations, vector clipping, attribute filtering
* **Conversion**: Format conversion, rasterization, vectorization
* **Data processing**: Time series conversion, outlier removal
* **Date utilities**: Julian dates, day-of-year, leap year checks
* **Geometry**: Hexagon calculations, buffer zones
* **Management**: Raster/vector merging, reprojection, resampling
* **Math**: Statistical operations, gap filling
* **Mesh**: Generate and manipulate meshes
* **Reader**: Configuration and XML parsing

System Module
-------------

The **system** module provides system-wide operations:

* Symbolic link creation
* Global variable definitions
* Filename utilities
* Python environment detection

Next Steps
==========

* Explore the :doc:`algorithm/gis/gis` for GIS operations
* Check out :doc:`algorithm/toolbox/toolbox` for toolbox functions
* See :doc:`algorithm/api` for complete API reference
* Read the :doc:`faq` for common issues and solutions

Related Packages
================

PyEarth is part of the EarthSuite ecosystem. For specialized tasks, consider:

* **pyearthviz** - 2D visualization with matplotlib/cartopy
* **pyearthviz3d** - 3D globe visualization with GeoVista
* **pyearthriver** - River network topology and graph algorithms
* **pyearthmesh** - Advanced mesh generation tools
* **pyearthhelp** - Data retrieval and HPC job management
