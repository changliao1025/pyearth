#############
System Module
#############

The System module provides system-level configuration, cross-platform file and path operations, dynamic linker management, and Python environment introspection utilities.

***********************************
Global Variables and Constants
***********************************

Physical and Unit Constants
===========================

The ``pyearth.system.define_global_variables`` module defines standard constants used throughout Earth science calculations:

* **earth_radius** - Standard Earth radius in meters (6,378,137.0 m, WGS-84 semi-major axis)
* **missing_value** - Standard missing/NoData fill value (-9999.0)
* **cms2cmd** - Velocity conversion factor: centimeters per second to centimeters per day (86,400)
* **mms2mmd** - Flux conversion factor: millimeters per second to millimeters per day (86,400)
* **feet2meter** - Length conversion: feet to meters (0.3048)
* **inch2mm** - Length conversion: inches to millimeters (25.4)
* **nmonth** - Number of months in a calendar year (12)
* **iMonth_start**, **iMonth_end** - Calendar start (1) and end (12) month indices

Platform and Workspace Information
==================================

* **sPlatform_os** - Operating system identifier (``'linux'``, ``'darwin'``, or ``'windows'``)
* **sMachine** - Hardware machine type (e.g. ``'x86_64'``, ``'arm64'``)
* **sUsername** - Current user login name
* **sWorkspace_home** - Current user home directory path
* **sWorkspace_scratch** - Temporary scratch workspace directory
* **slash** - System path separator (``'/'`` on Unix/macOS, ``'\\'`` on Windows)
* **sPath_proj_lib** - Detected PROJ data directory path
* **sConda_env_name**, **sConda_env_path**, **sConda_env_type** - Active Conda environment metadata

File Extensions
===============

Predefined standard file extensions with dot notation:

* **sExtension_txt** - ``'.txt'``
* **sExtension_envi** - ``'.dat'``
* **sExtension_tiff** - ``'.tif'``
* **sExtension_header** - ``'.hdr'``
* **sExtension_netcdf** - ``'.nc'``
* **sExtension_shapefile** - ``'.shp'``
* **sExtension_json** - ``'.geojson'``
* **sExtension_png** - ``'.png'``
* **sExtension_jpg** - ``'.jpg'``
* **sExtension_ps** - ``'.ps'``
* **sExtension_vtk** - ``'.vtk'``

Environment Summary
===================

* **print_environment_info** - Print formatted summary of the operating system, machine, user, and conda paths.

.. code-block:: python

   from pyearth.system.define_global_variables import (
       earth_radius,
       print_environment_info,
   )

   print(f"Earth radius: {earth_radius} m")
   print_environment_info()

***********************************
Dynamic Linker Utilities
***********************************

Dynamic Library Path Management
===============================

Earth science packages using GDAL and PROJ often encounter shared library loading conflicts when executed across different Linux clusters or HPC systems.

* **configure_linker_path** - Automatically detects Conda or system shared library directories and configures the ``LD_LIBRARY_PATH`` (or ``DYLD_LIBRARY_PATH`` on macOS) environment variable before C extensions or sub-processes are invoked.

.. code-block:: python

   from pyearth.system.set_linker_directory import configure_linker_path

   # Configure linker path for the current workspace
   configure_linker_path('/path/to/workspace')

***********************************
File and Path Utilities
***********************************

Symbolic Links
==============

* **create_symlink** - Safely create symbolic links with support for relative or absolute paths and cross-platform error handling.

.. code-block:: python

   from pyearth.system.create_symlink import create_symlink

   # Create a relative symbolic link
   create_symlink('source_file.tif', 'link_to_source.tif', iFlag_relative_in=1)

Filename and Extension Operations
=================================

* **get_extension_from_filename** - Extract extension from filename (e.g., ``'dem.tif'`` -> ``'.tif'``).
* **get_extension_from_path** - Extract extension from a complete path string.
* **get_filename_without_extension** - Strip extension from a filename (e.g., ``'dem.tif'`` -> ``'dem'``).
* **get_filename_from_path_without_extension** - Extract base filename without extension from a full path.
* **get_folder_path** - Extract parent folder path from a file path.

.. code-block:: python

   from pyearth.system.filename import (
       get_extension_from_filename,
       get_filename_without_extension,
       get_filename_from_path_without_extension,
       get_folder_path,
   )

   path = "/data/gis/watershed.geojson"
   print(get_folder_path(path))                             # "/data/gis"
   print(get_filename_from_path_without_extension(path))    # "watershed"
   print(get_extension_from_filename("watershed.geojson"))  # ".geojson"

***********************************
Python Environment Introspection
***********************************

* **get_python_environment** - Returns Python interpreter executable path, sys.version, and environment information.
* **get_conda_environment** - Inspects and returns the active Conda environment name and prefix.

.. code-block:: python

   from pyearth.system.python.get_python_environment import (
       get_python_environment,
       get_conda_environment,
   )

   py_path, py_ver = get_python_environment()
   conda_name, conda_prefix = get_conda_environment()
   print(f"Running Python {py_ver} in Conda env: {conda_name}")
