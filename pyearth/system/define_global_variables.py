"""Define global variables and constants for the PyEarth package.

This module provides system-wide constants, file extensions, unit conversions,
and environment configuration used throughout the PyEarth package.

Attributes
----------
sPlatform_os : str
    Operating system name ('Windows', 'Linux', 'Darwin')
sUsername : str
    Current username
sWorkspace_home : str
    User's home directory path
sWorkspace_scratch : str
    Scratch workspace directory path
slash : str
    OS-specific path separator
sMachine : str
    Machine type identifier
earth_radius : float
    Earth's radius in meters (WGS84 equatorial radius: 6,378,137 m)

File Extensions
---------------
sExtension_txt : str
    Text file extension (.txt)
sExtension_envi : str
    ENVI data file extension (.dat)
sExtension_tiff : str
    GeoTIFF file extension (.tif)
sExtension_header : str
    Header file extension (.hdr)
sExtension_netcdf : str
    NetCDF file extension (.nc)
sExtension_shapefile : str
    Shapefile extension (.shp)
sExtension_json : str
    JSON file extension (.json)
sExtension_png : str
    PNG image extension (.png)
sExtension_jpg : str
    JPEG image extension (.jpg)
sExtension_ps : str
    PostScript extension (.ps)
sExtension_vtk : str
    VTK file extension (.vtk)

Constants
---------
missing_value : float
    Standard missing/no-data value (-9999.0)
nmonth : int
    Number of months in a year (12)

Unit Conversions
----------------
mms2mmd : float
    Millimeters per second to millimeters per day (86400.0)
feet2meter : float
    Feet to meters conversion factor (0.3048)
inch2mm : float
    Inches to millimeters conversion factor (25.4)
cms2cmd : float
    Centimeters per second to centimeters per day (86400)

Time
----
iMonth_start : int
    First month index (1-based: January = 1)
iMonth_end : int
    Last month index (1-based: December = 12)

Environment
-----------
sConda_env_path : str
    Python environment root directory path
sConda_env_name : str
    Python environment name
sConda_env_type : str
    Python environment type ('conda', 'venv', 'virtualenv', 'pyenv', 'system')

Notes
-----
This module automatically configures PROJ_LIB environment variable by searching
for the PROJ data directory in standard conda/pip installation locations.

The module executes initialization code at import time to:
1. Detect the operating system and set platform-specific paths
2. Detect the Python environment (conda, venv, etc.)
3. Configure PROJ library data path for coordinate transformations

Warnings
--------
Modifying global variables in this module may affect behavior across the entire
PyEarth package. Exercise caution when changing values.
"""

import os
import platform
from pathlib import Path
import getpass
from typing import Optional
from pyearth.system.python.get_python_environment import get_python_environment


# Platform and user information
sPlatform_os = platform.system()
sUsername = getpass.getuser()
sWorkspace_home = str(Path.home())

# Platform-specific configuration
if sPlatform_os == "Windows":
    slash = "\\"
    sMachine = "windows"
    sWorkspace_scratch = "C:\\scratch"
else:  # Linux or Unix-like systems
    slash = "/"

    if sPlatform_os == "Linux":
        sMachine = "linux"
        sWorkspace_scratch = os.path.join(sWorkspace_home, "scratch")
    elif sPlatform_os == "Darwin":
        sMachine = "mac"
        sWorkspace_scratch = os.path.join(sWorkspace_home, "scratch")
    else:
        sMachine = "unix"
        sWorkspace_scratch = os.path.join(sWorkspace_home, "scratch")


# Data file type extensions
sExtension_txt = ".txt"
sExtension_envi = ".dat"
sExtension_tiff = ".tif"
sExtension_header = ".hdr"
sExtension_netcdf = ".nc"
sExtension_shapefile = ".shp"
sExtension_json = ".json"

# Graphics format extensions
sExtension_png = ".png"
sExtension_jpg = ".jpg"
sExtension_ps = ".ps"
sExtension_vtk = ".vtk"


# Constant values
missing_value = -9999.0
nmonth = 12  # Number of months in a year

# Time period defaults
iMonth_start = 1  # January (1-based indexing)
iMonth_end = 12  # December (1-based indexing)

# Unit conversion factors
mms2mmd = 24 * 3600.0  # mm/s to mm/day
feet2meter = 0.3048  # feet to meters
inch2mm = 25.4  # inches to millimeters
cms2cmd = 24 * 3600  # cm/s to cm/day

# Earth parameters (WGS84)
earth_radius = 6378137.0  # Equatorial radius in meters


# Python environment detection
sConda_env_path, sConda_env_name, sConda_env_type = get_python_environment()


def _configure_proj_lib() -> Optional[str]:
    """
    Configure PROJ_LIB environment variable for PROJ library data.

    Searches for the PROJ data directory in standard installation locations
    based on the detected Python environment. Sets the PROJ_LIB environment
    variable if found.

    Returns
    -------
    Optional[str]
        Path to PROJ data directory if found, None otherwise.

    Notes
    -----
    Search order:
    1. <env>/share/proj (Linux/Mac conda/venv)
    2. <env>/Library/share/proj (Windows conda-forge)
    3. <env>/Lib/site-packages/pyproj/proj_dir/share/proj (PyPI pyproj)

    This function is called automatically when the module is imported.
    """
    # Check for PROJ data directory in multiple possible locations
    proj_paths = [
        os.path.join(sConda_env_path, "share", "proj"),  # Linux/Mac typical
        os.path.join(
            sConda_env_path, "Library", "share", "proj"
        ),  # Windows conda-forge
        os.path.join(
            sConda_env_path,
            "Lib",
            "site-packages",
            "pyproj",
            "proj_dir",
            "share",
            "proj",
        ),  # PyPI pyproj
    ]

    for proj_path in proj_paths:
        if os.path.exists(proj_path) and os.path.isdir(proj_path):
            os.environ["PROJ_LIB"] = proj_path
            return proj_path

    return None


# Configure PROJ library on module import
_proj_lib_path = _configure_proj_lib()


def print_environment_info(verbose: bool = True) -> None:
    """
    Print system and environment configuration information.

    Displays platform details, Python environment information, and
    PROJ library configuration status.

    Parameters
    ----------
    verbose : bool, default=True
        If True, prints detailed information including all paths and
        constants. If False, prints only essential environment info.

    Examples
    --------
    >>> from pyearth.system.define_global_variables import print_environment_info
    >>> print_environment_info()
    === PyEarth Environment Configuration ===
    Platform: Linux (linux)
    User: username
    Home: /home/username
    ...

    >>> # Brief output
    >>> print_environment_info(verbose=False)
    Python Environment: conda (myenv)
    PROJ_LIB: /path/to/proj
    """
    print("=== PyEarth Environment Configuration ===")
    print(f"Platform: {sPlatform_os} ({sMachine})")
    print(f"User: {sUsername}")

    if verbose:
        print(f"Home: {sWorkspace_home}")
        print(f"Scratch: {sWorkspace_scratch}")
        print(f"Path separator: '{slash}'")

    print(f"\nPython Environment: {sConda_env_type} ({sConda_env_name})")

    if verbose:
        print(f"Environment Path: {sConda_env_path}")

    if _proj_lib_path:
        print(f"\n✓ PROJ_LIB: {_proj_lib_path}")
    else:
        print("\n⚠ PROJ_LIB: Not found (coordinate transformations may be limited)")

    if verbose:
        print(f"\nConstants:")
        print(f"  Earth radius: {earth_radius:,.0f} m")
        print(f"  Missing value: {missing_value}")
        print(f"\nSupported formats:")
        print(f"  Raster: {sExtension_tiff}, {sExtension_envi}, {sExtension_netcdf}")
        print(f"  Vector: {sExtension_shapefile}, {sExtension_json}")
        print(f"  Image: {sExtension_png}, {sExtension_jpg}")
