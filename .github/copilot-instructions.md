# PyEarth AI Coding Agent Instructions

## Project Overview
PyEarth is a lightweight Python package for Earth science applications, designed as a "leatherman pocket knife" for geospatial data processing. It provides functionality for GIS operations, data processing, visualization, and system utilities, inspired by IDL Coyote library and ArcGIS toolbox architecture.

## Architecture & Structure

### Four Main Components
1. **`pyearth/gis/`** - Spatial data operations (GDAL/OGR wrappers, coordinate systems, geometry)
2. **`pyearth/toolbox/`** - General utilities (data processing, math, date handling, file management)  
3. **`pyearth/visual/`** - Plotting and visualization (matplotlib/cartopy based)
4. **`pyearth/system/`** - System-wide configuration and platform utilities

### Key Patterns

#### Global Variables Pattern
- Import from `pyearth.system.define_global_variables` for constants, file extensions, units
- Use `missing_value = -9999.0`, `earth_radius`, platform-specific paths (`slash`, `sPlatform_os`)
- File extensions: `sExtension_tiff`, `sExtension_netcdf`, `sExtension_shapefile`, etc.

#### GDAL/OGR Integration
- All spatial functions use GDAL/OGR via `from osgeo import gdal, ogr`
- Type conversion helpers in `gis/gdal/gdal_to_numpy_datatype.py`
- File validation in `gis/gdal/gdal_check_file_type.py` (returns 'raster'/'vector'/'unknown')
- Functions return spatial reference info and handle coordinate transformations

#### Function Documentation Style
```python
def function_name(param: type) -> return_type:
    """Brief description.
    
    Parameters
    ----------
    param : type
        Description
        
    Returns
    -------
    type
        Description
        
    Notes
    -----
    Implementation details, mathematical methods, use cases
    """
```

#### Error Handling
- Validate file existence before GDAL operations
- Handle missing data with `missing_value` constant
- Return None for invalid inputs rather than raising exceptions
- Use type hints extensively

## Development Workflows

### Testing
- Examples in `tests/example/` (minimal test structure)
- Test imports with simple `import pyearth` in test files
- No comprehensive test suite - focus on functional examples

### Package Management
- Uses setuptools with `setup.py` (not pyproject.toml exclusively)
- Core dependencies: `numpy`, `gdal`, `matplotlib`, `cartopy`
- Optional extras: `netCDF4`, `pandas`, `scipy`, `statsmodels` (install with `pip install pyearth[statistics]`)

### Release Process
- Version in `setup.py` must match git tags
- GitHub Actions builds and publishes to PyPI on release
- Manual workflow dispatch supported

## Code Conventions

### Import Organization
```python
import os
from typing import Tuple, List, Optional, Union
import numpy as np
from osgeo import gdal, ogr
from pyearth.system.define_global_variables import *
```

### Function Naming
- Use descriptive names: `gdal_check_file_type`, `get_hydrosheds_continent_from_extent`
- Prefix with module context: `gdal_*`, `plot_*`, `remove_*`

### Data Type Patterns
- Use tuple for bounding boxes: `(minx, miny, maxx, maxy)`
- Handle coordinate system transformations explicitly
- Support both single values and arrays for most functions

### Visual Functions
- All plotting functions save to file (required `sFilename_out` parameter)
- Extensive customization parameters (colors, markers, labels, DPI)
- Import matplotlib as `mpl` and pyplot as `plt`

## Common Tasks

### Adding GIS Functions
1. Place in appropriate `gis/` subdirectory (`gdal/`, `geometry/`, `location/`)
2. Import GDAL/OGR and handle both raster/vector data
3. Use coordinate reference system validation
4. Return structured data with metadata

### Adding Visualization
1. Place in `visual/` subdirectory by plot type
2. Always require output filename parameter
3. Support customization via optional parameters
4. Use global color/style utilities from `visual/color/` and `visual/create_line_style.py`

### Adding Utilities
1. Place in appropriate `toolbox/` subdirectory
2. Focus on reusable, lightweight functions
3. Handle numpy arrays and missing data gracefully
4. Document mathematical methods and use cases extensively

When modifying existing functions, preserve the lightweight philosophy and avoid adding heavy dependencies.