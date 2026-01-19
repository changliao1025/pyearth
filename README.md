### PyEarth

[![DOI](https://zenodo.org/badge/67889575.svg)](https://zenodo.org/badge/latestdoi/67889575)
[![Downloads](https://static.pepy.tech/badge/pyearth)](https://pepy.tech/project/pyearth)

Every hiker/camper would love to have a leatherman pocket knife around because it is lightweight, versatile and reliable, sometimes life saving.

This is why I developed PyEarth, a lightweight Python package to support various Earth science tasks.
I use it in my daily work, and nearly all my research papers use it in some way.

**Note:** PyEarth has been restructured and is now part of the **PyEarthSuite** ecosystem. To keep each package lightweight and focused, PyEarth has been split into several specialized packages:

- **pyearth** (this package) - Core GIS operations, spatial toolbox, and system utilities
- **pyearthviz** - 2D visualization utilities
- **pyearthviz3d** - 3D visualization with GeoVista
- **pyearthriver** - River network graph algorithms and data structures
- **pyearthmesh** - Mesh generation and manipulation tools (currently not implemented)
- **pyearthhelp** - Helper utilities for data access (NWIS, NLDI, GSIM) and HPC operations

This modular approach allows you to install only what you need while maintaining the ability to use all packages together.

### About PyEarth

PyEarth is designed to be lightweight, so you won't be stopped by Conda install because of dependency issues. You can also clone it and just use the functions you need.

It is versatile, supporting various tasks such as GIS operations, data processing, spatial analysis, and more.

It is designed to be a general-purpose library inspired by the popular IDL Coyote library. Some of the code structure is inspired by the ArcGIS toolbox.

If you find this package useful, please cite it in your work.
You can also support it by [buying me a coffee](https://www.buymeacoffee.com/changliao) or sponsoring it on [GitHub].

### Dependency

PyEarth depends on the following packages:

**Required:**
1. `numpy` - numerical computing
2. `gdal` - geospatial data abstraction library
3. `netCDF4` - netCDF file support
4. `pandas` - data manipulation and analysis
5. `scipy` - scientific computing
6. `rtree` - spatial indexing
7. `geographiclib` - geodesic calculations

**Build-time dependencies:**
- `Cython` (>= 0.29.0) - for building Cython extensions
- `numpy` - required during build for C extensions

### Documentation

Please refer to the [documentation](https://pyearth.readthedocs.io) for details on how to get started using the PyEarth package.

### Installation

`PyEarth` depends on several other packages, including GDAL, which cannot be installed through `pip` easily. You are recommended to use `conda` to install dependencies:

```bash
# Install from conda (recommended)
conda install pyearth

# Or install from pip (requires GDAL pre-installed)
pip install pyearth
```

#### Building from Source

PyEarth now uses modern `pyproject.toml` configuration. To build from source with Cython extensions for maximum performance:

```bash
# Clone the repository
git clone https://github.com/changliao1025/pyearth.git
cd pyearth

# Install dependencies
conda install numpy gdal netcdf4 pandas scipy rtree geographiclib cython

# Build and install (modern way)
pip install -e .

# Or build a distributable package
pip install build
python -m build
```

For Cython compilation, you'll need a C compiler:
- On Linux/macOS: `gcc` or `clang` (usually pre-installed)
- On Windows: Install conda's compilers: `conda install -c conda-forge c-compiler`

**Note:** PyEarth has migrated from `setup.py` to the modern `pyproject.toml` standard. See [MIGRATION_TO_PYPROJECT.md](MIGRATION_TO_PYPROJECT.md) for details on why this is better.

### Content

PyEarth provides general-purpose functions organized into several categories:

#### 1. GIS Module
Provides comprehensive spatial data operations:
- **GDAL wrappers**: Read/write raster (GeoTIFF, ENVI, ASCII) and vector formats
- **Geometry operations**: Distance calculations, polygon area, convex hull checks, angle calculations
- **Location utilities**: Coordinate conversions, spatial indexing, lat/lon to 3D sphere
- **Spatial reference**: Projection handling, coordinate reprojection, UTM utilities
- **ENVI support**: ENVI header file operations

#### 2. Toolbox Module
Collection of utilities for common Earth science tasks:
- **Analysis**: Polygon intersection/difference, vector clipping, attribute filtering
- **Conversion**: Format conversion, rasterization, vectorization
- **Data processing**: Time series conversion, outlier removal, percentile calculations
- **Date utilities**: Julian dates, day-of-year, leap year checks, timers
- **Geometry**: Hexagon calculations, buffer zones
- **Management**: Raster/vector merging, reprojection, resampling
- **Math**: Statistical operations, gap filling, KDE, remapping
- **Mesh**: Generate and manipulate meshes (lat/lon, square, hexagon, circles)
- **Reader**: Configuration files, XML parsing, text processing

#### 3. System Module
System-wide operations:
- Symbolic link creation
- Global variable definitions
- Filename utilities
- Python environment detection

#### 4. Cython Extensions
Performance-critical geometry and location calculations are accelerated with Cython:
- `pyearth.gis.geometry.kernel` - Fast geometry operations
- `pyearth.gis.location.kernel` - Efficient location calculations

### Related Packages in EarthSuite

PyEarth works seamlessly with other EarthSuite packages:



- **[pyearthviz](../pyearthviz)** - 2D plotting and visualization
- **[pyearthviz3d](../pyearthviz3d)** - 3D globe visualization with GeoVista
- **[pyearthriver](../pyearthriver)** - River network topology and graph algorithms
- **[pyearthmesh](../pyearthmesh)** - Advanced mesh generation tools
- **[pyearthhelp](../pyearthhelp)** - Data retrieval and HPC job management

### Acknowledgment

This research was supported as part of the Next Generation Ecosystem Experiments-Tropics, funded by the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research at Pacific Northwest National Laboratory. The study was also partly supported by U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth and Environmental System Modeling program as part of the Energy Exascale Earth System Model (E3SM) project.

### License

Copyright © 2022, Battelle Memorial Institute

1. Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any person or entity lawfully obtaining a copy of this software and associated documentation files (hereinafter "the Software") to redistribute and use the Software in source and binary forms, with or without modification. Such person or entity may use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software and may permit others to do so, subject to the following conditions:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.

* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

* Other than as used herein, neither the name Battelle Memorial Institute or Battelle may be used in any form whatsoever without the express written consent of Battelle.

2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### References

Several publications describe the algorithms used in `PyEarth` in detail. If you make use of `PyEarth` in your work, please consider including a reference to the following:

* Chang Liao. (2022). PyEarth: A lightweight Python package for Earth science (Software). Zenodo. https://doi.org/10.5281/zenodo.6109987

PyEarth is also supporting several other Python packages/projects, including:

* Liao et al., (2023). pyflowline: a mesh-independent river network generator for hydrologic models. Journal of Open Source Software, 8(91), 5446, https://doi.org/10.21105/joss.05446

* Liao. C. (2022). HexWatershed: a mesh independent flow direction model for hydrologic models (0.1.1). Zenodo. https://doi.org/10.5281/zenodo.6425881
