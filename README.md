### PyEarth

[![DOI](https://zenodo.org/badge/67889575.svg)](https://zenodo.org/badge/latestdoi/67889575)
[![Downloads](https://static.pepy.tech/badge/pyearth)](https://pepy.tech/project/pyearth)

Every hiker/camper would love to have a leatherman pocket knife around because it is lightweight, versatile and reliable, sometimes life saving.

This is why I developed PyEarth, a lightweight Python package to support various Earth science tasks.
I use it in my daily work, and nearly all my research papers use it in some way.

It is supposed to be lightweight, so that you won't be stopped by Conda install because of dependency issues.
You can also clone it and just use the functions you need.

It is supposed to be versatile, so that you can use it for various tasks, such as GIS, data processing, plotting, etc.

It is designed to be a general-purpose library as it is inspired by the popular IDL Coyote library. Some of the code structure is inspired by the ArcGIS toolbox.

If you find this package useful, please cite it in your work.
You can also support it by [buying me a coffee](https://www.buymeacoffee.com/changliao) or sponsoring it on [GitHub].

### Dependency

PyEarth depends on the following packages

1. `numpy`
2. `gdal`
3. `matplotlib`
4. `cartopy`


PyEarth also has optional dependency packages for several functions:

1. `netCDF4` for netCDF support
2. `pandas` for pandas dataframes 
3. `scipy` for scientific computing
4. `statsmodels` for statistical analysis

### Documentation

Please refer to the [documentation](https://pyearth.readthedocs.io) for details on how to get started using the PyEarth package.

### Installation

`PyEarth` depends on several other packages, including gdal, which cannot be installed through `pip` easily. You are recommended to use `conda` to install dependency if necessary.

    conda install pyearth

### Content

PyEarth mainly provides many general-purpose funcations to support other libraries.
These functions are classified into several categories:
1. GIS: This component provides major spatial dataset operations.
2. Toolbox: This component provides many functions for data, date, math, etc.
3. Visual: This component provides a plotting function for time series, scatter, etc.
4. System: This component provides system-wide operations.

You can either call these functions through this package, or you can modify them for your own applications.

### Acknowledgment

This research was supported as part of the Next Generation Ecosystem Experiments-Tropics, funded by the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research at Pacific Northwest National Laboratory. The study was also partly supported by U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth and Environmental System Modeling program as part of the Energy Exascale Earth System Model (E3SM) project. 

### License

Copyright © 2022, Battelle Memorial Institute

1. Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any person or entity lawfully obtaining a copy of this software and associated documentation files (hereinafter “the Software”) to redistribute and use the Software in source and binary forms, with or without modification. Such person or entity may use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software and may permit others to do so, subject to the following conditions:

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


