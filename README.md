# PyEarth

[![DOI](https://zenodo.org/badge/67889575.svg)](https://zenodo.org/badge/latestdoi/67889575)

PyEarth is a lightweight Python package to support various Earth science tasks.
It is designed to be a general purpose library as it is inspired by the popular IDL Coyote libarary (http://www.idlcoyote.com/).

Some of the code structure is inspired by the ArcGIS toolbox.

# Content
PyEarth mainly provides many general purpose funcations to support other libaries.
These functions are classified into several categories:
1. GIS: This component provides major spatial dataset operations.
2. Toolbox: This component provides many fuctions for data, date, math, etc.
3. Visual: This component provides plotting function for time series, scatter, etc.
4. System: This component provides system wide operations.

You can either call these functions through this package, or you can modify them for your own applications.

# Acknowledgement
This research was supported as part of the Next Generation Ecosystem Experiments-Tropics, funded by the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research at Pacific Northwest National Laboratory. The study was also partly supported by U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth and Environmental System Modeling program as part of the Energy Exascale Earth System Model (E3SM) project. 

# Install
pyearth depends on several other packages including gdal, which cannot be installed through pip easily. You are recommended to use conda to install dependency if necessary.

1. conda install pyearth

# Contact
Please contact Chang Liao (changliao.climate@gmail.com) if you have any questions.