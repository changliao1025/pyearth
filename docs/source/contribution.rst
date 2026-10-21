============
Contributing
============

Contributions to PyEarth are welcome and greatly appreciated! Every little bit helps, and credit will always be given.

Types of Contributions
=======================

Report Bugs
-----------

Report bugs at https://github.com/changliao1025/pyearth/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting (Python version, GDAL version).
* Detailed steps to reproduce the bug.

Submit Feedback & Feature Requests
-----------------------------------

The best way to send feedback is to file an issue at https://github.com/changliao1025/pyearth/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible to make it easier to implement and maintain.
* Remember that this is an open-source project and contributions are welcome!

Improve Documentation
---------------------

PyEarth could always use more documentation, whether as part of the official PyEarth docs, in docstrings, or even on the web in blog posts and articles.

Write Code
----------

You can contribute bug fixes, new geoprocessing tools, or performance enhancements.

Get Started!
============

Ready to contribute? Here's how to set up `pyearth` for local development:

1. Fork the `pyearth` repository on GitHub.
2. Clone your fork locally:

   .. code-block:: bash

      git clone git@github.com:your_name_here/pyearth.git
      cd pyearth

3. Create a conda environment and install dependencies:

   .. code-block:: bash

      conda create -n earthenv python=3.12 gdal netcdf4 pandas scipy rtree geographiclib cython pytest
      conda activate earthenv

4. Install the package in editable mode with Cython extensions:

   .. code-block:: bash

      pip install -e .

5. Install pre-commit hooks:

   .. code-block:: bash

      pre-commit install

6. Create a branch for local development:

   .. code-block:: bash

      git checkout -b name-of-your-bugfix-or-feature

7. When you're done making changes, check that your changes pass tests and linting:

   .. code-block:: bash

      pytest
      pre-commit run --all-files

8. Commit your changes and push your branch to GitHub:

   .. code-block:: bash

      git add .
      git commit -m "Your detailed description of your changes."
      git push origin name-of-your-bugfix-or-feature

9. Submit a pull request through the GitHub website.

Coding Standards
================

* **Hungarian Variable Conventions**: PyEarth follows explicit variable and parameter naming prefixes:

  - ``sFilename_*``: String file paths
  - ``dLongitude_*``, ``dLatitude_*``: Double-precision floating point coordinates
  - ``aCoordinates_*``, ``aCell_*``: Array-like sequences
  - ``iFlag_*``: Integer configuration flags
  - ``bFlag_*``: Boolean flags
  - ``pSpatialRef_*``, ``pLayer_*``, ``pDataset_*``: OGR/GDAL object references

* **Dual Execution**: Performance-critical loops are implemented in Cython with graceful pure-Python/NumPy fallback.
* **Decoupled Architecture**: PyEarth focuses strictly on core GIS operations, spatial toolbox routines, and system utilities. Visualization belongs in ``pyearthviz`` / ``pyearthviz3d``, and river graph networks belong in ``pyearthriver``.
