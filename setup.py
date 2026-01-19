import io
import os

from setuptools import setup, find_packages, Extension

import numpy

# Try to import Cython, but don't fail if it's not available
try:
    from Cython.Build import cythonize
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

NAME = "pyearth"
DESCRIPTION = "Python for Earth Science."
AUTHOR = "Chang Liao"
AUTHOR_EMAIL = "changliao.climate@gmail.com"
URL = "https://github.com/changliao1025/pyearth"
VERSION = "0.2.1"
REQUIRES_PYTHON = ">=3.9.0"
KEYWORDS = "Earth Science"

REQUIRED = [
    "numpy",
    "gdal",
    "netCDF4",
    "pandas",
    "scipy",
    "rtree",
    "geographiclib",
]

CLASSIFY = [
    "Development Status :: 4 - Beta",
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Visualization",
]

HERE = os.path.abspath(os.path.dirname(__file__))


# Cython extensions (build from .pyx if Cython available, otherwise from .c files)
if HAVE_CYTHON:
    ext_sources = {
        "pyearth.gis.geometry.kernel": ["pyearth/gis/geometry/kernel.pyx"],
        "pyearth.gis.location.kernel": ["pyearth/gis/location/kernel.pyx"],
    }
else:
    ext_sources = {
        "pyearth.gis.geometry.kernel": ["pyearth/gis/geometry/kernel.c"],
        "pyearth.gis.location.kernel": ["pyearth/gis/location/kernel.c"],
    }

# Build extensions list; NumPy is required for header includes
extensions = [
    Extension(
        name,
        sources,
        include_dirs=[numpy.get_include()],
        libraries=[],
        library_dirs=[],
    )
    for name, sources in ext_sources.items()
]

try:
    with io.open(os.path.join(HERE, "README.md"), encoding="utf-8") as f:
        LONG_DESCRIPTION = "\n" + f.read()

except FileNotFoundError:
    LONG_DESCRIPTION = DESCRIPTION


setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    license="BSD-3-Clause",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    python_requires=REQUIRES_PYTHON,
    keywords=KEYWORDS,
    url=URL,
    packages=find_packages(),
    install_requires=REQUIRED,
    include_package_data=True,
    classifiers=CLASSIFY,
    ext_modules=cythonize(extensions, compiler_directives={'language_level': "3"}) if HAVE_CYTHON else extensions,
    setup_requires=["numpy", "Cython>=0.29.0"],
    extras_require={

    },
)
