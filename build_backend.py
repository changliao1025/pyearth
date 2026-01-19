"""Custom build configuration to handle Cython extensions.

This module provides a function to build Cython extensions that can be used
with setuptools. It tries to use .pyx files if Cython is available,
otherwise falls back to pre-compiled .c files.
"""
from setuptools import Extension


def get_extensions():
    """Build list of extension modules with proper numpy include paths.

    Returns:
        list: List of Extension objects for Cython modules
    """
    import numpy

    # Try to import Cython, but don't fail if it's not available
    try:
        from Cython.Build import cythonize
        HAVE_CYTHON = True
    except ImportError:
        HAVE_CYTHON = False
        cythonize = None

    # Define Cython extension sources
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

    # Build extensions list
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

    # Cythonize if Cython is available
    if HAVE_CYTHON:
        extensions = cythonize(
            extensions,
            compiler_directives={'language_level': "3"},
            force=False  # Only rebuild if source changed
        )

    return extensions
