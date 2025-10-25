from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np


extensions = [
    Extension(
        name="pyearth.gis.geometry.kernel",
        sources=["pyearth/gis/geometry/kernel.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O2"],
    ),
    Extension(
        name="pyearth.gis.location.kernel",
        sources=["pyearth/gis/location/kernel.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O2"],
    )
]

setup(
    name="pyearth-cython-debug",
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
    zip_safe=False,
)
