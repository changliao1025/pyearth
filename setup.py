"""Minimal setup.py for building Cython extensions.

This file is used by pyproject.toml to build Cython extensions.
All metadata is now in pyproject.toml, this is only for extension building.
"""
from setuptools import setup
from build_backend import get_extensions

# Build Cython extensions
extensions = get_extensions()

# Minimal setup call - all metadata comes from pyproject.toml
setup(
    ext_modules=extensions,
)
