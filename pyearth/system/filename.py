"""
File and directory path parsing utilities.

This module provides path extraction, extension inspection, and filename
utilities that handle both string and PathLike objects across platforms.
"""

import os
from typing import Union

PathLike = Union[str, os.PathLike]

__all__ = [
    "get_extension_from_filename",
    "get_extension_from_path",
    "get_filename_from_path_without_extension",
    "get_filename_without_extension",
    "get_folder_path",
]


def get_extension_from_filename(filename: PathLike) -> str:
    """
    Extract the file extension from a given filename.

    Parameters
    ----------
    filename : str or os.PathLike
        The name or path of the file.

    Returns
    -------
    str
        The file extension including the leading dot (e.g., '.txt'),
        or an empty string if no extension is present.
    """
    if filename is None:
        raise ValueError("filename cannot be None")
    _, ext = os.path.splitext(os.fspath(filename))
    return ext


def get_extension_from_path(filepath: PathLike) -> str:
    """
    Extract the file extension from a full file path.

    Parameters
    ----------
    filepath : str or os.PathLike
        The full or relative path to the file.

    Returns
    -------
    str
        The file extension including the leading dot (e.g., '.tif'),
        or an empty string if no extension is present.
    """
    if filepath is None:
        raise ValueError("filepath cannot be None")
    _, ext = os.path.splitext(os.fspath(filepath))
    return ext


def get_filename_from_path_without_extension(filepath: PathLike) -> str:
    """
    Extract the base filename from a file path, excluding its directory and extension.

    Parameters
    ----------
    filepath : str or os.PathLike
        The path to the file.

    Returns
    -------
    str
        The base filename without directory and extension.
    """
    if filepath is None:
        raise ValueError("filepath cannot be None")
    basename = os.path.basename(os.fspath(filepath))
    name, _ = os.path.splitext(basename)
    return name


def get_filename_without_extension(filename: PathLike) -> str:
    """
    Strip the extension from a filename or path.

    Parameters
    ----------
    filename : str or os.PathLike
        The filename or path to process.

    Returns
    -------
    str
        The filename or path stripped of its trailing extension.
    """
    if filename is None:
        raise ValueError("filename cannot be None")
    name, _ = os.path.splitext(os.fspath(filename))
    return name


def get_folder_path(filepath: PathLike) -> str:
    """
    Get the directory containing the file from a file path.

    Parameters
    ----------
    filepath : str or os.PathLike
        The path to the file.

    Returns
    -------
    str
        The parent folder path containing the file.
    """
    if filepath is None:
        raise ValueError("filepath cannot be None")
    return os.path.dirname(os.fspath(filepath))
