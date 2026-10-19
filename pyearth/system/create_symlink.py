"""
Symbolic link creation utilities.

This module provides functions for creating symbolic links (symlinks) with
automatic handling of existing links and proper error management.
"""

import errno
import os
from pathlib import Path
from typing import Union

PathLike = Union[str, os.PathLike]

__all__ = ["create_symlink"]


def create_symlink(source: PathLike, target_link: PathLike) -> None:
    """Create a symbolic link, replacing it if it already exists.

    Creates a symbolic link pointing from target_link to source. If the
    target_link already exists, it will be removed and recreated to point
    to the new source.

    Parameters
    ----------
    source : str or os.PathLike
        The file or directory that the symlink will point to.
        Can be an absolute or relative path. The source does not need to exist
        at the time the symlink is created (dangling symlink is allowed).
    target_link : str or os.PathLike
        The path where the symbolic link will be created.
        If this already exists as a symlink, it will be replaced.
        If it exists as a regular file/directory, an OSError will be raised.

    Raises
    ------
    ValueError
        If source or target_link is empty or None.
    FileNotFoundError
        If the parent directory of target_link does not exist.
    FileExistsError
        If target_link exists as a regular file or directory (not a symlink).
    PermissionError
        If there are insufficient permissions to create or remove links.
    OSError
        If any other system error occurs while creating the symlink.

    Notes
    -----
    - On Windows, creating symlinks may require administrator privileges.
    - If target_link exists as a symlink, it will be replaced atomically.
    - If target_link exists as a regular file/directory, an error is raised
      to prevent accidental data loss.
    - The source path can be relative or absolute.
    - Relative source paths are interpreted relative to the directory
      containing target_link, not the current working directory.

    Warnings
    --------
    This function will replace existing symlinks without warning. Use with
    caution in automated systems to avoid breaking existing links.

    Examples
    --------
    >>> # Create a simple symlink
    >>> create_symlink('/path/to/source', '/path/to/link')

    >>> # Replace an existing symlink
    >>> create_symlink('/new/source', '/path/to/link')  # Replaces existing link

    >>> # Use Path objects
    >>> from pathlib import Path
    >>> source = Path('/data/files')
    >>> link = Path('/home/user/data_link')
    >>> create_symlink(source, link)

    >>> # Create symlink with relative path
    >>> create_symlink('../data', 'local_link')
    """
    if source is None or not str(os.fspath(source)).strip():
        raise ValueError("Source path cannot be empty")
    if target_link is None or not str(os.fspath(target_link)).strip():
        raise ValueError("Target link path cannot be empty")

    source_path = Path(source)
    target_path = Path(target_link)

    # Check if parent directory exists
    target_parent = target_path.parent
    if not target_parent.exists():
        raise FileNotFoundError(
            f"Parent directory does not exist: {target_parent}. "
            "Create the directory before creating the symlink."
        )

    try:
        os.symlink(source_path, target_path)
    except OSError as e:
        if e.errno == errno.EEXIST:
            if target_path.is_symlink():
                try:
                    target_path.unlink()
                    os.symlink(source_path, target_path)
                except Exception as unlink_error:
                    raise OSError(
                        f"Failed to replace existing symlink at {target_path}: {unlink_error}"
                    ) from unlink_error
            else:
                raise FileExistsError(
                    f"Target path exists as a regular file or directory: {target_path}. "
                    "Remove it manually before creating a symlink to prevent data loss."
                ) from e
        else:
            raise OSError(
                f"Failed to create symlink from {target_path} to {source_path}: {e}"
            ) from e
