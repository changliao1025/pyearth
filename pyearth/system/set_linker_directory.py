"""Configure dynamic linker library path for python environments.

This module sets dynamic linker environment variables (LD_LIBRARY_PATH on Linux,
DYLD_LIBRARY_PATH on macOS, PATH on Windows) to include the python/conda environment's
library directories. This ensures that shared libraries from the environment are
available to Python extensions and compiled modules.

Notes
-----
This module has side effects when imported - it configures dynamic library paths
in os.environ. Import this module early in application startup if you need to ensure
environment libraries are available.

Warnings
--------
Modifying linker paths after program startup may not affect already-loaded
shared libraries. Import this module before importing packages that depend on
compiled extensions (e.g., GDAL, NumPy with MKL, etc.).
"""

import os
import platform
from pathlib import Path
from typing import List, Optional, Union

PathLike = Union[str, os.PathLike]

__all__ = ["configure_linker_path"]


def configure_linker_path(
    env_path: Optional[PathLike] = None,
    additional_paths: Optional[List[PathLike]] = None,
    prepend: bool = True,
    var_name: Optional[str] = None,
) -> str:
    """Configure dynamic linker library path for environment shared libraries.

    Sets LD_LIBRARY_PATH (and DYLD_LIBRARY_PATH on macOS, PATH on Windows) to
    include the environment's existing 'lib' and 'lib64' directories.

    Parameters
    ----------
    env_path : str or os.PathLike, optional
        Path to the conda or virtual environment. If None, automatically detects
        the current environment using get_python_environment().
    additional_paths : list of (str or os.PathLike), optional
        Additional library paths to include.
    prepend : bool, default=True
        If True, prepend detected paths to existing environment variable.
        If False, append them.
    var_name : str, optional
        Specific environment variable to modify. If None, determined
        automatically based on operating system.

    Returns
    -------
    str
        The updated environment variable string that was set.

    Raises
    ------
    ValueError
        If env_path does not exist or is not a directory.
    RuntimeError
        If unable to detect Python environment when env_path is None.

    Examples
    --------
    >>> ld_path = configure_linker_path()
    >>> print('lib' in ld_path)
    True
    """
    # Auto-detect Python environment if not provided
    if env_path is None:
        from pyearth.system.python.get_python_environment import get_python_environment

        try:
            detected_path, _, _ = get_python_environment()
            env_path = detected_path
        except (ValueError, RuntimeError) as e:
            raise RuntimeError(
                "Unable to configure linker path: could not detect Python environment"
            ) from e

    env_path_obj = Path(env_path).resolve()
    if not env_path_obj.exists():
        raise ValueError(f"Environment path does not exist: {env_path}")
    if not env_path_obj.is_dir():
        raise ValueError(f"Environment path is not a directory: {env_path}")

    # Build list of library paths to add - only include existing directories
    candidate_paths = [env_path_obj / "lib", env_path_obj / "lib64"]
    if additional_paths:
        candidate_paths.extend([Path(p) for p in additional_paths])

    new_paths: List[str] = []
    for cand in candidate_paths:
        try:
            cand_resolved = cand.resolve()
            if cand_resolved.is_dir():
                s = str(cand_resolved)
                if s not in new_paths:
                    new_paths.append(s)
        except OSError:
            continue

    # Determine environment variables to update
    os_name = platform.system()
    target_vars: List[str] = []
    if var_name:
        target_vars.append(var_name)
    else:
        if os_name == "Darwin":
            target_vars.extend(["DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"])
        elif os_name == "Windows":
            target_vars.append("PATH")
        else:
            target_vars.append("LD_LIBRARY_PATH")

    result_path = ""
    for target in target_vars:
        existing_val = os.environ.get(target, "").strip()
        existing_list = [p for p in existing_val.split(os.pathsep) if p]

        if prepend:
            merged = [p for p in new_paths if p not in existing_list] + existing_list
        else:
            merged = existing_list + [p for p in new_paths if p not in existing_list]

        combined = os.pathsep.join(merged)
        os.environ[target] = combined
        if not result_path:
            result_path = combined

    return result_path


# Automatically configure when module is imported
try:
    _configured_path = configure_linker_path()
except (ValueError, RuntimeError):
    pass

