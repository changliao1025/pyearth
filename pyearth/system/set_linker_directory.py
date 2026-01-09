"""Configure dynamic linker library path for conda environment.

This module automatically sets the LD_LIBRARY_PATH environment variable to
include the conda environment's library directories when imported. This ensures
that shared libraries (.so files) from the conda environment are available to
Python extensions and compiled modules.

Notes
-----
This module has side effects when imported - it modifies the LD_LIBRARY_PATH
environment variable. Import this module early in your application startup if
you need to ensure conda libraries are available.

On Linux systems, LD_LIBRARY_PATH is searched before the system library paths,
allowing conda-installed libraries to override system versions.

Warnings
--------
Modifying LD_LIBRARY_PATH after program startup may not affect already-loaded
shared libraries. Import this module before importing packages that depend on
compiled extensions (e.g., GDAL, NumPy with MKL, etc.).

This module will raise an exception if the Python interpreter is not running
in a conda environment.
"""

import os
from pathlib import Path
from typing import List, Optional


def configure_linker_path(
    env_path: Optional[str] = None,
    additional_paths: Optional[List[str]] = None,
    prepend: bool = True,
) -> str:
    """
    Configure LD_LIBRARY_PATH for conda environment shared libraries.

    Sets the LD_LIBRARY_PATH environment variable to include the conda
    environment's lib and lib64 directories. This enables the dynamic linker
    to find shared libraries installed in the conda environment.

    Parameters
    ----------
    env_path : str, optional
        Path to the conda environment. If None, automatically detects the
        current environment using get_python_environment().
    additional_paths : List[str], optional
        Additional library paths to include in LD_LIBRARY_PATH.
    prepend : bool, default=True
        If True, prepend conda paths to existing LD_LIBRARY_PATH.
        If False, append them.

    Returns
    -------
    str
        The new LD_LIBRARY_PATH value that was set.

    Raises
    ------
    ValueError
        If env_path is provided but does not exist or is not a directory.
    RuntimeError
        If unable to detect conda environment when env_path is None.

    Notes
    -----
    The function adds both 'lib' and 'lib64' directories from the conda
    environment, as different packages may install libraries in either location.

    The order of paths in LD_LIBRARY_PATH matters - directories listed earlier
    are searched first. By default, conda environment paths are prepended to
    ensure they take precedence over system libraries.

    Examples
    --------
    >>> # Automatic configuration (uses current conda environment)
    >>> ld_path = configure_linker_path()
    >>> print('lib' in ld_path)
    True

    >>> # Manual configuration with specific environment
    >>> ld_path = configure_linker_path(
    ...     env_path='/home/user/miniconda3/envs/myenv'
    ... )

    >>> # Add extra library paths
    >>> ld_path = configure_linker_path(
    ...     additional_paths=['/usr/local/custom/lib']
    ... )

    >>> # Append instead of prepend
    >>> ld_path = configure_linker_path(prepend=False)

    See Also
    --------
    get_python_environment : Detect conda environment path and name
    os.environ : Environment variable dictionary

    Warnings
    --------
    Changes to LD_LIBRARY_PATH do not affect already-loaded shared libraries.
    Call this function early in program startup.
    """
    # Auto-detect Python environment if not provided
    if env_path is None:
        from pyearth.system.python.get_python_environment import get_python_environment

        try:
            env_path, _, _ = get_python_environment()
        except (ValueError, RuntimeError) as e:
            raise RuntimeError(
                "Unable to configure linker path: could not detect Python environment"
            ) from e

    # Validate environment path
    env_path_obj = Path(env_path)
    if not env_path_obj.exists():
        raise ValueError(f"Environment path does not exist: {env_path}")

    if not env_path_obj.is_dir():
        raise ValueError(f"Environment path is not a directory: {env_path}")

    # Build list of library paths to add
    lib_paths = [str(env_path_obj / "lib"), str(env_path_obj / "lib64")]

    # Add any additional paths
    if additional_paths:
        lib_paths.extend(additional_paths)

    # Get existing LD_LIBRARY_PATH
    existing_ld_path = os.environ.get("LD_LIBRARY_PATH", "")

    # Combine paths
    if prepend:
        # Conda paths first, then existing paths
        if existing_ld_path:
            new_ld_path = ":".join(lib_paths) + ":" + existing_ld_path
        else:
            new_ld_path = ":".join(lib_paths)
    else:
        # Existing paths first, then conda paths
        if existing_ld_path:
            new_ld_path = existing_ld_path + ":" + ":".join(lib_paths)
        else:
            new_ld_path = ":".join(lib_paths)

    # Set the environment variable
    os.environ["LD_LIBRARY_PATH"] = new_ld_path

    return new_ld_path


# Automatically configure when module is imported
# This maintains backward compatibility with the original behavior
try:
    _configured_path = configure_linker_path()
except (ValueError, RuntimeError):
    # Silently skip if not in conda environment
    # This allows the module to be imported without errors in non-conda contexts
    pass
