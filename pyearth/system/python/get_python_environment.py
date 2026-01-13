"""Get Python environment information for various Python environments.

This module provides utilities to retrieve the Python environment path, name,
and type from the current Python interpreter. Supports conda, venv, virtualenv,
pyenv, and system Python installations.
"""

import sys
import os
from pathlib import Path
from typing import Tuple, Optional


# Type alias for environment types (compatible with Python < 3.8)
# Valid values: 'conda', 'venv', 'virtualenv', 'pyenv', 'system'


def get_python_environment() -> Tuple[str, str, str]:
    """
    Retrieve the Python environment path, name, and type from the interpreter.

    This function detects and extracts information about the current Python
    environment, supporting multiple environment types including conda, venv,
    virtualenv, pyenv, and system Python installations.

    Returns
    -------
    Tuple[str, str, str]
        A tuple containing:
        - env_path : str
            Absolute path to the Python environment directory (or executable
            directory for system Python)
        - env_name : str
            Name of the environment (directory name, or 'system' for system Python)
        - env_type : str
            Type of environment: 'conda', 'venv', 'virtualenv', 'pyenv', or 'system'

    Raises
    ------
    RuntimeError
        If unable to determine Python executable path, or if the detected
        environment path does not exist or is not a directory.

    Notes
    -----
    Environment detection logic:

    1. **Conda**: Checks for 'envs' directory in path or CONDA_DEFAULT_ENV variable
       Path structure: `.../envs/<env_name>/bin/python` or `.../miniconda3/bin/python`

    2. **venv**: Checks for pyvenv.cfg file in parent directories
       Path structure: `<env_path>/bin/python` with pyvenv.cfg in env_path

    3. **virtualenv**: Checks for VIRTUAL_ENV environment variable
       Path structure: `<env_path>/bin/python`

    4. **pyenv**: Checks for 'versions' directory in path or PYENV_VERSION variable
       Path structure: `.../pyenv/versions/<version>/bin/python`

    5. **System**: Default if no virtual environment is detected
       Uses parent directory of Python executable

    Examples
    --------
    >>> # In a conda environment named 'myenv'
    >>> env_path, env_name, env_type = get_python_environment()
    >>> print(f"{env_type}: {env_name}")
    'conda: myenv'
    >>> print(env_path)
    '/home/user/miniconda3/envs/myenv'

    >>> # In a venv environment
    >>> env_path, env_name, env_type = get_python_environment()
    >>> print(env_type)
    'venv'

    >>> # System Python
    >>> env_path, env_name, env_type = get_python_environment()
    >>> print(env_name)
    'system'

    See Also
    --------
    sys.executable : Path to the Python interpreter
    sys.prefix : Python installation prefix
    os.environ : Environment variable dictionary

    Warnings
    --------
    For system Python installations, env_path points to the bin/Scripts directory
    containing the Python executable, not a virtual environment root.
    """
    python_path = sys.executable

    if not python_path:
        raise RuntimeError("Unable to determine Python executable path")

    python_path_obj = Path(python_path).resolve()

    # Check for conda environment
    # Method 1: Check for 'envs' directory in path
    components = str(python_path_obj).split(os.sep)
    if "envs" in components:
        try:
            envs_index = components.index("envs")
            if envs_index + 1 < len(components):
                env_name = components[envs_index + 1]
                env_path = os.sep.join(components[: envs_index + 2])
                env_type = "conda"

                # Validate path
                env_path_obj = Path(env_path)
                if not env_path_obj.exists() or not env_path_obj.is_dir():
                    raise RuntimeError(
                        f"Detected conda environment path '{env_path}' does not exist or is not a directory"
                    )

                return env_path, env_name, env_type
        except (ValueError, IndexError):
            pass  # Fall through to other checks

    # Method 2: Check CONDA_DEFAULT_ENV variable (for conda base or named envs)
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    if conda_env:
        # Use sys.prefix as the environment path for conda
        env_path = sys.prefix
        env_name = conda_env
        env_type = "conda"

        env_path_obj = Path(env_path)
        if env_path_obj.exists() and env_path_obj.is_dir():
            return env_path, env_name, env_type

    # Check for venv (has pyvenv.cfg file)
    # Search up from python executable for pyvenv.cfg
    current = python_path_obj.parent
    for _ in range(3):  # Check up to 3 levels up
        pyvenv_cfg = current / "pyvenv.cfg"
        if pyvenv_cfg.exists():
            env_path = str(current)
            env_name = current.name
            env_type = "venv"
            return env_path, env_name, env_type
        current = current.parent

    # Check for virtualenv (VIRTUAL_ENV variable)
    virtual_env = os.environ.get("VIRTUAL_ENV")
    if virtual_env:
        env_path = virtual_env
        env_name = Path(virtual_env).name
        env_type = "virtualenv"

        env_path_obj = Path(env_path)
        if env_path_obj.exists() and env_path_obj.is_dir():
            return env_path, env_name, env_type

    # Check for pyenv
    # Method 1: Check for 'versions' directory in path
    if "versions" in components:
        try:
            versions_index = components.index("versions")
            if versions_index + 1 < len(components):
                env_name = components[versions_index + 1]
                env_path = os.sep.join(components[: versions_index + 2])
                env_type = "pyenv"

                env_path_obj = Path(env_path)
                if env_path_obj.exists() and env_path_obj.is_dir():
                    return env_path, env_name, env_type
        except (ValueError, IndexError):
            pass

    # Method 2: Check PYENV_VERSION variable
    pyenv_version = os.environ.get("PYENV_VERSION")
    if pyenv_version:
        env_path = sys.prefix
        env_name = pyenv_version
        env_type = "pyenv"

        env_path_obj = Path(env_path)
        if env_path_obj.exists() and env_path_obj.is_dir():
            return env_path, env_name, env_type

    # Default to system Python
    # Use the directory containing the Python executable
    env_path = str(python_path_obj.parent)
    env_name = "system"
    env_type = "system"

    env_path_obj = Path(env_path)
    if not env_path_obj.exists() or not env_path_obj.is_dir():
        raise RuntimeError(
            f"Python executable directory '{env_path}' does not exist or is not a directory"
        )

    return env_path, env_name, env_type


def get_conda_environment() -> Tuple[str, str]:
    """
    Get conda-specific environment information (legacy compatibility).

    This function only works with conda environments and will raise an
    exception for other environment types. For a more flexible solution
    that works with all environment types, use get_python_environment().

    Returns
    -------
    Tuple[str, str]
        A tuple containing:
        - env_path : str
            Absolute path to the conda environment directory
        - env_name : str
            Name of the conda environment

    Raises
    ------
    ValueError
        If the current Python environment is not a conda environment.
    RuntimeError
        If the conda environment path does not exist or is not a directory.

    See Also
    --------
    get_python_environment : Flexible function supporting all environment types

    Examples
    --------
    >>> # In a conda environment
    >>> env_path, env_name = get_conda_environment()
    >>> print(env_name)
    'myenv'
    """
    env_path, env_name, env_type = get_python_environment()

    if env_type != "conda":
        raise ValueError(
            f"Current Python environment is '{env_type}', not 'conda'. "
            f"Use get_python_environment() for non-conda environments."
        )

    return env_path, env_name
