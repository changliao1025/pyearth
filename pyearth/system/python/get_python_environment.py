"""Get Python environment information for various Python environments.

This module provides utilities to retrieve the Python environment path, name,
and type from the current Python interpreter. Supports conda, venv, virtualenv,
pyenv, and system Python installations.
"""

import os
import sys
from pathlib import Path
from typing import Tuple

__all__ = ["get_python_environment", "get_conda_environment"]


def get_python_environment() -> Tuple[str, str, str]:
    """Retrieve the Python environment path, name, and type from the interpreter.

    This function detects and extracts information about the current Python
    environment, supporting multiple environment types including conda, venv,
    virtualenv, pyenv, and system Python installations.

    Returns
    -------
    Tuple[str, str, str]
        A tuple containing:
        - env_path : str
            Absolute path to the Python environment directory (or prefix
            directory for system Python).
        - env_name : str
            Name of the environment (directory name, or 'system' for system Python).
        - env_type : str
            Type of environment: 'conda', 'venv', 'virtualenv', 'pyenv', or 'system'.

    Raises
    ------
    RuntimeError
        If unable to determine Python executable path, or if the detected
        environment path does not exist or is not a directory.
    """
    python_path = sys.executable
    if not python_path:
        raise RuntimeError("Unable to determine Python executable path")

    python_path_obj = Path(python_path).resolve()
    prefix_obj = Path(sys.prefix).resolve()

    # 1. Conda detection
    # Check 1a: conda-meta directory in sys.prefix
    if (prefix_obj / "conda-meta").is_dir():
        env_path = str(prefix_obj)
        if "envs" in prefix_obj.parts:
            env_name = prefix_obj.name
        elif os.environ.get("CONDA_DEFAULT_ENV"):
            env_name = os.environ["CONDA_DEFAULT_ENV"]
        else:
            env_name = prefix_obj.name
        return env_path, env_name, "conda"

    # Check 1b: conda-meta in parent path of python executable
    parts = python_path_obj.parts
    if "envs" in parts:
        try:
            envs_index = parts.index("envs")
            if envs_index + 1 < len(parts):
                candidate_path = Path(*parts[: envs_index + 2])
                if (candidate_path / "conda-meta").is_dir():
                    return str(candidate_path), candidate_path.name, "conda"
        except (ValueError, IndexError):
            pass

    # Check 1c: CONDA_DEFAULT_ENV or CONDA_PREFIX environment variable
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix and (Path(conda_prefix) / "conda-meta").is_dir():
        cand = Path(conda_prefix).resolve()
        return str(cand), conda_env or cand.name, "conda"
    if conda_env and prefix_obj.is_dir():
        return str(prefix_obj), conda_env, "conda"

    # 2. Virtual environment (venv / PEP 405)
    # Check 2a: pyvenv.cfg in sys.prefix or executable parent tree
    for check_dir in [prefix_obj, python_path_obj.parent, python_path_obj.parent.parent]:
        if (check_dir / "pyvenv.cfg").exists():
            return str(check_dir), check_dir.name, "venv"

    # Check 2b: sys.base_prefix != sys.prefix
    base_prefix = getattr(sys, "base_prefix", sys.prefix)
    if prefix_obj != Path(base_prefix).resolve():
        return str(prefix_obj), prefix_obj.name, "venv"

    # 3. Legacy virtualenv
    if hasattr(sys, "real_prefix"):
        return str(prefix_obj), prefix_obj.name, "virtualenv"

    virtual_env = os.environ.get("VIRTUAL_ENV")
    if virtual_env:
        vpath = Path(virtual_env).resolve()
        if vpath.is_dir():
            return str(vpath), vpath.name, "virtualenv"

    # 4. Pyenv detection
    if "versions" in parts and any(".pyenv" in part for part in parts):
        try:
            versions_index = parts.index("versions")
            if versions_index + 1 < len(parts):
                cand = Path(*parts[: versions_index + 2])
                if cand.is_dir():
                    return str(cand), cand.name, "pyenv"
        except (ValueError, IndexError):
            pass

    pyenv_version = os.environ.get("PYENV_VERSION")
    if pyenv_version and prefix_obj.is_dir():
        return str(prefix_obj), pyenv_version, "pyenv"

    # 5. Default: System Python
    env_path = str(prefix_obj) if prefix_obj.is_dir() else str(python_path_obj.parent)
    env_name = "system"
    env_type = "system"

    env_path_obj = Path(env_path)
    if not env_path_obj.exists() or not env_path_obj.is_dir():
        raise RuntimeError(
            f"Python environment directory '{env_path}' does not exist or is not a directory"
        )

    return env_path, env_name, env_type


def get_conda_environment() -> Tuple[str, str]:
    """Get conda-specific environment information.

    Returns
    -------
    Tuple[str, str]
        A tuple containing:
        - env_path : str
            Absolute path to the conda environment directory.
        - env_name : str
            Name of the conda environment.

    Raises
    ------
    ValueError
        If the current Python environment is not a conda environment.
    RuntimeError
        If the conda environment path does not exist or is not a directory.
    """
    env_path, env_name, env_type = get_python_environment()

    if env_type != "conda":
        raise ValueError(
            f"Current Python environment is '{env_type}', not 'conda'. "
            "Use get_python_environment() for non-conda environments."
        )

    return env_path, env_name

