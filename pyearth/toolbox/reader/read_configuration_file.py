"""Read text-based configuration file utility."""

import os
from typing import Dict


def read_configuration_file(sFilename_configuration_in: str) -> Dict[str, str]:
    """
    Read a text-based comma-separated key,value configuration file.

    Parameters
    ----------
    sFilename_configuration_in : str
        Path to the configuration file.

    Returns
    -------
    dict
        Dictionary mapping keys to values.
    """
    if not os.path.exists(sFilename_configuration_in):
        raise FileNotFoundError(
            f"The configuration file does not exist: {sFilename_configuration_in}"
        )

    config_out: Dict[str, str] = {}
    with open(sFilename_configuration_in, "r", encoding="utf-8", errors="ignore") as ifs:
        for line in ifs:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            parts = line_str.split(",", 1)
            if len(parts) == 2:
                key = parts[0].strip()
                val = parts[1].strip()
                config_out[key] = val

    return config_out
