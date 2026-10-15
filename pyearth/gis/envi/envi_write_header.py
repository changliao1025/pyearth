"""ENVI raster header generation utilities."""

import os
from pathlib import Path
from typing import Any, Dict, Optional, Union

PathLike = Union[str, os.PathLike]

__all__ = ["envi_write_header"]

_REQUIRED_KEYS = {
    "sFilename",
    "ncolumn",
    "nrow",
    "nband",
    "offset",
    "data_type",
    "byte_order",
    "missing_value",
    "ULlon",
    "ULlat",
    "pixelSize",
}


def envi_write_header(
    filename: Optional[PathLike] = None,
    header_metadata: Optional[Dict[str, Any]] = None,
    **kwargs,
) -> None:
    """Write an ENVI header (.hdr) file from a metadata dictionary.

    The routine supports WGS-84 map info and writes standard ENVI-format
    header text files for binary raster datasets.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the header file to create (e.g. 'raster.hdr').
    header_metadata : dict
        Dictionary containing header metadata. Required keys include:
        'sFilename', 'ncolumn', 'nrow', 'nband', 'offset', 'data_type',
        'byte_order', 'missing_value', 'ULlon', 'ULlat', 'pixelSize',
        and either 'interleave' or 'bsq'.

    Raises
    ------
    ValueError
        If filename is empty, header_metadata is missing required keys,
        or interleave layout is not specified.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")
    if header_metadata is None:
        header_metadata = kwargs.get("aHeader_in")

    if filename is None:
        raise ValueError("Filename must be provided.")
    if header_metadata is None or not isinstance(header_metadata, dict):
        raise ValueError("header_metadata must be a dictionary.")

    header_path = Path(filename).resolve()
    header_path.parent.mkdir(parents=True, exist_ok=True)

    missing_keys = sorted(_REQUIRED_KEYS.difference(header_metadata.keys()))
    if missing_keys:
        missing = ", ".join(missing_keys)
        raise ValueError(f"ENVI header dictionary is missing required keys: {missing}")

    interleave = header_metadata.get("interleave", header_metadata.get("bsq"))
    if interleave is None:
        raise ValueError(
            "ENVI header dictionary must include an 'interleave' or 'bsq' entry"
        )

    def _as_str(key: str) -> str:
        try:
            return str(header_metadata[key])
        except KeyError as exc:
            raise ValueError(f"Missing required header value for '{key}'") from exc

    map_info = (
        "{Geographic Lat/Lon, 1.000, 1.000, "
        f"{_as_str('ULlon')}, {_as_str('ULlat')}, {_as_str('pixelSize')}, {_as_str('pixelSize')}, "
        "WGS-84, units = Degrees}"
    )

    header_lines = [
        "ENVI",
        f"description = {_as_str('sFilename')}",
        f"samples = {_as_str('ncolumn')}",
        f"lines = {_as_str('nrow')}",
        f"bands = {_as_str('nband')}",
        f"header offset = {_as_str('offset')}",
        f"data type = {_as_str('data_type')}",
        f"interleave = {interleave}",
        "sensor type = Unknown",
        f"byte order = {_as_str('byte_order')}",
        f"data ignore value = {_as_str('missing_value')}",
        f"map info = {map_info}",
        "wavelength units = Unknown",
    ]

    with open(header_path, "w", encoding="utf-8") as header_file:
        header_file.write("\n".join(header_lines) + "\n")

