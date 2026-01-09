import os
from pathlib import Path


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


def envi_write_header(sFilename_in, aHeader_in):
    """
    Write an ENVI header file from a metadata dictionary.

    The routine currently supports WGS-84 map info and expects the caller to supply
    a dictionary containing ENVI-required fields.

    Parameters
    ----------
    sFilename_in : str
        Path to the header file that will be written (for example ``"foo.hdr"``).
    aHeader_in : dict
        Dictionary containing header metadata. The following keys are required:
        ``sFilename``, ``ncolumn``, ``nrow``, ``nband``, ``offset``, ``data_type``,
        ``byte_order``, ``missing_value``, ``ULlon``, ``ULlat``, ``pixelSize``.
        Provide either ``interleave`` or ``bsq`` to define the interleave layout.
    """

    header_path = Path(sFilename_in)
    header_path.parent.mkdir(parents=True, exist_ok=True)

    missing_keys = sorted(_REQUIRED_KEYS.difference(aHeader_in.keys()))
    if missing_keys:
        missing = ", ".join(missing_keys)
        raise ValueError(f"ENVI header dictionary is missing required keys: {missing}")

    interleave = aHeader_in.get("interleave", aHeader_in.get("bsq"))
    if interleave is None:
        raise ValueError(
            "ENVI header dictionary must include an 'interleave' or 'bsq' entry"
        )

    def _as_str(key):
        try:
            value = aHeader_in[key]
        except KeyError as exc:
            raise ValueError(f"Missing required header value for '{key}'") from exc
        return str(value)

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

    with header_path.open("w", encoding="utf-8") as header_file:
        header_file.write("\n".join(header_lines) + "\n")
