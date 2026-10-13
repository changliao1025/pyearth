"""XML file parser for namelist configurations."""

import os
import xml.etree.ElementTree as ET
from typing import Any, Dict
import numpy as np


def parse_xml_file(sFilename_xml_in: str) -> Dict[str, Any]:
    """
    Parse an E3SM/CESM style XML configuration/namelist file.

    Parameters
    ----------
    sFilename_xml_in : str
        Path to the XML file.

    Returns
    -------
    dict
        Parsed namelist configuration dictionary.
    """
    if not os.path.exists(sFilename_xml_in):
        raise FileNotFoundError(f"The XML file does not exist: {sFilename_xml_in}")

    tree = ET.parse(sFilename_xml_in)
    root = tree.getroot()
    namelist: Dict[str, Any] = {}
    keys_str = []

    for entry in root.findall("./group/entry"):
        key = entry.get("id")
        type_elem = entry.find("type")
        dtype = type_elem.text.strip() if type_elem is not None and type_elem.text else ""

        if dtype == "char":
            keys_str.append(key)

        if "value" in entry.attrib:
            val_raw = entry.get("value")
            if dtype == "integer":
                value = int(val_raw)
            elif dtype == "real":
                value = float(val_raw)
            elif dtype == "logical":
                value = val_raw.upper() in ("TRUE", "1", "T")
            else:
                value = val_raw
            namelist[key] = value
        else:
            values = []
            for val_elem in entry.findall("./values"):
                text = val_elem.text.strip() if val_elem.text else ""
                if dtype == "integer":
                    values.append(int(text))
                elif dtype == "real":
                    values.append(float(text))
                elif dtype == "logical":
                    values.append(text.upper() in ("TRUE", "1", "T"))
                else:
                    values.append(text)
            namelist[key] = np.array(values, dtype=np.float64, order="F")

    # Substitute variable references (e.g. $VAR)
    for key in keys_str:
        if key not in namelist:
            continue
        val_sub = str(namelist[key])
        for other_key in keys_str:
            if other_key != key and other_key in namelist:
                val_other = namelist[other_key]
                if isinstance(val_other, str) and f"${key}" in val_other:
                    namelist[other_key] = val_other.replace(f"${key}", val_sub)

    return namelist
