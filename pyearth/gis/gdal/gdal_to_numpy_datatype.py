"""Type conversion utilities between GDAL raster data types and NumPy data types."""

from typing import Any, Optional, Type, Union

import numpy as np
from osgeo import gdalconst, ogr

__all__ = [
    "gdal_to_numpy_datatype",
    "numpy_dtype_to_gdal",
    "numpy_to_gdal_type",
]

# Mapping GDAL data type constants to NumPy types
_GDAL_TO_NUMPY = {
    gdalconst.GDT_Byte: np.uint8,
    gdalconst.GDT_UInt16: np.uint16,
    gdalconst.GDT_Int16: np.int16,
    gdalconst.GDT_UInt32: np.uint32,
    gdalconst.GDT_Int32: np.int32,
    gdalconst.GDT_Float32: np.float32,
    gdalconst.GDT_Float64: np.float64,
}

if hasattr(gdalconst, "GDT_Int8"):
    _GDAL_TO_NUMPY[gdalconst.GDT_Int8] = np.int8
if hasattr(gdalconst, "GDT_Int64"):
    _GDAL_TO_NUMPY[gdalconst.GDT_Int64] = np.int64
if hasattr(gdalconst, "GDT_UInt64"):
    _GDAL_TO_NUMPY[gdalconst.GDT_UInt64] = np.uint64


def gdal_to_numpy_datatype(gdal_dtype: int) -> Type[Any]:
    """Convert a GDAL data type constant to the corresponding NumPy scalar data type.

    Parameters
    ----------
    gdal_dtype : int
        The GDAL data type constant (e.g. gdalconst.GDT_Byte).

    Returns
    -------
    type
        The corresponding NumPy scalar data type.

    Raises
    ------
    ValueError
        If the GDAL data type is not recognized.
    """
    if gdal_dtype in _GDAL_TO_NUMPY:
        return _GDAL_TO_NUMPY[gdal_dtype]
    raise ValueError(f"GDAL data type {gdal_dtype} not recognized")


def numpy_dtype_to_gdal(numpy_dtype: Union[Type[np.number], np.dtype, str]) -> int:
    """Convert a NumPy data type to the corresponding GDAL data type constant.

    Parameters
    ----------
    numpy_dtype : type, np.dtype, or str
        The NumPy data type (e.g. np.uint8, np.dtype('float32'), or 'int16').

    Returns
    -------
    int
        The corresponding GDAL data type constant (e.g. gdalconst.GDT_Byte).

    Raises
    ------
    ValueError
        If the NumPy data type is not recognized.
    """
    try:
        dt = np.dtype(numpy_dtype)
    except TypeError as e:
        raise ValueError(f"Invalid NumPy data type specification: {numpy_dtype}") from e

    kind = dt.kind
    itemsize = dt.itemsize

    if kind == "u":
        if itemsize == 1:
            return gdalconst.GDT_Byte
        elif itemsize == 2:
            return gdalconst.GDT_UInt16
        elif itemsize == 4:
            return gdalconst.GDT_UInt32
        elif itemsize == 8 and hasattr(gdalconst, "GDT_UInt64"):
            return gdalconst.GDT_UInt64
    elif kind == "i":
        if itemsize == 1 and hasattr(gdalconst, "GDT_Int8"):
            return gdalconst.GDT_Int8
        elif itemsize == 2:
            return gdalconst.GDT_Int16
        elif itemsize == 4:
            return gdalconst.GDT_Int32
        elif itemsize == 8 and hasattr(gdalconst, "GDT_Int64"):
            return gdalconst.GDT_Int64
    elif kind == "f":
        if itemsize == 4:
            return gdalconst.GDT_Float32
        elif itemsize == 8:
            return gdalconst.GDT_Float64
    elif kind == "b":
        return gdalconst.GDT_Byte

    raise ValueError(f"NumPy data type {numpy_dtype} has no matching GDAL type")


def numpy_to_gdal_type(
    numpy_value: Any, target_type: Optional[int] = None
) -> Union[int, float, str, bool, None]:
    """Convert a NumPy scalar or value to a native Python type suitable for GDAL/OGR fields.

    Parameters
    ----------
    numpy_value : Any
        The value to convert.
    target_type : int, optional
        Target OGR field type constant (e.g. ogr.OFTInteger, ogr.OFTReal, ogr.OFTString).

    Returns
    -------
    int, float, str, bool, or None
        The converted Python native value.
    """
    is_nan = False
    if numpy_value is not None:
        try:
            if np.isnan(numpy_value):
                is_nan = True
        except (TypeError, ValueError):
            pass

    if numpy_value is None or is_nan:
        if target_type in (ogr.OFTInteger, ogr.OFTInteger64):
            return 0
        elif target_type == ogr.OFTReal:
            return 0.0
        elif target_type == ogr.OFTString:
            return ""
        return None

    if hasattr(numpy_value, "item"):
        try:
            numpy_value = numpy_value.item()
        except Exception:
            pass

    if target_type is not None:
        if target_type in (ogr.OFTInteger, ogr.OFTInteger64):
            return int(numpy_value)
        elif target_type == ogr.OFTReal:
            return float(numpy_value)
        elif target_type == ogr.OFTString:
            return str(numpy_value)

    if isinstance(numpy_value, (int, np.integer)):
        return int(numpy_value)
    elif isinstance(numpy_value, (float, np.floating)):
        return float(numpy_value)
    elif isinstance(numpy_value, (bool, np.bool_)):
        return bool(numpy_value)
    elif isinstance(numpy_value, (str, np.str_)):
        return str(numpy_value)

    return numpy_value
