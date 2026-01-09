import numpy as np
from osgeo import gdalconst, ogr
from typing import Any, Type, Union


def gdal_to_numpy_datatype(gdal_dtype: int) -> Type[np.number]:
    """
    Convert a GDAL data type to a NumPy data type.

    Parameters:
    gdal_dtype (int): The GDAL data type (e.g., gdalconst.GDT_Byte).

    Returns:
    Type[np.number]: The corresponding NumPy data type.
    """
    if gdal_dtype == gdalconst.GDT_Byte:
        return np.uint8
    elif gdal_dtype == gdalconst.GDT_UInt16:
        return np.uint16
    elif gdal_dtype == gdalconst.GDT_Int16:
        return np.int16
    elif gdal_dtype == gdalconst.GDT_UInt32:
        return np.uint32
    elif gdal_dtype == gdalconst.GDT_Int32:
        return np.int32
    elif gdal_dtype == gdalconst.GDT_Float32:
        return np.float32
    elif gdal_dtype == gdalconst.GDT_Float64:
        return np.float64
    elif gdal_dtype == gdalconst.GDT_Int8:
        return np.int8
    else:
        raise ValueError(f"GDAL data type {gdal_dtype} not recognized")


def numpy_dtype_to_gdal(numpy_dtype: Type[np.number]) -> int:
    """
    Convert a NumPy data type to a GDAL data type.

    Parameters:
    numpy_dtype (Type[np.number]): The NumPy data type (e.g., np.uint8).

    Returns:
    int: The corresponding GDAL data type.
    """
    if numpy_dtype == np.uint8:
        return gdalconst.GDT_Byte
    elif numpy_dtype == np.uint16:
        return gdalconst.GDT_UInt16
    elif numpy_dtype == np.int16:
        return gdalconst.GDT_Int16
    elif numpy_dtype == np.uint32:
        return gdalconst.GDT_UInt32
    elif numpy_dtype == np.int32:
        return gdalconst.GDT_Int32
    elif numpy_dtype == np.float32:
        return gdalconst.GDT_Float32
    elif numpy_dtype == np.float64:
        return gdalconst.GDT_Float64
    elif numpy_dtype == np.int8:
        return gdalconst.GDT_Int8
    else:
        raise ValueError(f"Numpy data type {numpy_dtype} not recognized")


def numpy_to_gdal_type(
    numpy_value: Any, target_type: int = None
) -> Union[int, float, str, bool, None]:
    """
    Convert a NumPy value to an appropriate type for GDAL/OGR functions.

    Parameters:
    numpy_value (Any): The NumPy value to convert.
    target_type (int, optional): Target OGR field type (e.g., ogr.OFTInteger). Defaults to None.

    Returns:
    Union[int, float, str, bool, None]: The value converted to a Python native type that GDAL/OGR can handle.
    """
    if numpy_value is None or (hasattr(numpy_value, "dtype") and np.isnan(numpy_value)):
        if target_type in (ogr.OFTInteger, ogr.OFTInteger64):
            return 0
        elif target_type == ogr.OFTReal:
            return 0.0
        elif target_type == ogr.OFTString:
            return ""
        else:
            return None

    if target_type is not None:
        if target_type in (ogr.OFTInteger, ogr.OFTInteger64):
            return int(numpy_value)
        elif target_type == ogr.OFTReal:
            return float(numpy_value)
        elif target_type == ogr.OFTString:
            return str(numpy_value)

    if isinstance(
        numpy_value,
        (
            np.integer,
            np.int_,
            np.int8,
            np.int16,
            np.int32,
            np.int64,
            np.uint,
            np.uint8,
            np.uint16,
            np.uint32,
            np.uint64,
        ),
    ):
        return int(numpy_value)
    elif isinstance(numpy_value, (np.float_, np.float16, np.float32, np.float64)):
        return float(numpy_value)
    elif isinstance(numpy_value, np.bool_):
        return bool(numpy_value)
    elif isinstance(numpy_value, (np.string_, np.unicode_)):
        return str(numpy_value)
    else:
        return numpy_value
