import numpy as np
from osgeo import gdalconst

def gdal_to_numpy_datatype(gdal_dtype):
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

# Usage:
#numpy_dtype = gdal_to_numpy(gdalconst.GDT_Int16)