"""ENVI raster writing utilities using GDAL."""

import os
from pathlib import Path
from typing import Any, Optional, Tuple, Union

import numpy as np
from osgeo import gdal

from pyearth.gis.gdal.gdal_to_numpy_datatype import numpy_dtype_to_gdal

PathLike = Union[str, os.PathLike]

__all__ = [
    "gdal_write_envi_file",
    "gdal_write_envi_file_multiple_band",
]


def _resolve_projection_wkt(projection: Any) -> str:
    """Normalize a spatial reference to a WKT string."""
    if projection is None:
        return ""
    if isinstance(projection, str):
        return projection
    if hasattr(projection, "ExportToWkt"):
        return projection.ExportToWkt()
    return str(projection)


def gdal_write_envi_file(
    filename: Optional[PathLike] = None,
    data: Optional[np.ndarray] = None,
    pixel_width: Optional[float] = None,
    origin_x: Optional[float] = None,
    origin_y: Optional[float] = None,
    missing_value: Optional[float] = None,
    spatial_reference: Optional[Any] = None,
    datatype: Optional[int] = None,
    **kwargs,
) -> Tuple[Tuple[float, ...], str]:
    """Write a single-band ENVI format raster file.

    Parameters
    ----------
    filename : str or os.PathLike
        Output file path.
    data : np.ndarray
        2D array of raster data.
    pixel_width : float
        Pixel resolution.
    origin_x : float
        Upper-left X coordinate.
    origin_y : float
        Upper-left Y coordinate.
    missing_value : float
        NoData value.
    spatial_reference : str or osr.SpatialReference
        Spatial reference system.
    datatype : int, optional
        GDAL data type constant. If None, auto-detected from data.dtype.

    Returns
    -------
    tuple
        (geotransform, projection_wkt)
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")
    if data is None:
        data = kwargs.get("aData_in")
    if pixel_width is None:
        pixel_width = kwargs.get("dPixelWidth_in")
    if origin_x is None:
        origin_x = kwargs.get("dOriginX_in")
    if origin_y is None:
        origin_y = kwargs.get("dOriginY_in")
    if missing_value is None:
        missing_value = kwargs.get("dMissing_value_in")
    if spatial_reference is None:
        spatial_reference = kwargs.get("pSpatial_reference_in")

    if filename is None or data is None:
        raise ValueError("Both filename and data must be provided.")

    filepath = Path(filename).resolve()
    filepath.parent.mkdir(parents=True, exist_ok=True)

    if filepath.exists():
        filepath.unlink()

    driver = gdal.GetDriverByName("ENVI")
    if driver is None:
        raise RuntimeError("GDAL ENVI driver is not available.")

    data_arr = np.asarray(data)
    if data_arr.ndim != 2:
        raise ValueError(f"Single-band raster expects 2D array, got {data_arr.ndim}D.")

    nrow, ncolumn = data_arr.shape

    if datatype is None:
        try:
            datatype = numpy_dtype_to_gdal(data_arr.dtype)
        except Exception:
            datatype = gdal.GDT_Float32

    dataset = driver.Create(str(filepath), ncolumn, nrow, 1, datatype)
    if dataset is None:
        raise RuntimeError(f"Unable to create ENVI raster: {filepath}")

    band = None
    try:
        pw = float(pixel_width if pixel_width is not None else 1.0)
        ox = float(origin_x if origin_x is not None else 0.0)
        oy = float(origin_y if origin_y is not None else 0.0)

        geotransform = [ox, pw, 0.0, oy, 0.0, -pw]
        dataset.SetGeoTransform(geotransform)

        proj_wkt = _resolve_projection_wkt(spatial_reference)
        if proj_wkt:
            dataset.SetProjection(proj_wkt)

        band = dataset.GetRasterBand(1)
        band.WriteArray(data_arr)
        if missing_value is not None:
            band.SetNoDataValue(float(missing_value))

        dataset.FlushCache()
        out_gt = dataset.GetGeoTransform()
        out_proj = dataset.GetProjection()

        return out_gt, out_proj
    finally:
        band = None
        dataset = None


def gdal_write_envi_file_multiple_band(
    filename: Optional[PathLike] = None,
    data: Optional[np.ndarray] = None,
    pixel_width: Optional[float] = None,
    origin_x: Optional[float] = None,
    origin_y: Optional[float] = None,
    missing_value: Optional[float] = None,
    spatial_reference: Optional[Any] = None,
    datatype: Optional[int] = None,
    **kwargs,
) -> Tuple[str, Tuple[float, ...], str]:
    """Write a multi-band ENVI format raster file.

    Parameters
    ----------
    filename : str or os.PathLike
        Output file path.
    data : np.ndarray
        3D array of raster data (nband, nrow, ncolumn).
    pixel_width : float
        Pixel resolution.
    origin_x : float
        Upper-left X coordinate.
    origin_y : float
        Upper-left Y coordinate.
    missing_value : float
        NoData value.
    spatial_reference : str or osr.SpatialReference
        Spatial reference system.
    datatype : int, optional
        GDAL data type constant. If None, auto-detected from data.dtype.

    Returns
    -------
    tuple
        (filename, geotransform, projection_wkt)
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")
    if data is None:
        data = kwargs.get("aData_in")
    if pixel_width is None:
        pixel_width = kwargs.get("dPixelWidth_in")
    if origin_x is None:
        origin_x = kwargs.get("dOriginX_in")
    if origin_y is None:
        origin_y = kwargs.get("dOriginY_in")
    if missing_value is None:
        missing_value = kwargs.get("dMissing_value_in")
    if spatial_reference is None:
        spatial_reference = kwargs.get("pSpatial_reference_in")

    if filename is None or data is None:
        raise ValueError("Both filename and data must be provided.")

    filepath = Path(filename).resolve()
    filepath.parent.mkdir(parents=True, exist_ok=True)

    if filepath.exists():
        filepath.unlink()

    driver = gdal.GetDriverByName("ENVI")
    if driver is None:
        raise RuntimeError("GDAL ENVI driver is not available.")

    data_arr = np.asarray(data)
    if data_arr.ndim != 3:
        raise ValueError(f"Multi-band raster expects 3D array, got {data_arr.ndim}D.")

    nband, nrow, ncolumn = data_arr.shape

    if datatype is None:
        try:
            datatype = numpy_dtype_to_gdal(data_arr.dtype)
        except Exception:
            datatype = gdal.GDT_Float32

    dataset = driver.Create(str(filepath), ncolumn, nrow, nband, datatype)
    if dataset is None:
        raise RuntimeError(f"Unable to create multi-band ENVI raster: {filepath}")

    band = None
    try:
        pw = float(pixel_width if pixel_width is not None else 1.0)
        ox = float(origin_x if origin_x is not None else 0.0)
        oy = float(origin_y if origin_y is not None else 0.0)

        geotransform = [ox, pw, 0.0, oy, 0.0, -pw]
        dataset.SetGeoTransform(geotransform)

        proj_wkt = _resolve_projection_wkt(spatial_reference)
        if proj_wkt:
            dataset.SetProjection(proj_wkt)

        for i_band in range(nband):
            band = dataset.GetRasterBand(i_band + 1)
            band.WriteArray(data_arr[i_band, :, :])
            if missing_value is not None:
                band.SetNoDataValue(float(missing_value))

        dataset.FlushCache()
        out_gt = dataset.GetGeoTransform()
        out_proj = dataset.GetProjection()

        return str(filepath), out_gt, out_proj
    finally:
        band = None
        dataset = None
