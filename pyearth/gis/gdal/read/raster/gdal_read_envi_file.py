"""ENVI raster reading utilities using GDAL."""

import os
from typing import Any, Optional, Tuple, Union

import numpy as np
from osgeo import gdal

gdal.UseExceptions()

PathLike = Union[str, os.PathLike]

__all__ = [
    "RasterDataset",
    "MultiBandRasterDataset",
    "gdal_read_envi_file",
    "gdal_read_envi_file_multiple_band",
]


class RasterDataset:
    """Container for single-band raster data and spatial metadata.

    Supports tuple unpacking:
    (data, pixel_width, origin_x, origin_y, nrow, ncolumn, missing_value, geotransform, projection_wkt)
    as well as dictionary-style and attribute access.
    """

    def __init__(
        self,
        data: Any,
        pixel_width: float,
        origin_x: float,
        origin_y: float,
        nrow: int,
        ncolumn: int,
        missing_value: Optional[float],
        geotransform: Tuple[float, ...],
        projection_wkt: str,
    ) -> None:
        self.data = data
        self.dataOut = data
        self.pixel_width = pixel_width
        self.pixelWidth = pixel_width
        self.origin_x = origin_x
        self.originX = origin_x
        self.origin_y = origin_y
        self.originY = origin_y
        self.nrow = nrow
        self.ncolumn = ncolumn
        self.missing_value = missing_value
        self.missingValue = missing_value
        self.geotransform = geotransform
        self.projection = projection_wkt
        self.projection_wkt = projection_wkt

    def __iter__(self):
        return iter(
            (
                self.data,
                self.pixel_width,
                self.origin_x,
                self.origin_y,
                self.nrow,
                self.ncolumn,
                self.missing_value,
                self.geotransform,
                self.projection_wkt,
            )
        )

    def __getitem__(self, item: Union[int, str]) -> Any:
        if isinstance(item, int):
            return tuple(self)[item]
        if isinstance(item, str) and hasattr(self, item):
            return getattr(self, item)
        raise KeyError(item)

    def __len__(self) -> int:
        return 9


class MultiBandRasterDataset:
    """Container for multi-band raster data and spatial metadata.

    Supports tuple unpacking:
    (data, pixel_width, origin_x, origin_y, nband, nrow, ncolumn, missing_value, geotransform, projection_wkt)
    as well as dictionary-style and attribute access.
    """

    def __init__(
        self,
        data: np.ndarray,
        pixel_width: float,
        origin_x: float,
        origin_y: float,
        nband: int,
        nrow: int,
        ncolumn: int,
        missing_value: Optional[float],
        geotransform: Tuple[float, ...],
        projection_wkt: str,
    ) -> None:
        self.data = data
        self.dataOut = data
        self.pixel_width = pixel_width
        self.pixelWidth = pixel_width
        self.origin_x = origin_x
        self.originX = origin_x
        self.origin_y = origin_y
        self.originY = origin_y
        self.nband = nband
        self.nrow = nrow
        self.ncolumn = ncolumn
        self.missing_value = missing_value
        self.missingValue = missing_value
        self.geotransform = geotransform
        self.projection = projection_wkt
        self.projection_wkt = projection_wkt

    def __iter__(self):
        return iter(
            (
                self.data,
                self.pixel_width,
                self.origin_x,
                self.origin_y,
                self.nband,
                self.nrow,
                self.ncolumn,
                self.missing_value,
                self.geotransform,
                self.projection_wkt,
            )
        )

    def __getitem__(self, item: Union[int, str]) -> Any:
        if isinstance(item, int):
            return tuple(self)[item]
        if isinstance(item, str) and hasattr(self, item):
            return getattr(self, item)
        raise KeyError(item)

    def __len__(self) -> int:
        return 10


def gdal_read_envi_file(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> RasterDataset:
    """Read a single-band ENVI raster and return data with metadata.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the ENVI raster file.

    Returns
    -------
    RasterDataset
        Container unpackable as a 9-tuple:
        (data, pixel_width, origin_x, origin_y, nrow, ncolumn, missing_value, geotransform, projection_wkt)
        and accessible via dictionary keys or attributes.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    RuntimeError
        If GDAL cannot open the file.
    ValueError
        If metadata is missing or invalid.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")

    if filename is None:
        raise ValueError("Filename must be provided.")

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist.")

    dataset = gdal.Open(filepath, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError(f"Unable to open raster {filepath}.")

    band = None
    try:
        geotransform = dataset.GetGeoTransform()
        if not geotransform:
            raise ValueError(
                f"Raster {filepath} does not provide geotransform metadata."
            )

        origin_x = geotransform[0]
        origin_y = geotransform[3]
        pixel_width = geotransform[1]

        ncolumn = dataset.RasterXSize
        nrow = dataset.RasterYSize

        band = dataset.GetRasterBand(1)
        if band is None:
            raise ValueError(f"Raster {filepath} does not contain a raster band.")

        missing_value = band.GetNoDataValue()
        data = band.ReadAsArray(0, 0, ncolumn, nrow)
        if data is None:
            raise RuntimeError(f"Failed to read raster data from {filepath}.")

        projection_wkt = dataset.GetProjectionRef() or dataset.GetProjection() or ""

        return RasterDataset(
            data=data,
            pixel_width=pixel_width,
            origin_x=origin_x,
            origin_y=origin_y,
            nrow=nrow,
            ncolumn=ncolumn,
            missing_value=missing_value,
            geotransform=geotransform,
            projection_wkt=projection_wkt,
        )
    finally:
        band = None
        dataset = None


def gdal_read_envi_file_multiple_band(
    filename: Optional[PathLike] = None,
    **kwargs,
) -> MultiBandRasterDataset:
    """Read a multi-band ENVI raster and return stacked data with metadata.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the ENVI raster file.

    Returns
    -------
    MultiBandRasterDataset
        Container unpackable as a 10-tuple:
        (data, pixel_width, origin_x, origin_y, nband, nrow, ncolumn, missing_value, geotransform, projection_wkt)
        and accessible via dictionary keys or attributes.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")

    if filename is None:
        raise ValueError("Filename must be provided.")

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist.")

    dataset = gdal.Open(filepath, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError(f"Unable to open raster {filepath}.")

    band = None
    try:
        geotransform = dataset.GetGeoTransform()
        if not geotransform:
            raise ValueError(
                f"Raster {filepath} does not provide geotransform metadata."
            )

        origin_x = geotransform[0]
        origin_y = geotransform[3]
        pixel_width = geotransform[1]

        ncolumn = dataset.RasterXSize
        nrow = dataset.RasterYSize
        nband = dataset.RasterCount

        if nband < 1:
            raise ValueError(f"Raster {filepath} does not contain any bands.")

        band = dataset.GetRasterBand(1)
        if band is None:
            raise ValueError(f"Raster {filepath} does not contain raster bands.")

        missing_value = band.GetNoDataValue()
        first_data = band.ReadAsArray(0, 0, ncolumn, nrow)
        if first_data is None:
            raise RuntimeError(f"Failed to read raster data from {filepath}.")

        data = np.empty((nband, nrow, ncolumn), dtype=first_data.dtype)
        data[0, :, :] = first_data

        for i_band in range(1, nband):
            band = dataset.GetRasterBand(i_band + 1)
            if band is None:
                raise ValueError(f"Raster {filepath} is missing band {i_band + 1}.")

            band_data = band.ReadAsArray(0, 0, ncolumn, nrow)
            if band_data is None:
                raise RuntimeError(
                    f"Failed to read raster data from band {i_band + 1} of {filepath}."
                )

            data[i_band, :, :] = band_data

        projection_wkt = dataset.GetProjectionRef() or dataset.GetProjection() or ""

        return MultiBandRasterDataset(
            data=data,
            pixel_width=pixel_width,
            origin_x=origin_x,
            origin_y=origin_y,
            nband=nband,
            nrow=nrow,
            ncolumn=ncolumn,
            missing_value=missing_value,
            geotransform=geotransform,
            projection_wkt=projection_wkt,
        )
    finally:
        band = None
        dataset = None
