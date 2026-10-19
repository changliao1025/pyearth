"""Raster validity verification utilities using GDAL band checksums."""

import os
from typing import Optional, Union

from osgeo import gdal

PathLike = Union[str, os.PathLike]

__all__ = ["gdal_check_raster_valid"]


def gdal_check_raster_valid(
    filename: Optional[PathLike] = None,
    print_err: bool = False,
    **kwargs,
) -> bool:
    """Check if a raster file is valid by computing checksums for all raster bands.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the raster file to validate.
    print_err : bool, default=False
        If True, print GDAL error messages if validation fails.

    Returns
    -------
    bool
        True if all raster bands can be opened and checksums computed successfully,
        False otherwise.
    """
    if filename is None:
        filename = kwargs.get("sFilename_in")

    if filename is None:
        return False

    filepath = os.fspath(filename)
    if not os.path.exists(filepath):
        if print_err:
            print(f"File does not exist: {filepath}")
        return False

    dataset = None
    try:
        dataset = gdal.Open(filepath, gdal.GA_ReadOnly)
        if dataset is None:
            if print_err:
                error_msg = gdal.GetLastErrorMsg()
                print(
                    f"GDAL error opening file: {error_msg if error_msg else 'Unknown error'}"
                )
            return False

        band_count = dataset.RasterCount
        if band_count < 1:
            if print_err:
                print(f"Raster has no bands: {filepath}")
            return False

        for i in range(band_count):
            band = dataset.GetRasterBand(i + 1)
            if band is None:
                if print_err:
                    print(f"Unable to access band {i + 1}")
                return False
            band.Checksum()
            band = None

        return True

    except Exception as e:
        if print_err:
            print(f"Error validating raster: {e}")
        return False

    finally:
        dataset = None
