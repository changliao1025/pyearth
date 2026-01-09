import os
import logging
from typing import Union
from pathlib import Path
from osgeo import gdal


def gdal_check_raster_valid(
    sFilename_in: Union[str, Path], print_err: bool = False
) -> bool:
    """Check if a raster file is valid by computing band checksums.

    This function attempts to open a raster and compute checksums for all bands
    as a quick validity test. While not infallible, it should capture most cases
    of corrupted or invalid raster files.

    Parameters
    ----------
    sFilename_in : str or Path
        Path to the raster file to validate.
    print_err : bool, optional
        If True, print GDAL error messages to stdout. Default is False.

    Returns
    -------
    bool
        True if the raster is valid (all band checksums computed successfully),
        False otherwise.

    Notes
    -----
    Adapted from https://lists.osgeo.org/pipermail/gdal-dev/2013-November/037520.html
    This check reads through all raster bands to compute checksums, which may be
    slow for large files. Returns False if the file doesn't exist or can't be opened.
    """

    sFilename_in = str(sFilename_in)

    if not os.path.exists(sFilename_in):
        if print_err:
            print(f"File does not exist: {sFilename_in}")
        return False

    dataset = None
    try:
        dataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
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
                print(f"Raster has no bands: {sFilename_in}")
            return False

        for i in range(band_count):
            band = dataset.GetRasterBand(i + 1)
            if band is None:
                if print_err:
                    print(f"Unable to access band {i + 1}")
                return False
            band.Checksum()

        return True

    except RuntimeError as e:
        if print_err:
            error_msg = gdal.GetLastErrorMsg()
            print(f"GDAL error: {error_msg if error_msg else str(e)}")
        return False

    except Exception as e:
        if print_err:
            print(f"Unexpected error validating raster: {str(e)}")
        return False

    finally:
        dataset = None
