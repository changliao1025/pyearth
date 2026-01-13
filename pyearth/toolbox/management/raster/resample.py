"""
Raster Resampling Module
=========================

This module provides functionality to resample raster datasets to a specified
resolution and projection using GDAL's warping capabilities. It supports various
resampling algorithms, handles missing data values, preserves multi-band data,
and writes compressed outputs.

The module is intended for geospatial preprocessing where consistent grid
resolution and spatial reference are required across datasets.

Main Function
--------------
resample_raster : Resample a raster file to a target resolution and projection

Features
--------
- Multi-band support (preserve number of bands)
- Flexible resampling algorithms (MODE, NEAREST, BILINEAR, CUBIC, etc.)
- NoData value mapping from source to target
- Optional overwrite and verbose logging
- Outputs compressed GeoTIFF by default (DEFLATE + PREDICTOR=2)

Dependencies
------------
- numpy : Array operations and data type conversion
- osgeo.gdal : Raster I/O, warping, and driver management
- osgeo.osr : Spatial reference system handling
- pyearth.gis.gdal : Custom GDAL utility functions (format support, type mapping)

Typical Workflow
----------------
1. Provide an input raster and desired resolution (x/y)
2. Optionally provide a target projection WKT or use the input raster's
     projection
3. Call `resample_raster` to create the resampled output

Notes
-----
- The function uses GDAL Warp in-memory and writes out a compressed raster.
- If GDAL/NumPy are not available in the environment, static analysis may show
    unresolved import warnings even though the code is correct.

See Also
--------
osgeo.gdal.Warp, osgeo.gdal.WarpOptions
pyearth.gis.gdal.gdal_raster_format_support

Examples
--------
See `resample_raster` docstring for usage examples.

Author
------
PyEarth Development Team

License
-------
Part of the PyEarth package.
"""

import os
import numpy as np
from osgeo import gdal, osr, gdalconst
import logging

from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.gdal.read.raster.gdal_get_raster_spatial_reference import (
    gdal_get_raster_spatial_reference_wkt,
)
from pyearth.gis.gdal.gdal_to_numpy_datatype import gdal_to_numpy_datatype
from pyearth.gis.gdal.gdal_raster_format_support import (
    get_raster_driver_from_filename,
)

gdal.UseExceptions()
# set default spatial reference


def resample_raster(
    sFilename_in,
    sFilename_out,
    dResolution_x,
    dResolution_y,
    pProjection_target_in=None,
    iData_type=gdalconst.GDT_Int16,
    sResampleAlg="MODE",
    dMissing_value_source=255,
    dMissing_value_target=-9999,
    iFlag_overwrite=True,
    iFlag_verbose=False,
):
    """Resample a raster to a target resolution and projection using GDAL Warp.

    This function performs an in-memory GDAL warp followed by writing a
    compressed output raster. It preserves multi-band inputs and supports
    flexible resampling, NoData mapping, and optional overwrite/verbosity.

    Parameters
    ----------
    sFilename_in : str
        Path to the input raster file.
    sFilename_out : str
        Path for the output resampled raster. The output driver is inferred
        from the file extension; GeoTIFF is used as a fallback.
    dResolution_x, dResolution_y : float
        Desired pixel size in X and Y (units are in the target spatial
        reference's coordinate units).
    pProjection_target_in : str, optional
        Target projection WKT or EPSG string (e.g., 'EPSG:4326'). If None,
        the input raster's projection is used.
    iData_type : int, optional
        GDAL data type constant for the output (e.g., ``gdal.GDT_Int16``).
    sResampleAlg : str, optional
        Resampling algorithm to use. Typical values: 'MODE', 'near', 'bilinear',
        'cubic', 'lanczos'. Default 'MODE'.
    dMissing_value_source : numeric, optional
        Value in the source raster that represents missing data.
    dMissing_value_target : numeric, optional
        Value to use as NoData in the output raster.
    iFlag_overwrite : bool, optional
        If True, overwrite an existing output file. Default True.
    iFlag_verbose : bool, optional
        If True, enable informational logging. Default False.

    Returns
    -------
    dict
        Metadata about the created raster: keys include 'success',
        'output_file', 'output_width', 'output_height', 'output_bands',
        'output_projection', 'resampling_algorithm', and 'driver'.

    Notes
    -----
    - Uses an in-memory GDAL Warp (``format='MEM'``) to resample efficiently.
    - Writes output with DEFLATE compression and predictor optimization.
    - Comparisons for floating-point NoData values use ``numpy.isclose`` to
      avoid precision issues.

    Examples
    --------
    Basic resampling to 0.01-degree pixels::

        resample_raster('in.tif', 'out.tif', 0.01, 0.01)

    Resample and reproject to WebMercator with overwrite enabled::

        resample_raster('in.tif', 'out.tif', 30, 30, 'EPSG:3857',
                        sResampleAlg='bilinear', iFlag_overwrite=True)

    """
    logger = logging.getLogger(__name__)
    if iFlag_verbose:
        logging.basicConfig(
            level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
        )

    # defaults and validation
    default_srs = osr.SpatialReference()
    default_srs.ImportFromEPSG(4326)
    sProjection_default = default_srs.ExportToWkt()

    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"Input raster not found: {sFilename_in}")

    # handle output overwrite
    if os.path.exists(sFilename_out):
        if os.path.getsize(sFilename_out) == 0:
            if iFlag_verbose:
                logger.info("Removing empty output file: %s", sFilename_out)
            os.remove(sFilename_out)
        else:
            if iFlag_overwrite:
                if iFlag_verbose:
                    logger.info("Overwriting output file: %s", sFilename_out)
                os.remove(sFilename_out)
            else:
                raise FileExistsError(f"Output exists: {sFilename_out}")
    else:
        if iFlag_verbose:
            logger.info("Creating output file: %s", sFilename_out)

    # open input and get projection
    src_ds = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    if src_ds is None:
        raise RuntimeError(f"Cannot open input raster: {sFilename_in}")

    pProjection_source = gdal_get_raster_spatial_reference_wkt(sFilename_in)
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(pProjection_source)

    pProjection_target = (
        pProjection_target_in
        if pProjection_target_in is not None
        else pProjection_source
    )
    tgt_srs = osr.SpatialReference()
    tgt_srs.ImportFromWkt(pProjection_target)
    if int(gdal.__version__[0]) >= 3:
        src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        tgt_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    # warp in-memory
    warp_opts = gdal.WarpOptions(
        cropToCutline=False,
        xRes=dResolution_x,
        yRes=dResolution_y,
        dstSRS=tgt_srs,
        format="MEM",
        resampleAlg=sResampleAlg,
        srcNodata=dMissing_value_source,  # Source NoData value
        dstNodata=dMissing_value_target,  # Target NoData value
    )
    warped = gdal.Warp("", src_ds, options=warp_opts)
    if warped is None:
        raise RuntimeError("GDAL warp failed")

    gt = warped.GetGeoTransform()
    arr = warped.ReadAsArray()
    if arr is None:
        raise RuntimeError("Failed to read warped data")

    # handle multi-band
    if arr.ndim == 3:
        out_bands = arr.shape[0]
        out_h = arr.shape[1]
        out_w = arr.shape[2]
    else:
        out_bands = 1
        out_h, out_w = arr.shape

    # choose driver from extension
    try:
        driver = get_raster_driver_from_filename(sFilename_out)
    except Exception:
        driver = gdal.GetDriverByName("GTiff")
    if driver is None:
        raise RuntimeError("No suitable GDAL driver found for output")

    creation_options = ["COMPRESS=DEFLATE", "PREDICTOR=2"]
    out_ds = driver.Create(
        sFilename_out,
        out_w,
        out_h,
        out_bands,
        eType=iData_type,
        options=creation_options,
    )
    out_ds.SetGeoTransform(gt)
    out_ds.SetProjection(tgt_srs.ExportToWkt())

    # write bands
    for ib in range(1, out_bands + 1):
        if out_bands == 1:
            band_arr = arr
        else:
            band_arr = arr[ib - 1]

        # cast and nodata handling
        try:
            np_dtype = gdal_to_numpy_datatype(iData_type)
            band_arr = band_arr.astype(np_dtype)
        except Exception:
            pass

        if dMissing_value_source is not None:
            try:
                if np.issubdtype(band_arr.dtype, np.floating):
                    mask = np.isclose(
                        band_arr, dMissing_value_source, rtol=1e-9, atol=1e-9
                    )
                else:
                    mask = band_arr == dMissing_value_source
            except Exception:
                mask = band_arr == dMissing_value_source
            band_arr[mask] = dMissing_value_target
            out_ds.GetRasterBand(ib).SetNoDataValue(dMissing_value_target)

        out_ds.GetRasterBand(ib).WriteArray(band_arr)
        out_ds.GetRasterBand(ib).FlushCache()

    out_ds.FlushCache()
    # cleanup
    warped = None
    out_ds = None
    src_ds = None

    return {
        "success": True,
        "output_file": sFilename_out,
        "output_width": out_w,
        "output_height": out_h,
        "output_bands": out_bands,
        "output_projection": tgt_srs.ExportToWkt(),
        "resampling_algorithm": sResampleAlg,
        "driver": driver.ShortName if hasattr(driver, "ShortName") else None,
    }
