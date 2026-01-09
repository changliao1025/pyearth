"""
Raster Reprojection Module
==========================

This module provides helpers to reproject raster datasets between spatial
reference systems using GDAL's warping functionality. It uses an in-memory warp
to perform coordinate transformations and resampling, and writes compressed
output files.

Main Function
--------------
reproject_raster : Reproject a raster file to a target SRS using GDAL Warp

Features
--------
- Multi-band support (preserves band count)
- Optional forced output resolution (xRes/yRes)
- Flexible resampling algorithm selection (near, bilinear, cubic, mode, etc.)
- NoData value mapping and min/max clamping
- Compressed GeoTIFF output by default (DEFLATE + PREDICTOR=2)

Dependencies
------------
- numpy : Array operations and data type conversion
- osgeo.gdal : Raster I/O and warp operations
- osgeo.osr : Spatial reference handling
- pyearth.gis.gdal : Utility helpers (driver detection, type mapping)

Notes
-----
- The module uses GDAL Warp with ``format='MEM'`` for intermediate processing.
- Static linters may report unresolved optional imports if GDAL/NumPy are not
    available in the analysis environment.

See Also
--------
osgeo.gdal.Warp, osgeo.gdal.WarpOptions
pyearth.gis.gdal.gdal_raster_format_support

Examples
--------
        reproject_raster('input.tif', 'out.tif', 'EPSG:3857')

Author
------
PyEarth Development Team

License
-------
Part of the PyEarth package.
"""

import os
import numpy as np
import osgeo
import logging
from osgeo import gdal, osr
from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.gdal.gdal_to_numpy_datatype import gdal_to_numpy_datatype
from pyearth.gis.gdal.gdal_raster_format_support import (
    get_raster_driver_from_filename,
)


def reproject_raster(
    sFilename_raster_in,
    sFilename_raster_out,
    pProjection_target,
    xRes=None,
    yRes=None,
    sResampleAlg="near",
    iFlag_force_resolution_in=0,
    iData_type=None,
    iFlag_overwrite=True,
    iFlag_verbose=False,
):
    """Reproject a raster to a target spatial reference using GDAL Warp.

    This function is a robust wrapper around `osgeo.gdal.Warp`. It performs an
    in-memory warp and writes the result to disk with compression. The function
    preserves the number of bands from the warped result, maps NoData values
    (with float/int-aware comparisons), optionally forces output resolution,
    and allows specifying an output GDAL data type.

    Parameters
    ----------
    sFilename_raster_in : str
        Path to the input raster file.
    sFilename_raster_out : str
        Path where the reprojected raster will be written. The GDAL driver is
        inferred from the file extension; GeoTIFF is used as a fallback.
    pProjection_target : str
        Target projection as a WKT string or an EPSG code string like
        ``'EPSG:3857'``.
    xRes, yRes : float, optional
        Optional forced pixel resolution for output (units depend on target SRS).
        Used only when ``iFlag_force_resolution_in==1``.
    sResampleAlg : str, optional
        Resampling algorithm name passed to GDAL (e.g. 'near', 'bilinear',
        'cubic', 'mode'). Default is 'near'.
    iFlag_force_resolution_in : int, optional
        Set to 1 to apply the provided ``xRes``/``yRes``. Default 0.
    iData_type : int, optional
        GDAL data type constant (e.g. ``gdal.GDT_Byte``, ``gdal.GDT_Int16``).
        If ``None`` the source dataset first band's type is used.
    iFlag_overwrite : bool, optional
        If True, overwrite existing output file. Default True.
    iFlag_verbose : bool, optional
        If True, enable informational logging to stdout.

    Returns
    -------
    dict
        Metadata about the output including keys: 'success', 'output_file',
        'output_width', 'output_height', 'output_bands', 'output_projection',
        'resampling_algorithm', and 'driver'.

    Notes
    -----
    - The function uses an in-memory Warp (``format='MEM'``) for performance
      and then writes the final file with DEFLATE compression and predictor
      optimization.
    - NoData value comparisons use ``numpy.isclose`` for floating point arrays
      to avoid precision issues.
    - If GDAL or NumPy are not available in the execution environment, the
      module will fail to import or run; static linters may also complain about
      unresolved imports even if the runtime has the dependencies available.

    Examples
    --------
    Basic usage to reproject a file to WebMercator::

        reproject_raster('input.tif', 'out.tif', 'EPSG:3857')

    Forcing resolution (in target SRS units)::

        reproject_raster('in.tif', 'out.tif', 'EPSG:3857', xRes=30, yRes=30,
                         iFlag_force_resolution_in=1, sResampleAlg='bilinear')

    """
    logger = logging.getLogger(__name__)
    if iFlag_verbose:
        logging.basicConfig(
            level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
        )

    # validation
    if not isinstance(sFilename_raster_in, str) or not os.path.exists(
        sFilename_raster_in
    ):
        raise FileNotFoundError(f"Input raster not found: {sFilename_raster_in}")
    if not isinstance(sFilename_raster_out, str):
        raise TypeError("sFilename_raster_out must be a string")
    if not isinstance(pProjection_target, str):
        raise TypeError("pProjection_target must be a WKT or EPSG string")

    # handle overwrite
    if os.path.exists(sFilename_raster_out):
        if iFlag_overwrite:
            try:
                os.remove(sFilename_raster_out)
            except Exception as e:
                raise OSError(f"Failed to remove existing output: {e}")
        else:
            raise FileExistsError(f"Output exists: {sFilename_raster_out}")

    # open source
    src_ds = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)
    if src_ds is None:
        raise RuntimeError(f"GDAL cannot open input raster: {sFilename_raster_in}")

    src_band1 = src_ds.GetRasterBand(1)
    src_nodata = src_band1.GetNoDataValue()
    src_min = src_band1.GetMinimum()
    src_max = src_band1.GetMaximum()
    src_dtype = src_band1.DataType
    band_count = src_ds.RasterCount or 1

    if iData_type is None:
        iData_type_out = src_dtype
    else:
        iData_type_out = iData_type

    # prepare target SRS
    target_srs = osr.SpatialReference()
    if isinstance(pProjection_target, str) and pProjection_target.upper().startswith(
        "EPSG:"
    ):
        target_srs.ImportFromEPSG(int(pProjection_target.split(":", 1)[1]))
    else:
        target_srs.ImportFromWkt(pProjection_target)
    if int(osgeo.__version__[0]) >= 3:
        target_srs.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    # build warp options
    resample = sResampleAlg if isinstance(sResampleAlg, str) else "near"
    warp_kwargs = dict(
        dstSRS=target_srs.ExportToWkt(), format="MEM", resampleAlg=resample
    )
    if iFlag_force_resolution_in == 1 and xRes is not None and yRes is not None:
        warp_kwargs["xRes"] = float(xRes)
        warp_kwargs["yRes"] = float(yRes)

    if iFlag_verbose:
        logger.info("Starting gdal.Warp with options: %s", warp_kwargs)

    warp_opts = gdal.WarpOptions(**warp_kwargs)
    warped = gdal.Warp("", src_ds, options=warp_opts)
    if warped is None:
        raise RuntimeError("GDAL warp failed")

    out_x = warped.RasterXSize
    out_y = warped.RasterYSize
    out_bands = warped.RasterCount or 1

    # pick driver from extension (fallback to GTiff)
    try:
        driver = get_raster_driver_from_filename(sFilename_raster_out)
    except Exception:
        driver = gdal.GetDriverByName("GTiff")
        if driver is None:
            raise RuntimeError("Cannot find a suitable GDAL driver for output")

    creation_options = ["COMPRESS=DEFLATE", "PREDICTOR=2"]
    out_ds = driver.Create(
        sFilename_raster_out,
        out_x,
        out_y,
        out_bands,
        eType=iData_type_out,
        options=creation_options,
    )
    if out_ds is None:
        raise RuntimeError("Failed to create output dataset")

    out_ds.SetGeoTransform(warped.GetGeoTransform())
    out_ds.SetProjection(target_srs.ExportToWkt())

    # write bands
    for ib in range(1, out_bands + 1):
        arr = warped.GetRasterBand(ib).ReadAsArray()
        if arr is None:
            raise RuntimeError(f"Failed to read band {ib} from warped dataset")
        # handle nodata mapping
        try:
            np_dtype = gdal_to_numpy_datatype(iData_type_out)
            arr = arr.astype(np_dtype)
        except Exception:
            # fallback: keep array as-is
            pass

        if src_nodata is not None:
            try:
                if np.issubdtype(arr.dtype, np.floating):
                    mask = np.isclose(arr, src_nodata, rtol=1e-9, atol=1e-9)
                else:
                    mask = arr == src_nodata
            except Exception:
                mask = arr == src_nodata
            # write target nodata same as source nodata
            out_nodata = src_nodata
            arr[mask] = out_nodata
            out_ds.GetRasterBand(ib).SetNoDataValue(out_nodata)

        # clamp to source min/max if available
        if src_min is not None:
            arr[arr < src_min] = out_nodata if src_nodata is not None else src_min
        if src_max is not None:
            arr[arr > src_max] = out_nodata if src_nodata is not None else src_max

        out_ds.GetRasterBand(ib).WriteArray(arr)
        out_ds.GetRasterBand(ib).FlushCache()

    # cleanup
    warped = None
    out_ds.FlushCache()
    out_ds = None
    src_ds = None

    return {
        "success": True,
        "output_file": sFilename_raster_out,
        "output_width": out_x,
        "output_height": out_y,
        "output_bands": out_bands,
        "output_projection": target_srs.ExportToWkt(),
        "resampling_algorithm": resample,
        "driver": driver.ShortName if hasattr(driver, "ShortName") else None,
    }
