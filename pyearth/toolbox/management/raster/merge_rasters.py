"""
Raster Merging and Mosaicking Module
====================================

This module provides functionality for merging multiple raster datasets into a single
continuous raster using GDAL's warping capabilities. It supports various resampling
algorithms, handles missing data values, and uses a reference raster or mask to define
the output extent and spatial reference system.

The module is designed for geospatial data processing workflows where multiple
overlapping or adjacent rasters need to be combined into a single dataset with
consistent spatial properties.

Main Functions
--------------
merge_rasters : Merge multiple raster files into a single output raster

Features
--------
- Multiple raster format support (GeoTIFF, JPEG2000, HDF, NetCDF, etc.)
- Flexible resampling algorithm selection (nearest neighbor, bilinear, cubic, etc.)
- Missing value handling with source-to-target value mapping
- Reference mask-based extent and projection definition
- Automatic coordinate system transformation
- Compressed output with DEFLATE compression and predictor optimization
- Support for various GDAL data types (byte, int16, float32, etc.)

Dependencies
------------
- numpy : Array operations and data type conversion
- osgeo.gdal : Raster I/O, warping, and driver management
- osgeo.osr : Spatial reference system handling
- pyearth.gis.gdal : Custom GDAL utility functions

Typical Workflow
----------------
1. Prepare list of input raster files to merge
2. Define or create a reference mask raster with desired extent and projection
3. Specify resolution, resampling algorithm, and data type
4. Call merge_rasters to create merged output

Notes
-----
- All input rasters are reprojected to match the reference mask's spatial reference
- The output extent is determined by the reference mask, not the input rasters
- Missing values are consistently handled across all input rasters
- Output is always compressed GeoTIFF with DEFLATE compression for efficiency

See Also
--------
osgeo.gdal.Warp : GDAL warping functionality
osgeo.gdal.WarpOptions : Warp configuration options
pyearth.gis.gdal.gdal_raster_format_support : Raster format detection and driver utilities
pyearth.gis.gdal.read.raster.gdal_get_raster_extent : Extract raster spatial extent
pyearth.gis.gdal.read.raster.gdal_read_geotiff_file : Read GeoTIFF metadata and data

References
----------
.. [1] GDAL Documentation: https://gdal.org/
.. [2] GDAL Warp: https://gdal.org/programs/gdalwarp.html
.. [3] GeoTIFF Specification: https://www.ogc.org/standards/geotiff

Examples
--------
Basic usage examples are provided in the merge_rasters function docstring.

Author
------
PyEarth Development Team

License
-------
This module is part of the PyEarth package.
"""

import os
import logging
from typing import List, Union, Optional, Literal
import numpy as np
from osgeo import gdal, osr, gdalconst
from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.gdal.gdal_to_numpy_datatype import gdal_to_numpy_datatype
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.gis.gdal.gdal_raster_format_support import (
    get_raster_driver_from_extension,
    get_raster_format_from_extension,
)


def merge_rasters(
    aFilename_rasters: List[str],
    sFilename_merge_out: str,
    dResolution_x: float,
    dResolution_y: float,
    sFilename_mask: str,
    sResampleAlg: Literal[
        "NEAREST",
        "BILINEAR",
        "CUBIC",
        "CUBICSPLINE",
        "LANCZOS",
        "AVERAGE",
        "MODE",
        "MAX",
        "MIN",
        "MED",
        "Q1",
        "Q3",
    ] = "NEAREST",
    dMissing_value_source: Union[int, float] = -9999,
    dMissing_value_target: Union[int, float] = -9999,
    iData_type: int = gdal.GDT_Int16,
    iFlag_overwrite: bool = True,
    iFlag_verbose: bool = False,
) -> dict:
    """Implementation: validate inputs, warp+merge rasters, write output, return stats.

    The implementation intentionally keeps the warped result in-memory and writes
    a single-band output GeoTIFF (or other supported format) with compression.
    """

    logger = logging.getLogger(__name__)
    if iFlag_verbose:
        logging.basicConfig(
            level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
        )

    # ----------------------- Input checks -----------------------
    if not isinstance(aFilename_rasters, list):
        raise TypeError("aFilename_rasters must be a list of file paths")
    if len(aFilename_rasters) == 0:
        raise ValueError("aFilename_rasters cannot be empty")
    for f in aFilename_rasters:
        if not isinstance(f, str):
            raise TypeError("Each element in aFilename_rasters must be a string")
        if not os.path.exists(f):
            raise ValueError(f"Input raster not found: {f}")
        ds = gdal.Open(f, gdal.GA_ReadOnly)
        if ds is None:
            raise ValueError(f"GDAL cannot open input raster: {f}")
        ds = None

    if not isinstance(sFilename_mask, str) or not os.path.exists(sFilename_mask):
        raise ValueError("sFilename_mask must be a valid path to a raster file")

    if not isinstance(sFilename_merge_out, str):
        raise TypeError("sFilename_merge_out must be a string path")

    # Ensure positive resolutions
    if dResolution_x <= 0 or dResolution_y <= 0:
        raise ValueError("Resolutions dResolution_x and dResolution_y must be positive")

    # Validate resampling algorithm (case-insensitive)
    sResampleAlg = sResampleAlg.upper()
    valid_algorithms = [
        "NEAREST",
        "BILINEAR",
        "CUBIC",
        "CUBICSPLINE",
        "LANCZOS",
        "AVERAGE",
        "MODE",
        "MAX",
        "MIN",
        "MED",
        "Q1",
        "Q3",
    ]
    if sResampleAlg not in valid_algorithms:
        raise ValueError(f"Unsupported resampling algorithm: {sResampleAlg}")

    # ----------------------- Prepare output -----------------------
    # Determine output format and driver
    try:
        pDriver_out = get_raster_driver_from_extension(sFilename_merge_out)
    except Exception:
        # fallback to GTiff for safety
        pDriver_out = gdal.GetDriverByName("GTiff")
        if pDriver_out is None:
            raise RuntimeError("Cannot obtain GDAL driver for output format")

    # Remove existing file if overwrite requested
    if os.path.exists(sFilename_merge_out):
        if iFlag_overwrite:
            if iFlag_verbose:
                logger.info(f"Removing existing output file: {sFilename_merge_out}")
            try:
                os.remove(sFilename_merge_out)
            except Exception as e:
                raise OSError(f"Failed to remove existing output: {e}")
        else:
            raise FileExistsError(f"Output file exists: {sFilename_merge_out}")

    # Get extent and grid from mask
    dLon_min, dLon_max, dLat_min, dLat_max = gdal_get_raster_extent(sFilename_mask)
    mask_meta = gdal_read_geotiff_file(sFilename_mask)
    pGeoTransform = mask_meta["geotransform"]
    pProjection_target = mask_meta["projection"]
    ncolumn = mask_meta["ncolumn"]
    nrow = mask_meta["nrow"]

    iNewWidth = int(ncolumn)
    iNewHeight = int(nrow)

    # Create output dataset
    options = ["COMPRESS=DEFLATE", "PREDICTOR=2"]
    pDataset_out = pDriver_out.Create(
        sFilename_merge_out, iNewWidth, iNewHeight, 1, eType=iData_type, options=options
    )
    if pDataset_out is None:
        raise RuntimeError(f"Failed to create output dataset: {sFilename_merge_out}")

    pDataset_out.SetGeoTransform(pGeoTransform)
    pDataset_out.SetProjection(pProjection_target)
    pDataset_out.GetRasterBand(1).SetNoDataValue(dMissing_value_target)

    # ----------------------- Warping / merging -----------------------
    if iFlag_verbose:
        logger.info("Starting gdal.Warp for merging/reprojection")

    # Use target projection WKT as dstSRS
    pWarpOptions = gdal.WarpOptions(
        cropToCutline=False,
        width=iNewWidth,
        height=iNewHeight,
        outputBounds=[dLon_min, dLat_min, dLon_max, dLat_max],
        dstSRS=pProjection_target,
        format="MEM",
        resampleAlg=sResampleAlg,
    )

    pDataset_warped = gdal.Warp("", aFilename_rasters, options=pWarpOptions)
    if pDataset_warped is None:
        pDataset_out = None
        raise RuntimeError("GDAL warp operation failed")

    aData_merged = pDataset_warped.ReadAsArray()
    if aData_merged is None:
        pDataset_warped = None
        pDataset_out = None
        raise RuntimeError("Failed to read warped dataset into array")

    # Handle multi-band -> single band (use first band)
    if aData_merged.ndim == 3:
        # (bands, rows, cols) -> pick first band, warn user
        if iFlag_verbose:
            logger.info(
                f"Warped data contains {aData_merged.shape[0]} bands; using band 1"
            )
        aData_merged = aData_merged[0, :, :]

    # Convert nodata values (handle floats and ints)
    try:
        if np.issubdtype(aData_merged.dtype, np.floating):
            mask_missing = np.isclose(
                aData_merged, dMissing_value_source, rtol=1e-9, atol=1e-9
            )
        else:
            mask_missing = aData_merged == dMissing_value_source
    except Exception:
        # Fallback: direct equality
        mask_missing = aData_merged == dMissing_value_source

    aData_merged[mask_missing] = dMissing_value_target

    # Convert to desired numpy dtype
    iData_type_numpy = gdal_to_numpy_datatype(iData_type)
    aData_merged = aData_merged.astype(iData_type_numpy)

    # Statistics
    num_total = aData_merged.size
    num_nodata = int(np.sum(aData_merged == dMissing_value_target))
    num_valid = int(num_total - num_nodata)
    percent_valid = (num_valid / num_total * 100.0) if num_total > 0 else 0.0

    if iFlag_verbose:
        logger.info(f"Writing merged raster to {sFilename_merge_out}")
        logger.info(f"Valid pixels: {num_valid} ({percent_valid:.2f}%)")

    # Write data to output dataset
    pDataset_out.GetRasterBand(1).WriteArray(aData_merged)
    pDataset_out.GetRasterBand(1).FlushCache()

    # Clean up
    pDataset_warped = None
    pDataset_out = None

    # Map GDAL data type to name
    data_type_names = {
        gdal.GDT_Byte: "Byte",
        gdal.GDT_Int16: "Int16",
        gdal.GDT_UInt16: "UInt16",
        gdal.GDT_Int32: "Int32",
        gdal.GDT_UInt32: "UInt32",
        gdal.GDT_Float32: "Float32",
        gdal.GDT_Float64: "Float64",
    }
    data_type_name = data_type_names.get(iData_type, "Unknown")

    return {
        "success": True,
        "output_file": sFilename_merge_out,
        "num_input_files": len(aFilename_rasters),
        "output_width": iNewWidth,
        "output_height": iNewHeight,
        "output_extent": (dLon_min, dLon_max, dLat_min, dLat_max),
        "output_projection": pProjection_target,
        "resampling_algorithm": sResampleAlg,
        "data_type": data_type_name,
        "nodata_value": dMissing_value_target,
        "num_valid_pixels": num_valid,
        "num_nodata_pixels": num_nodata,
        "percent_valid": float(percent_valid),
    }
