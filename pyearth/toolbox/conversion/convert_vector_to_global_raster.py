"""
Global Raster Conversion for Vector Data

This module provides functionality for converting vector datasets (polygons, polylines)
into global-extent rasters with specified resolution. It is optimized for creating
worldwide raster representations of geographic features while maintaining geographic
coordinate system (latitude/longitude) alignment.

Main Functions
--------------
convert_vector_to_global_raster : Convert vector data to global raster with WGS84 extent

Key Features
------------
- Global extent rasterization (-180 to 180°, -90 to 90°)
- Configurable spatial resolution (latitude/longitude degrees)
- Boundary-only or fill rasterization modes
- Automatic pixel alignment to global grid
- Support for all GDAL-compatible vector formats
- Output directory auto-creation
- Proper handling of edge cases at global boundaries

Use Cases
---------
1. **Global Datasets**: Create worldwide rasters from country boundaries, coastlines
2. **Climate Modeling**: Rasterize climate zones for global climate models
3. **Land Cover**: Convert global land use polygons to raster format
4. **Ocean Features**: Rasterize marine protected areas, shipping routes
5. **Administrative Boundaries**: Create global country/region masks
6. **Data Harmonization**: Convert vector reference data to raster for model input

Technical Details
-----------------
The function ensures that output pixels align with a global grid anchored at
(-180°, 90°). Pixel boundaries are computed using:

    x_min = -180 + floor((lon_min - (-180)) / resolution_x) * resolution_x
    x_max = -180 + ceil((lon_max - (-180)) / resolution_x) * resolution_x
    y_min = 90 - ceil((90 - lat_min) / resolution_y) * resolution_y
    y_max = 90 - floor((90 - lat_max) / resolution_y) * resolution_y

This ensures consistent pixel registration across different datasets at the same
resolution, critical for multi-layer analysis and data fusion.

Resolution Considerations
-------------------------
- High resolution (e.g., 0.001°): ~111m at equator, creates large files
- Medium resolution (e.g., 0.01°): ~1.1km at equator, balance of detail/size
- Low resolution (e.g., 0.1°): ~11km at equator, suitable for global analysis
- Note: Resolution in degrees means actual ground distance varies with latitude

Performance Characteristics
---------------------------
- Time Complexity: O(N*M) where N=features, M=pixels per feature
- Space Complexity: O(rows*columns) for output raster
- Typical file sizes:
  * 0.1° resolution: ~360 x 180 = 64,800 pixels (~63 KB for Byte type)
  * 0.01° resolution: ~3,600 x 1,800 = 6,480,000 pixels (~6.2 MB)
  * 0.001° resolution: ~36,000 x 18,000 = 648,000,000 pixels (~618 MB)

Dependencies
------------
- GDAL/OGR: Vector I/O and rasterization
- NumPy: Grid calculations and array operations
- rasterize_vector: Core rasterization functionality

See Also
--------
- rasterize_vector: General vector to raster conversion with custom extents
- convert_vector_to_raster: Vector to raster with extent from input data
"""

import os
import logging
from typing import Optional
import numpy as np
from osgeo import ogr, osr, gdal

from pyearth.toolbox.conversion.rasterize_vector import rasterize_vector

# Configure logging
logger = logging.getLogger(__name__)


def convert_vector_to_global_raster(
    sFilename_vector_in: str,
    sFilename_tif_out: str,
    dResolution_x_in: float,
    dResolution_y_in: float,
    iFlag_boundary_only_in: int = 0,
    dFill_value_in: float = 2,
) -> None:
    """
    Convert vector dataset to global-extent raster with WGS84 geographic coordinates.

    This function rasterizes vector data (polygons or polylines) to a raster with
    global extent (-180° to 180° longitude, -90° to 90° latitude) at specified
    resolution. The output pixels are automatically aligned to a global grid anchored
    at (-180°, 90°), ensuring consistent pixel registration across datasets.

    Unlike general rasterization functions that use the input data extent, this
    function always creates a global raster covering the entire Earth, making it
    ideal for worldwide datasets and global modeling applications.

    Parameters
    ----------
    sFilename_vector_in : str
        Absolute path to input vector file. Supports all GDAL-compatible formats
        including Shapefile (.shp), GeoJSON (.geojson), GeoPackage (.gpkg), etc.
        File must use geographic coordinates (latitude/longitude). For projected
        data, reproject to EPSG:4326 (WGS84) first.
    sFilename_tif_out : str
        Absolute path for output GeoTIFF file. Must end with '.tif' or '.tiff'.
        Parent directory will be created if it doesn't exist. Existing files are
        overwritten.
    dResolution_x_in : float
        Pixel width in longitudinal degrees. Common values:
        - 0.001 (1/1000°): ~111m at equator, high detail
        - 0.01 (1/100°): ~1.1km at equator, medium detail
        - 0.1 (1/10°): ~11km at equator, coarse detail
        - 1.0 (1°): ~111km at equator, very coarse
        Smaller values create larger files with more detail.
    dResolution_y_in : float
        Pixel height in latitudinal degrees. Typically same as dResolution_x_in
        for square pixels, though rectangular pixels are supported. Positive value
        representing north-south extent.
    iFlag_boundary_only_in : int, optional
        Rasterization mode flag. Default is 0.
        - 0: Fill polygons (interior and boundary pixels set to dFill_value_in)
        - 1: Boundary only (only polygon outline pixels set to dFill_value_in)
        Boundary-only mode useful for visualizing boundaries without masking interior.
    dFill_value_in : float, optional
        Pixel value for rasterized features. Default is 2. All pixels touched by
        vector features are set to this value. Background pixels remain at 0
        (the missing value). Common usage:
        - 1: Binary mask (feature present/absent)
        - 2+: Multiple feature classes can be distinguished
        Must fit within Byte datatype (0-255).

    Returns
    -------
    None
        Output written to `sFilename_tif_out`. Function returns None upon completion.

    Raises
    ------
    FileNotFoundError
        If input vector file does not exist.
    RuntimeError
        If vector file cannot be opened or rasterization fails.
    ValueError
        If resolution values are invalid (<=0 or unreasonably large).
    OSError
        If output directory cannot be created or file cannot be written.

    Notes
    -----
    1. **Global Extent**: The function always creates a raster with extent
       (-180°, -90°) to (180°, 90°), regardless of input data extent. Features
       outside this range are clipped. This ensures global coverage and consistent
       grid alignment across datasets.

    2. **Pixel Alignment**: Output pixels are aligned to a global grid with origin
       at (-180°, 90°). Grid alignment ensures that rasters at the same resolution
       from different vector sources have identical pixel boundaries, critical for
       overlay analysis and data fusion.

    3. **Resolution vs File Size**: File size grows quadratically with resolution:
       - 0.1° resolution: ~64 KB (360 x 180 pixels)
       - 0.01° resolution: ~6.2 MB (3,600 x 1,800 pixels)
       - 0.001° resolution: ~618 MB (36,000 x 18,000 pixels)
       Choose resolution based on analysis needs and available storage.

    4. **Data Type**: Output uses Byte datatype (8-bit unsigned integer, 0-255).
       Suitable for categorical data, masks, and simple classifications. For
       continuous data or more categories, modify the underlying rasterize_vector
       call to use different datatypes.

    5. **Coordinate System**: Function assumes WGS84 geographic coordinates
       (EPSG:4326). Input data should be in lat/lon degrees. Projected data
       must be reprojected before use.

    6. **Boundary Handling**: At global boundaries (-180°/180° and -90°/90°),
       pixel calculations are clamped to valid ranges. Features extending beyond
       global extent are clipped.

    7. **Performance**: Processing time depends on feature complexity and resolution.
       Global 0.01° raster from complex polygons may take several minutes. Use
       lower resolution for initial testing.

    8. **Empty Features**: If input vector file is empty or contains no valid
       geometries, an empty raster (all zeros) is created at global extent.

    9. **Output Format**: Output is always GeoTIFF with LZW compression (if
       supported by rasterize_vector). GeoTIFF includes georeferencing metadata
       for direct use in GIS applications.

    10. **Memory Usage**: The function may require significant memory for high
        resolution rasters. Memory requirement approximately:
        rows * columns * 8 bytes (for processing arrays)

    Examples
    --------
    Create global country mask at 0.1° resolution:

    >>> convert_vector_to_global_raster(
    ...     sFilename_vector_in='/data/countries.shp',
    ...     sFilename_tif_out='/output/global_countries_01deg.tif',
    ...     dResolution_x_in=0.1,
    ...     dResolution_y_in=0.1,
    ...     iFlag_boundary_only_in=0,
    ...     dFill_value_in=1
    ... )
    # Output: 360x180 raster, country polygons filled with value 1

    Create global coastline at 0.01° resolution (boundary only):

    >>> convert_vector_to_global_raster(
    ...     sFilename_vector_in='/data/coastlines.geojson',
    ...     sFilename_tif_out='/output/global_coastlines_001deg.tif',
    ...     dResolution_x_in=0.01,
    ...     dResolution_y_in=0.01,
    ...     iFlag_boundary_only_in=1,
    ...     dFill_value_in=255
    ... )
    # Output: 3600x1800 raster, only coastline boundaries drawn with value 255

    Create high-resolution global marine protected areas:

    >>> convert_vector_to_global_raster(
    ...     sFilename_vector_in='/data/marine_protected_areas.gpkg',
    ...     sFilename_tif_out='/output/global_mpa_0001deg.tif',
    ...     dResolution_x_in=0.001,
    ...     dResolution_y_in=0.001,
    ...     iFlag_boundary_only_in=0,
    ...     dFill_value_in=1
    ... )
    # Output: 36000x18000 raster (~618 MB), MPA polygons filled with value 1
    # WARNING: Large file size, processing may take several minutes

    Expected processing workflow:
    ```
    INFO: Creating output directory: /output
    INFO: Processing vector file: /data/countries.shp
    INFO: Input extent: (-180.0, -90.0) to (180.0, 90.0)
    INFO: Grid extent: (-180.0, -90.0) to (180.0, 90.0)
    INFO: Output dimensions: 360 columns x 180 rows
    INFO: Rasterizing with resolution 0.1° x 0.1°
    INFO: Rasterization completed successfully
    INFO: Output written to: /output/global_countries_01deg.tif
    ```

    See Also
    --------
    rasterize_vector : General vector to raster conversion with custom extents
    convert_vector_to_raster : Vector to raster using input data extent
    """
    # Validate input parameters
    if dResolution_x_in <= 0 or dResolution_y_in <= 0:
        raise ValueError(
            f"Resolution must be positive. Got x={dResolution_x_in}, y={dResolution_y_in}"
        )

    if dResolution_x_in > 10 or dResolution_y_in > 10:
        logger.warning(
            f"Very coarse resolution specified: {dResolution_x_in}° x {dResolution_y_in}°"
        )

    # Validate input file exists
    if not os.path.exists(sFilename_vector_in):
        raise FileNotFoundError(f"Input vector file not found: {sFilename_vector_in}")

    # Create output directory if needed
    sFolder = os.path.dirname(sFilename_tif_out)
    if sFolder and not os.path.exists(sFolder):
        os.makedirs(sFolder)
        logger.info(f"Created output directory: {sFolder}")

    # Remove existing output file if present
    if os.path.exists(sFilename_tif_out):
        # Check if file is empty (incomplete previous run)
        if os.path.getsize(sFilename_tif_out) == 0:
            logger.info(f"Removing empty file: {sFilename_tif_out}")
        else:
            logger.info(f"Removing existing file: {sFilename_tif_out}")
        os.remove(sFilename_tif_out)

    # Open input vector dataset
    pDataSource_clip = ogr.Open(sFilename_vector_in)
    if pDataSource_clip is None:
        raise RuntimeError(f"Could not open vector file: {sFilename_vector_in}")

    # Get layer and verify it contains features
    pLayer_clip = pDataSource_clip.GetLayer()
    if pLayer_clip is None:
        raise RuntimeError("Could not access layer in vector file")

    nFeature_count = pLayer_clip.GetFeatureCount()
    logger.info(
        f"Processing vector file with {nFeature_count} features: {sFilename_vector_in}"
    )

    # Get data extent (informational only - we'll use global extent)
    dLon_min, dLon_max, dLat_min, dLat_max = pLayer_clip.GetExtent()
    logger.info(
        f"Input extent: ({dLon_min:.4f}, {dLat_min:.4f}) to ({dLon_max:.4f}, {dLat_max:.4f})"
    )

    # Force global extent (WGS84 geographic bounds)
    # This ensures consistent grid alignment for global datasets
    dLon_min = -180.0
    dLon_max = 180.0
    dLat_min = -90.0
    dLat_max = 90.0

    logger.info(f"Grid extent: ({dLon_min}, {dLat_min}) to ({dLon_max}, {dLat_max})")

    # Calculate pixel-aligned grid boundaries
    # Align to global grid with origin at (-180, 90) to ensure consistent registration
    nleft = np.floor((dLon_min - (-180)) / dResolution_x_in)
    nright = np.ceil((dLon_max - (-180)) / dResolution_x_in)
    ntop = np.floor((90 - dLat_max) / dResolution_y_in)
    nbot = np.ceil((90 - dLat_min) / dResolution_y_in)

    # Calculate aligned grid coordinates
    dMin_x = -180 + nleft * dResolution_x_in
    dMax_x = -180 + nright * dResolution_x_in
    dMin_y = 90 - nbot * dResolution_y_in
    dMax_y = 90 - ntop * dResolution_y_in

    # Clamp to global extent boundaries (handle floating point edge cases)
    if dMin_x < -180:
        dMin_x = -180.0
    if dMax_x > 180:
        dMax_x = 180.0
    if dMin_y < -90:
        dMin_y = -90.0
    if dMax_y > 90:
        dMax_y = 90.0

    # Calculate output raster dimensions
    nrow = int(nbot - ntop)
    ncolumn = int(nright - nleft)

    logger.info(f"Output dimensions: {ncolumn} columns x {nrow} rows")
    logger.info(
        f"Rasterizing with resolution {dResolution_x_in}° x {dResolution_y_in}°"
    )

    # Estimate output file size
    estimated_size_mb = (nrow * ncolumn) / (1024 * 1024)
    if estimated_size_mb > 100:
        logger.warning(f"Large output file expected: ~{estimated_size_mb:.1f} MB")

    # Perform rasterization using the general rasterize_vector function
    # with global extent and calculated grid parameters
    try:
        rasterize_vector(
            sFilename_vector_in,
            sFilename_tif_out,
            dResolution_x_in,
            dResolution_y_in,
            dMissing_value_in=0,
            iDataType_out=gdal.GDT_Byte,
            dMin_x_in=dMin_x,
            dMax_x_in=dMax_x,
            dMin_y_in=dMin_y,
            dMax_y_in=dMax_y,
            nRow_in=nrow,
            nColumn_in=ncolumn,
            iFlag_boundary_only_in=iFlag_boundary_only_in,
            dFill_value_in=dFill_value_in,
        )
        logger.info("Rasterization completed successfully")
        logger.info(f"Output written to: {sFilename_tif_out}")
    except Exception as e:
        logger.error(f"Rasterization failed: {e}")
        raise RuntimeError(f"Failed to rasterize vector to global raster: {e}")
    finally:
        # Clean up data source handle
        pDataSource_clip = None

    return None
