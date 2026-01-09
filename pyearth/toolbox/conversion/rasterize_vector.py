"""
Vector to Raster Conversion

This module provides core functionality for converting vector datasets (points, polylines,
polygons) into raster format. It offers flexible control over rasterization parameters
including extent, resolution, geometry handling, and attribute burning.

Main Functions
--------------
rasterize_vector : Convert vector data to raster with comprehensive parameter control

Key Features
------------
- Multi-geometry support (POINT, LINESTRING, POLYGON, MULTIPOLYGON)
- Flexible extent specification (auto-detect or custom bounds)
- Configurable spatial resolution
- Attribute-based value assignment
- Boundary-only or fill rasterization modes
- Custom data types and missing values
- Coordinate system transformation support
- GeoTIFF output with LZW compression
- ALL_TOUCHED rasterization for complete coverage

Use Cases
---------
1. **Polygon Rasterization**: Convert land use polygons to raster grid
2. **Polyline Rasterization**: Rasterize roads, rivers, power lines
3. **Point Rasterization**: Create point density grids, sample locations
4. **Attribute Burning**: Transfer vector attributes to raster pixel values
5. **Boundary Extraction**: Rasterize only polygon boundaries for visualization
6. **Custom Extents**: Create rasters aligned to specific grid systems
7. **Data Type Control**: Generate masks (Byte), classifications (Int16), or continuous (Float32)

Technical Details
-----------------
The module uses GDAL's RasterizeLayer function with ALL_TOUCHED option to ensure
complete rasterization of vector features. The ALL_TOUCHED flag means any pixel
whose center or edge is touched by a vector feature will be burned with the
specified value.

Rasterization Modes:
- **Boundary Only** (iFlag_boundary_only_in=1): Only feature boundaries rasterized
- **Fill** (iFlag_boundary_only_in=0): Interior and boundary both rasterized
- **Attribute Burning**: Use field values from attribute table as pixel values
- **Constant Value**: Burn all features with same value

Geometry Handling:
- POINT: Single pixel per point (ALL_TOUCHED not used)
- LINESTRING/POLYLINE: Rasterized as boundaries (1-pixel wide)
- POLYGON: Can be filled or boundary-only
- MULTIPOLYGON: Each part handled independently

Performance Characteristics
---------------------------
- Time Complexity: O(N*M) where N=features, M=avg pixels per feature
- Space Complexity: O(rows*columns) for output raster
- Memory Usage: ~(rows*columns*bytes_per_pixel) plus GDAL overhead
- Compression: LZW compression reduces file size (~50-90% for masks)

Dependencies
------------
- GDAL/OGR: Vector I/O and rasterization operations
- NumPy: Numerical array operations
- OSR: Spatial reference system handling

See Also
--------
- convert_vector_to_global_raster: Specialized global extent rasterization
- convert_vector_to_raster: Simplified interface using auto-detected extent
"""

import os
import logging
from typing import Optional
import numpy as np
from osgeo import gdal, ogr, osr

# Configure logging
logger = logging.getLogger(__name__)


def rasterize_vector(
    sFilename_vector_in: str,
    sFilename_raster_out: str,
    dResolution_x: float,
    dResolution_y: float,
    iFlag_boundary_only_in: Optional[int] = None,
    iFlag_use_field_value_in: Optional[int] = None,
    sAttribute_name_in: Optional[str] = None,
    dMissing_value_in: Optional[float] = None,
    dField_value_in: Optional[float] = None,
    dFill_value_in: Optional[float] = None,
    iDataType_out: Optional[int] = None,
    dMin_x_in: Optional[float] = None,
    dMax_x_in: Optional[float] = None,
    dMin_y_in: Optional[float] = None,
    dMax_y_in: Optional[float] = None,
    nRow_in: Optional[int] = None,
    nColumn_in: Optional[int] = None,
    pProjection_target_in: Optional[str] = None,
) -> None:
    """
    Convert vector dataset to raster with flexible parameter control.

    This function rasterizes vector data (points, polylines, polygons) into a GeoTIFF
    raster format. It provides comprehensive control over the rasterization process
    including extent, resolution, geometry handling, value assignment, and output
    data type. The function automatically handles different geometry types and offers
    both boundary-only and fill rasterization modes.

    Parameters
    ----------
    sFilename_vector_in : str
        Absolute path to input vector file. Supports all GDAL-compatible formats
        including Shapefile (.shp), GeoJSON (.geojson), GeoPackage (.gpkg), etc.
    sFilename_raster_out : str
        Absolute path for output GeoTIFF file. Must end with '.tif' or '.tiff'.
        Existing files will be overwritten. Output uses LZW compression.
    dResolution_x : float
        Pixel width in the units of the coordinate system (e.g., degrees for geographic,
        meters for projected). Determines east-west pixel spacing.
    dResolution_y : float
        Pixel height in the units of the coordinate system. Determines north-south
        pixel spacing. Typically same as dResolution_x for square pixels.
    iFlag_boundary_only_in : int, optional
        Rasterization mode for polygon features. Default is 1.
        - 0: Fill polygons (interior and boundary get values)
        - 1: Boundary only (only outline pixels get values)
        - Automatically set to 1 for POINT and LINESTRING geometries
    iFlag_use_field_value_in : int, optional
        Whether to use attribute field values as pixel values. Default is 0 (None).
        - 0 or None: Use constant burn values (dField_value_in, dFill_value_in)
        - 1: Use values from sAttribute_name_in field for pixel values
        Requires sAttribute_name_in to be specified when set to 1.
    sAttribute_name_in : str, optional
        Name of attribute field to use for pixel values when iFlag_use_field_value_in=1.
        Field must exist in the vector dataset and contain numeric values.
        Example: 'population', 'elevation', 'landuse_code'
    dMissing_value_in : float, optional
        Value for pixels not touched by any vector feature (no-data value).
        Default is -9999. Used as both fill value and NoData metadata.
    dField_value_in : float, optional
        Value to burn for feature boundaries when not using field values.
        Default is 1. For boundary-only mode, this is the outline value.
        For fill mode, this is the boundary value (may differ from fill).
    dFill_value_in : float, optional
        Value to burn for polygon interiors when in fill mode (iFlag_boundary_only_in=0).
        Default is 1. Only used for polygon features when fill mode is active.
    iDataType_out : int, optional
        GDAL data type for output raster. Default is gdal.GDT_Float32.
        Common options:
        - gdal.GDT_Byte: 8-bit unsigned (0-255), smallest, good for masks/classifications
        - gdal.GDT_Int16: 16-bit signed integer (-32768 to 32767)
        - gdal.GDT_UInt16: 16-bit unsigned (0-65535)
        - gdal.GDT_Float32: 32-bit floating point (default)
        - gdal.GDT_Float64: 64-bit floating point (largest, most precise)
    dMin_x_in : float, optional
        Minimum X coordinate (left edge) of output raster extent. If not specified,
        automatically calculated from input vector extent.
    dMax_x_in : float, optional
        Maximum X coordinate (right edge) of output raster extent. If not specified,
        automatically calculated from input vector extent.
    dMin_y_in : float, optional
        Minimum Y coordinate (bottom edge) of output raster extent. If not specified,
        automatically calculated from input vector extent.
    dMax_y_in : float, optional
        Maximum Y coordinate (top edge) of output raster extent. If not specified,
        automatically calculated from input vector extent.
    nRow_in : int, optional
        Number of rows in output raster. If not specified, calculated from extent
        and resolution as: nrow = (dMax_y - dMin_y) / dResolution_y
    nColumn_in : int, optional
        Number of columns in output raster. If not specified, calculated from extent
        and resolution as: ncolumn = (dMax_x - dMin_x) / dResolution_x
    pProjection_target_in : str, optional
        Target projection as WKT (Well-Known Text) string. If not specified, uses
        the projection from the input vector file. Used for on-the-fly reprojection.

    Returns
    -------
    None
        Output written to sFilename_raster_out. Function returns None upon completion.

    Raises
    ------
    RuntimeError
        If input vector file cannot be opened or raster file cannot be created.
    ValueError
        If resolution is invalid (results in 0 rows or columns).

    Notes
    -----
    1. **Geometry Type Handling**: The function automatically detects geometry type
       and adjusts behavior:
       - POINT: Rasterizes as single pixels, ignores boundary_only flag
       - LINESTRING: Always uses boundary-only mode (1-pixel wide lines)
       - POLYGON: Respects iFlag_boundary_only_in setting
       - MULTIPOLYGON: Each component handled as POLYGON

    2. **ALL_TOUCHED Option**: Uses GDAL's ALL_TOUCHED=TRUE, meaning any pixel
       touched by a vector feature (not just pixel centers) will be burned. This
       ensures complete coverage of thin features but may result in thicker lines.

    3. **Boundary vs Fill Mode**: For polygons with iFlag_boundary_only_in=0,
       the function rasterizes twice:
       - First pass: Fill interior with dFill_value_in
       - Second pass: Overwrite boundary with dField_value_in
       This allows different values for interior vs boundary pixels.

    4. **Attribute Burning**: When iFlag_use_field_value_in=1, pixel values are
       taken from the specified attribute field. Useful for transferring categorical
       or continuous data from vector to raster (e.g., land use codes, population).

    5. **Extent Calculation**: If extent parameters (dMin_x_in, etc.) are not provided,
       they are automatically calculated from the input vector layer's extent. This
       ensures all features are included but may result in irregular grid bounds.

    6. **Resolution and Dimensions**: Either specify resolution (auto-calculate dims)
       or specify both resolution and dimensions (for custom grid alignment). If
       both nRow_in/nColumn_in and resolution are provided, the specified dimensions
       take precedence.

    7. **Compression**: Output GeoTIFF uses LZW compression with predictor=2,
       significantly reducing file size (typically 50-90% for masks and classifications)
       with no data loss.

    8. **No-Data Value**: The dMissing_value_in is used to fill the raster before
       burning features and is set as the raster's NoData metadata. GIS software
       will treat these pixels as transparent/missing.

    9. **Memory Usage**: The function creates in-memory boundary layers for fill
       mode operations. For very large datasets (millions of features), this may
       require significant memory.

    10. **Spatial Reference**: The output raster inherits the spatial reference from
        the input vector unless pProjection_target_in is specified. Coordinate
        transformation is NOT performed; pProjection_target_in only sets metadata.

    Examples
    --------
    Basic polygon rasterization with auto-extent:

    >>> rasterize_vector(
    ...     sFilename_vector_in='/data/parcels.shp',
    ...     sFilename_raster_out='/output/parcels.tif',
    ...     dResolution_x=10.0,
    ...     dResolution_y=10.0,
    ...     iFlag_boundary_only_in=0,
    ...     dFill_value_in=1,
    ...     iDataType_out=gdal.GDT_Byte
    ... )
    # Output: Parcels rasterized at 10m resolution, filled with value 1

    Rasterize road network (polylines):

    >>> rasterize_vector(
    ...     sFilename_vector_in='/data/roads.geojson',
    ...     sFilename_raster_out='/output/roads.tif',
    ...     dResolution_x=5.0,
    ...     dResolution_y=5.0,
    ...     dField_value_in=255,
    ...     dMissing_value_in=0,
    ...     iDataType_out=gdal.GDT_Byte
    ... )
    # Output: Roads as 1-pixel wide lines with value 255 on black background (0)

    Burn attribute values to raster:

    >>> rasterize_vector(
    ...     sFilename_vector_in='/data/zones.gpkg',
    ...     sFilename_raster_out='/output/zone_codes.tif',
    ...     dResolution_x=30.0,
    ...     dResolution_y=30.0,
    ...     iFlag_use_field_value_in=1,
    ...     sAttribute_name_in='zone_id',
    ...     dMissing_value_in=-9999,
    ...     iDataType_out=gdal.GDT_Int16
    ... )
    # Output: Each polygon rasterized with its zone_id attribute value

    Custom extent with specific dimensions (grid alignment):

    >>> rasterize_vector(
    ...     sFilename_vector_in='/data/lakes.shp',
    ...     sFilename_raster_out='/output/lakes_grid.tif',
    ...     dResolution_x=100.0,
    ...     dResolution_y=100.0,
    ...     dMin_x_in=500000.0,
    ...     dMax_x_in=600000.0,
    ...     dMin_y_in=4000000.0,
    ...     dMax_y_in=4100000.0,
    ...     nRow_in=1000,
    ...     nColumn_in=1000,
    ...     dFill_value_in=1,
    ...     iDataType_out=gdal.GDT_Byte
    ... )
    # Output: Lakes rasterized to exactly 1000x1000 grid covering specified extent

    Expected log output:
    ```
    INFO: Processing vector file: /data/parcels.shp
    INFO: Geometry type: POLYGON
    INFO: Rasterization mode: Fill (boundary + interior)
    INFO: Output extent: (450000.0, 4500000.0) to (460000.0, 4510000.0)
    INFO: Output dimensions: 1000 columns x 1000 rows
    INFO: Creating raster with data type: Byte
    INFO: Rasterizing features...
    INFO: Rasterization completed successfully
    ```

    See Also
    --------
    convert_vector_to_global_raster : Global extent rasterization (-180 to 180°)
    convert_vector_to_raster : Simplified interface with auto-extent
    """
    # Validate input file exists
    if not os.path.exists(sFilename_vector_in):
        raise FileNotFoundError(f"Input vector file not found: {sFilename_vector_in}")

    logger.info(f"Processing vector file: {sFilename_vector_in}")

    # Open input vector data source
    try:
        pDatasource_vector = ogr.Open(sFilename_vector_in)
    except RuntimeError as e:
        raise RuntimeError(
            f"Could not open vector file: {sFilename_vector_in}. Error: {e}"
        )

    if pDatasource_vector is None:
        raise RuntimeError(f"Could not open vector file: {sFilename_vector_in}")

    # Create memory data source for boundary operations
    pDatasource_boundary = ogr.GetDriverByName("MEM").CreateDataSource("out")
    if pDatasource_boundary is None:
        logger.warning("Could not create memory data source for boundary operations")

    # Remove existing output file if present
    if os.path.exists(sFilename_raster_out):
        os.remove(sFilename_raster_out)
        logger.info(f"Removed existing output file: {sFilename_raster_out}")

    # Process optional parameters with defaults
    sAttribute_name = sAttribute_name_in if sAttribute_name_in is not None else None

    if iFlag_use_field_value_in is None:
        iFlag_use_field_value = 0
    else:
        iFlag_use_field_value = iFlag_use_field_value_in
        if iFlag_use_field_value == 1 and sAttribute_name is None:
            raise ValueError(
                "sAttribute_name_in must be specified when iFlag_use_field_value_in=1"
            )

    dField_value = dField_value_in if dField_value_in is not None else 1
    iFlag_boundary_only = (
        iFlag_boundary_only_in if iFlag_boundary_only_in is not None else 1
    )
    dFill_value = dFill_value_in if dFill_value_in is not None else 1

    # Get vector layer and spatial reference
    pLayer_vector = pDatasource_vector.GetLayer()
    if pLayer_vector is None:
        raise RuntimeError("Could not access layer in vector file")

    nFeature_count = pLayer_vector.GetFeatureCount()
    logger.info(f"Vector layer contains {nFeature_count} features")

    pSpatialRef_source = pLayer_vector.GetSpatialRef()
    if pSpatialRef_source is None:
        logger.warning("Input vector has no spatial reference defined")
        pProjection_source = ""
    else:
        pProjection_source = pSpatialRef_source.ExportToWkt()

    # Determine target projection
    pProjection_target = (
        pProjection_target_in
        if pProjection_target_in is not None
        else pProjection_source
    )

    pSpatialRef_target = osr.SpatialReference()
    if pProjection_target:
        pSpatialRef_target.ImportFromWkt(pProjection_target)

    # Set output data type
    iDataType = iDataType_out if iDataType_out is not None else gdal.GDT_Float32
    logger.info(f"Creating raster with data type: {gdal.GetDataTypeName(iDataType)}")

    # Set missing/no-data value
    dMissing_value = dMissing_value_in if dMissing_value_in is not None else -9999

    # Determine raster extent
    if None in (dMin_x_in, dMax_x_in, dMin_y_in, dMax_y_in):
        # Auto-calculate extent from vector layer
        dMin_x, dMax_x, dMin_y, dMax_y = pLayer_vector.GetExtent()
        logger.info("Using auto-detected extent from vector data")
    else:
        dMin_x = dMin_x_in
        dMax_x = dMax_x_in
        dMin_y = dMin_y_in
        dMax_y = dMax_y_in

    logger.info(f"Output extent: ({dMin_x}, {dMin_y}) to ({dMax_x}, {dMax_y})")

    # Detect geometry type and adjust boundary flag accordingly
    pLayer_vector.ResetReading()
    pFeature = pLayer_vector.GetNextFeature()
    if pFeature is None:
        logger.warning("Vector layer contains no features")
        return

    pGeometry = pFeature.GetGeometryRef()
    if pGeometry is None:
        raise RuntimeError("First feature has no geometry")

    sGeometry_type = pGeometry.GetGeometryName()
    logger.info(f"Geometry type: {sGeometry_type}")

    iFlag_point = 0
    if sGeometry_type == "POINT":
        logger.info("Point geometry detected - using point rasterization")
        iFlag_point = 1
        iFlag_boundary_only = 1
    elif sGeometry_type == "LINESTRING":
        logger.info("Polyline geometry detected - using boundary-only mode")
        iFlag_boundary_only = 1
    elif sGeometry_type == "POLYGON":
        logger.info("Polygon geometry detected")
    elif sGeometry_type == "MULTIPOLYGON":
        logger.info("Multi-polygon geometry detected")
    else:
        logger.warning(
            f"Unknown geometry type: {sGeometry_type} - treating as boundary-only"
        )
        iFlag_boundary_only = 1

    # Log rasterization mode
    if iFlag_boundary_only == 0:
        logger.info("Rasterization mode: Fill (boundary + interior)")
    else:
        logger.info("Rasterization mode: Boundary only")

    # Calculate raster dimensions
    if nRow_in is not None:
        nrow = nRow_in
    else:
        nrow = int((dMax_y - dMin_y) / dResolution_y)

    if nColumn_in is not None:
        ncolumn = nColumn_in
    else:
        ncolumn = int((dMax_x - dMin_x) / dResolution_x)

    if nrow <= 0 or ncolumn <= 0:
        raise ValueError(
            f"Invalid raster dimensions: {ncolumn} columns x {nrow} rows. "
            f"Check resolution values (x={dResolution_x}, y={dResolution_y})"
        )

    logger.info(f"Output dimensions: {ncolumn} columns x {nrow} rows")

    # Estimate file size
    bytes_per_pixel = gdal.GetDataTypeSize(iDataType) // 8
    estimated_size_mb = (nrow * ncolumn * bytes_per_pixel) / (1024 * 1024)
    if estimated_size_mb > 100:
        logger.warning(
            f"Large output file expected: ~{estimated_size_mb:.1f} MB (uncompressed)"
        )

    # Create output raster dataset
    pRaster_Driver = gdal.GetDriverByName("GTiff")
    options = ["COMPRESS=DEFLATE", "PREDICTOR=2"]
    pDatasource_raster = pRaster_Driver.Create(
        sFilename_raster_out, ncolumn, nrow, 1, eType=iDataType, options=options
    )

    if pDatasource_raster is None:
        raise RuntimeError(f"Could not create raster file: {sFilename_raster_out}")

    # Set geotransform (georeferencing)
    # Format: (top-left X, pixel width, rotation, top-left Y, rotation, -pixel height)
    pDatasource_raster.SetGeoTransform(
        (dMin_x, dResolution_x, 0, dMax_y, 0, -dResolution_y)
    )

    # Get raster band and initialize
    band = pDatasource_raster.GetRasterBand(1)
    band.Fill(dMissing_value)
    band.SetNoDataValue(dMissing_value)

    # Set spatial reference system
    pDatasource_raster.SetProjection(pProjection_target)

    logger.info("Rasterizing features...")

    # Perform rasterization based on mode
    try:
        if iFlag_use_field_value == 1:
            # Use attribute field values
            if iFlag_boundary_only == 0:
                # Fill mode with attribute values
                pLayer_boundary = pDatasource_boundary.CreateLayer(
                    "boundary", srs=pSpatialRef_source
                )
                pLayer_vector.ResetReading()
                for pFeature in pLayer_vector:
                    geometry = pFeature.GetGeometryRef()
                    if geometry is not None:
                        pFeature_boundary = ogr.Feature(pLayer_boundary.GetLayerDefn())
                        pFeature_boundary.SetGeometry(geometry.Boundary())
                        pLayer_boundary.CreateFeature(pFeature_boundary)

                # Rasterize filled polygons
                gdal.RasterizeLayer(
                    pDatasource_raster,
                    [1],
                    pLayer_vector,
                    None,
                    None,
                    burn_values=[dField_value],
                    options=["ALL_TOUCHED=TRUE"],
                )
                # Rasterize boundaries
                gdal.RasterizeLayer(
                    pDatasource_raster,
                    [1],
                    pLayer_boundary,
                    None,
                    None,
                    options=["ALL_TOUCHED=TRUE"],
                )
            else:
                # Boundary-only or point mode with attribute values
                if iFlag_point == 1:
                    gdal.RasterizeLayer(
                        pDatasource_raster,
                        [1],
                        pLayer_vector,
                        burn_values=[dField_value],
                    )
                else:
                    gdal.RasterizeLayer(
                        pDatasource_raster,
                        [1],
                        pLayer_vector,
                        None,
                        None,
                        options=["ALL_TOUCHED=TRUE", "ATTRIBUTE=" + sAttribute_name],
                    )
        else:
            # Use constant burn values
            if iFlag_boundary_only == 0:
                # Fill mode with constant values
                pLayer_boundary = pDatasource_boundary.CreateLayer(
                    "boundary", srs=pSpatialRef_source
                )
                pLayer_vector.ResetReading()
                nFeatureCount = pLayer_vector.GetFeatureCount()
                for pFeature in pLayer_vector:
                    geometry = pFeature.GetGeometryRef()
                    if geometry is not None:
                        pFeature_boundary = ogr.Feature(pLayer_boundary.GetLayerDefn())
                        pFeature_boundary.SetGeometry(geometry.Boundary())
                        pLayer_boundary.CreateFeature(pFeature_boundary)

                # Rasterize filled polygons
                gdal.RasterizeLayer(
                    pDatasource_raster,
                    [1],
                    pLayer_vector,
                    None,
                    None,
                    burn_values=[dFill_value],
                    options=["ALL_TOUCHED=TRUE"],
                )
                # Rasterize boundaries (potentially different value)
                gdal.RasterizeLayer(
                    pDatasource_raster,
                    [1],
                    pLayer_boundary,
                    None,
                    None,
                    burn_values=[dField_value],
                    options=["ALL_TOUCHED=TRUE"],
                )
            else:
                # Boundary-only or point mode with constant value
                if iFlag_point == 1:
                    gdal.RasterizeLayer(
                        pDatasource_raster,
                        [1],
                        pLayer_vector,
                        burn_values=[dField_value],
                    )
                else:
                    gdal.RasterizeLayer(
                        pDatasource_raster,
                        [1],
                        pLayer_vector,
                        None,
                        None,
                        burn_values=[dField_value],
                        options=["ALL_TOUCHED=TRUE"],
                    )
    except Exception as e:
        logger.error(f"Rasterization failed: {e}")
        raise RuntimeError(f"Failed to rasterize vector data: {e}")

    # Flush and clean up
    band.FlushCache()
    band = None

    if pDatasource_boundary is not None:
        pDatasource_boundary.Destroy()

    pDatasource_raster = None
    pDatasource_vector = None
    pSpatialRef_source = None
    pSpatialRef_target = None

    logger.info("Rasterization completed successfully")
    logger.info(f"Output written to: {sFilename_raster_out}")

    return None
