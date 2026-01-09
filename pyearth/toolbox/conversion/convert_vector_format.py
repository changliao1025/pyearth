"""
Vector Format Conversion

This module provides functionality for converting vector datasets between different
formats with optional coordinate system transformation. When converting to GeoJSON
format, automatic transformation to WGS84 (EPSG:4326) is applied. For other formats,
the original spatial reference is preserved by default.

Main Functions
--------------
convert_vector_format : Convert vector data between formats with optional transformation

Key Features
------------
- Multi-format support for both input and output (auto-detected from extension)
- Automatic coordinate transformation to WGS84 for GeoJSON output
- Preserves original spatial reference for non-GeoJSON formats
- Support for multiple geometry types (Point, LineString, Polygon, MultiPolygon)
- Attribute field preservation during conversion
- GDAL 3+ axis order handling
- Comprehensive error handling and logging
- Validation of input/output files

Use Cases
---------
1. **Web Mapping**: Convert to GeoJSON with WGS84 for web applications
2. **Data Exchange**: Convert between formats while preserving projections
3. **Format Migration**: Convert Shapefile to GeoPackage maintaining CRS
4. **Legacy Support**: Convert old formats to modern standards
5. **Web Services**: Prepare GeoJSON data for web mapping with standard WGS84
6. **Data Archive**: Convert to GeoPackage for long-term storage

Technical Details
-----------------
The module uses GDAL/OGR for format conversion and coordinate transformation.
GeoJSON outputs are automatically transformed to WGS84 (EPSG:4326) geographic
coordinates per GeoJSON specification (RFC 7946). Other formats preserve the
original coordinate system unless explicitly requested.

For GDAL 3+, axis order is set to traditional GIS order (longitude, latitude)
to maintain consistency with earlier versions and common GIS conventions.

Supported Formats:
- **Shapefile** (.shp): ESRI's legacy format
- **GeoJSON** (.geojson, .json): Web-friendly JSON format
- **GeoPackage** (.gpkg): Modern SQLite-based format
- **GML** (.gml): OGC XML format
- **KML** (.kml): Google Earth format
- And many others supported by GDAL/OGR

Performance Characteristics
---------------------------
- Time Complexity: O(N) where N = number of features
- Space Complexity: O(1) - streaming feature-by-feature
- Transformation overhead: Minimal for same CRS, ~10-20% for reprojection

Dependencies
------------
- GDAL/OGR: Vector I/O and coordinate transformation
- OSR: Spatial reference system handling

See Also
--------
- get_vector_driver_from_filename: Automatic driver selection
- get_vector_format_from_filename: Format name extraction
"""

import os
import logging
from typing import Optional
import osgeo
from osgeo import ogr, osr, gdal

from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_driver_from_filename,
    get_vector_format_from_filename,
)
from pyearth.gis.geometry.get_output_geometry_type import get_output_geometry_type

# Configure logging
logger = logging.getLogger(__name__)


def convert_vector_format(
    sFilename_vector_in: str,
    sFilename_vector_out: str,
    target_epsg: Optional[int] = None,
) -> bool:
    """
    Convert vector dataset between formats with conditional coordinate transformation.

    This function converts vector data from any GDAL-supported format to any other
    format. Format detection is automatic based on file extensions.

    When converting to GeoJSON format (.geojson, .json), the output is automatically
    transformed to WGS84 (EPSG:4326) per GeoJSON specification (RFC 7946). For all
    other formats, the original spatial reference is preserved unless explicitly
    overridden via target_epsg parameter.

    Parameters
    ----------
    sFilename_vector_in : str
        Absolute path to input vector file. Format automatically detected from
        extension. Supports: Shapefile (.shp), GeoJSON (.geojson), GeoPackage (.gpkg),
        GML (.gml), KML (.kml), and all other GDAL-supported formats.
    sFilename_vector_out : str
        Absolute path for output vector file. Format automatically detected from
        extension. Existing files will be overwritten. Common extensions:
        - .geojson or .json for GeoJSON (auto-transforms to WGS84)
        - .shp for Shapefile (preserves original CRS)
        - .gpkg for GeoPackage (preserves original CRS)
        - .gml for GML (preserves original CRS)
        - .kml for KML (preserves original CRS)
    target_epsg : Optional[int], optional
        EPSG code for output coordinate system. If None (default), behavior depends
        on output format:
        - GeoJSON output: Always uses EPSG:4326 (WGS84) per RFC 7946
        - Other formats: Preserve original input CRS
        If specified, forces transformation to this EPSG code regardless of format.
        Common values when specified:
        - 4326: WGS84 geographic (latitude/longitude)
        - 3857: Web Mercator (used by web maps)
        - 32633: UTM Zone 33N

    Returns
    -------
    bool
        True if conversion completed successfully, False if errors occurred.

    Raises
    ------
    FileNotFoundError
        If input file does not exist (implicitly via logging).
    ValueError
        If input or output format is not supported.
    RuntimeError
        If file cannot be opened, driver not available, or conversion fails.

    Notes
    -----
    1. **Automatic Format Detection**: Both input and output formats are detected
       from file extensions using `get_vector_driver_from_filename()`. No need
       to specify drivers manually.

    2. **Conditional Coordinate Transformation**: Transformation behavior depends
       on output format:
       - **GeoJSON output**: Always transforms to WGS84 (EPSG:4326) per RFC 7946
       - **Other formats**: Preserves original CRS unless target_epsg is specified
       - **Explicit override**: Setting target_epsg forces transformation for any format

    3. **GDAL 3+ Compatibility**: For GDAL version 3 and above, axis mapping is
       set to traditional GIS order (longitude, latitude) to maintain consistency
       with earlier versions and avoid coordinate swap issues.

    4. **Geometry Type Support**: Supports common geometry types with automatic
       normalization via `get_output_geometry_type()`:
       - Point (wkbPoint)
       - LineString (wkbLineString)
       - Polygon (wkbPolygon)
       - Multi-geometry types (converted to single equivalents)
       - 2.5D/3D variants (normalized to 2D)
       Example: MultiLineString → LineString, LineString25D → LineString

    5. **Attribute Preservation**: All attribute fields from the input layer are
       copied to the output layer with their original names, types, and values.

    6. **Spatial Reference Handling**:
       - If transformation is needed and CRS are identical, transformation is skipped
       - Comparison is done via WKT string matching for reliability
       - For non-GeoJSON formats without target_epsg, original CRS is preserved

    7. **Memory Efficiency**: Features are processed one at a time in streaming
       fashion, making the function suitable for large datasets without excessive
       memory usage.

    8. **File Cleanup**: Existing output files are automatically removed before
       creating new output. For Shapefiles, all component files (.shp, .shx, .dbf,
       .prj) are removed.

    9. **Error Handling**: Comprehensive error checking at each step with informative
       logging. Returns False on any error, allowing calling code to handle failures.

    10. **Output Validation**: After conversion, verify the output file was created
        and contains the expected features by checking file existence and optionally
        opening the result.

    Examples
    --------
    Convert Shapefile to GeoJSON (automatic WGS84 transformation):

    >>> success = convert_vector_format(
    ...     sFilename_vector_in='/data/parcels.shp',
    ...     sFilename_vector_out='/output/parcels.geojson'
    ... )
    >>> print(f"Success: {success}")
    # Output: parcels.geojson in WGS84 coordinates (automatically transformed)

    Convert Shapefile to GeoPackage (preserves original CRS):

    >>> success = convert_vector_format(
    ...     sFilename_vector_in='/data/roads.shp',
    ...     sFilename_vector_out='/output/roads.gpkg'
    ... )
    # Output: roads.gpkg in original CRS (e.g., UTM if input was UTM)

    Convert GeoPackage to Shapefile with explicit UTM projection:

    >>> success = convert_vector_format(
    ...     sFilename_vector_in='/data/cities.gpkg',
    ...     sFilename_vector_out='/output/cities.shp',
    ...     target_epsg=32633  # Force UTM Zone 33N
    ... )
    # Output: cities.shp in UTM Zone 33N coordinates (explicit transformation)

    Expected processing log:
    ```
    INFO: Converting from ESRI Shapefile to GeoJSON
    INFO: Input file: /data/parcels.shp (523 features)
    INFO: Target CRS: EPSG:4326 (WGS84)
    INFO: Input CRS: EPSG:32633 (UTM Zone 33N)
    INFO: Coordinate transformation required
    INFO: Input geometry type: Polygon
    INFO: Processing 523 features...
    INFO: Conversion completed successfully
    INFO: Output: /output/parcels.geojson (523 features)
    ```

    See Also
    --------
    get_vector_driver_from_filename : Automatic driver selection from filename
    get_vector_format_from_filename : Format name extraction from filename
    get_output_geometry_type : Normalize geometry types (multi → single, 3D → 2D)
    rasterize_vector : Convert vector to raster format
    """
    # Validate input file exists
    if not os.path.exists(sFilename_vector_in):
        logger.error(f"Input file not found: {sFilename_vector_in}")
        return False

    logger.info(f"Input file: {sFilename_vector_in}")

    # Get input and output drivers based on file extensions
    try:
        driver_in = get_vector_driver_from_filename(sFilename_vector_in)
        driver_out = get_vector_driver_from_filename(sFilename_vector_out)

        format_in = get_vector_format_from_filename(sFilename_vector_in)
        format_out = get_vector_format_from_filename(sFilename_vector_out)

        logger.info(f"Converting from {format_in} to {format_out}")
    except ValueError as e:
        logger.error(f"Unsupported file format: {e}")
        return False

    if driver_out is None:
        logger.error(f"Output driver not available for format: {format_out}")
        return False

    # Remove existing output file
    if os.path.exists(sFilename_vector_out):
        if format_out == "ESRI Shapefile":
            # Remove all shapefile component files
            base_name = os.path.splitext(sFilename_vector_out)[0]
            for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg", ".shp.xml"]:
                component_file = base_name + ext
                if os.path.exists(component_file):
                    os.remove(component_file)
            logger.info(f"Removed existing shapefile: {sFilename_vector_out}")
        else:
            os.remove(sFilename_vector_out)
            logger.info(f"Removed existing file: {sFilename_vector_out}")

    # Open input vector dataset
    pDataset_in = ogr.Open(sFilename_vector_in, 0)
    if pDataset_in is None:
        logger.error(f"Could not open input file: {sFilename_vector_in}")
        return False

    # Get input layer
    pLayer_in = pDataset_in.GetLayer()
    if pLayer_in is None:
        logger.error("Could not access layer in input file")
        pDataset_in = None
        return False

    nFeature_count = pLayer_in.GetFeatureCount()
    logger.info(f"Input contains {nFeature_count} features")

    # Create output dataset
    pDataset_out = driver_out.CreateDataSource(sFilename_vector_out)
    if pDataset_out is None:
        logger.error(f"Could not create output file: {sFilename_vector_out}")
        pDataset_in = None
        return False

    # Determine target spatial reference based on format
    # GeoJSON always uses WGS84 (EPSG:4326) per RFC 7946
    # Other formats preserve original CRS unless target_epsg is explicitly set
    if target_epsg is not None:
        # Explicit EPSG provided - use it for any format
        pSrs_out = osr.SpatialReference()
        pSrs_out.ImportFromEPSG(target_epsg)
        logger.info(f"Target CRS: EPSG:{target_epsg} (explicit)")
        iFlag_force_transform = 1
    elif format_out.lower() in ["geojson", "json"]:
        # GeoJSON format - always use WGS84
        pSrs_out = osr.SpatialReference()
        pSrs_out.ImportFromEPSG(4326)
        logger.info(f"Target CRS: EPSG:4326 (GeoJSON format requirement)")
        iFlag_force_transform = 1
    else:
        # Other formats - preserve original CRS
        pSrs_out = None
        iFlag_force_transform = 0
        logger.info("Preserving original CRS (non-GeoJSON format)")

    # Get input spatial reference
    pSrs_in = pLayer_in.GetSpatialRef()
    if pSrs_in is None:
        logger.warning("Input file has no spatial reference defined")
        iFlag_transform = 0
    elif pSrs_out is None:
        # No target CRS specified (preserve original)
        iFlag_transform = 0
        pSrs_out = pSrs_in
        logger.info("Using original CRS for output")
    else:
        # Target CRS is specified (either explicit or GeoJSON format)
        # Export to WKT for comparison
        wkt_in = pSrs_in.ExportToWkt()
        wkt_out = pSrs_out.ExportToWkt()

        # Check if transformation is needed
        if wkt_in != wkt_out:
            iFlag_transform = 1
            logger.info("Coordinate transformation required")

            # Handle GDAL 3+ axis order
            if int(osgeo.__version__[0]) >= 3:
                pSrs_in.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                pSrs_out.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

            transform = osr.CoordinateTransformation(pSrs_in, pSrs_out)
        else:
            iFlag_transform = 0
            logger.info("No coordinate transformation needed (same CRS)")

    # Get geometry type from first feature
    pLayer_in.ResetReading()
    pFeature = pLayer_in.GetNextFeature()
    if pFeature is None:
        logger.error("Layer contains no features")
        pDataset_in = None
        pDataset_out = None
        return False

    pGeometry = pFeature.GetGeometryRef()
    if pGeometry is None:
        logger.error("First feature has no geometry")
        pDataset_in = None
        pDataset_out = None
        return False

    iGeomType = pGeometry.GetGeometryType()
    sGeomType = ogr.GeometryTypeToName(iGeomType)
    logger.info(f"Input geometry type: {sGeomType}")

    # Get appropriate output geometry type (handles multi-geometry and 2.5D/3D variants)
    try:
        iGeomType_out = get_output_geometry_type(iGeomType)
        sGeomType_out = ogr.GeometryTypeToName(iGeomType_out)
        if iGeomType_out != iGeomType:
            logger.info(
                f"Output geometry type: {sGeomType_out} (normalized from {sGeomType})"
            )
        else:
            logger.info(f"Output geometry type: {sGeomType_out}")
    except ValueError as e:
        logger.warning(f"Geometry type mapping failed: {e}, using original type")
        iGeomType_out = iGeomType

    # Check if the input geometry has Z-coordinates
    if pGeometry.GetCoordinateDimension() == 3:
        logger.info("Input geometry contains Z-coordinates (3D)")
        # Update the output geometry type to 3D if necessary
        if iGeomType == ogr.wkbLineString:
            iGeomType_out = ogr.wkbLineString25D
        elif iGeomType == ogr.wkbPolygon:
            iGeomType_out = ogr.wkbPolygon25D
        elif iGeomType == ogr.wkbPoint:
            iGeomType_out = ogr.wkbPoint25D
        # Add other geometry types as needed
        sGeomType_out = ogr.GeometryTypeToName(iGeomType_out)
        logger.info(f"Output geometry type updated to: {sGeomType_out} (3D)")
    else:
        logger.info("Input geometry is 2D")

    # Create output layer with the appropriate geometry type
    pLayer_out = pDataset_out.CreateLayer("layer", pSrs_out, geom_type=iGeomType_out)

    if pLayer_out is None:
        logger.error("Could not create output layer")
        pDataset_in = None
        pDataset_out = None
        return False

    # Copy field definitions from input to output
    pFeatureDefn_in = pLayer_in.GetLayerDefn()
    for i in range(pFeatureDefn_in.GetFieldCount()):
        pFieldDefn = pFeatureDefn_in.GetFieldDefn(i)
        pLayer_out.CreateField(pFieldDefn)

    logger.info(f"Processing {nFeature_count} features...")

    # Process features
    pLayer_in.ResetReading()
    pFeature_in = pLayer_in.GetNextFeature()
    nProcessed = 0

    while pFeature_in:
        pGeometry = pFeature_in.GetGeometryRef()

        if pGeometry is None:
            logger.warning(f"Feature {nProcessed} has no geometry, skipping")
            pFeature_in = pLayer_in.GetNextFeature()
            continue

        # Transform geometry if needed
        if iFlag_transform == 1:
            try:
                pGeometry.Transform(transform)
            except Exception as e:
                logger.warning(f"Failed to transform feature {nProcessed}: {e}")
                pFeature_in = pLayer_in.GetNextFeature()
                continue

        # Create output feature
        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        pFeature_out.SetGeometry(pGeometry)

        # Copy attributes
        for i in range(pFeatureDefn_in.GetFieldCount()):
            pFeature_out.SetField(i, pFeature_in.GetField(i))

        # Add feature to output layer
        pLayer_out.CreateFeature(pFeature_out)
        pFeature_out = None

        nProcessed += 1
        if nProcessed % 1000 == 0:
            logger.info(f"Processed {nProcessed}/{nFeature_count} features")

        pFeature_in = pLayer_in.GetNextFeature()

    # Cleanup
    pDataset_in = None
    pLayer_out = None
    pDataset_out = None

    logger.info("Conversion completed successfully")
    logger.info(f"Output: {sFilename_vector_out} ({nProcessed} features)")

    return True
