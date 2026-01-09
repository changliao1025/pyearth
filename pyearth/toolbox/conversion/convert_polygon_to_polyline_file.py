"""
Polygon to Polyline Conversion

This module provides functionality for converting polygon geometries to polyline
(linestring) geometries by extracting the exterior ring boundaries. The conversion
supports multiple vector formats with conditional coordinate transformation based
on output format requirements.

Main Functions
--------------
convert_polygon_to_polyline_file : Convert polygon boundaries to polylines

Key Features
------------
- Multi-format support for input and output (auto-detected from extension)
- Extracts exterior ring from polygon geometries
- Conditional coordinate transformation (GeoJSON → WGS84, others preserve CRS)
- Handles both simple polygons and multi-polygons
- Preserves feature IDs and optional attributes
- Comprehensive error handling and logging
- Progress tracking for large datasets

Use Cases
---------
1. **Boundary Extraction**: Convert area features to their perimeter lines
2. **Network Analysis**: Create road/river networks from watershed polygons
3. **Cartographic Simplification**: Show boundaries without fill for clarity
4. **Data Exchange**: Convert polygon data for line-based analysis
5. **Web Visualization**: Extract building footprints as outlines for web maps
6. **Topological Analysis**: Analyze shared boundaries between adjacent polygons

Technical Details
-----------------
The module extracts the exterior ring (boundary) from each polygon geometry and
creates a new linestring feature. For multi-polygon geometries, each polygon's
exterior ring is extracted separately, creating multiple linestring features.

Coordinate transformation follows the same rules as convert_vector_format:
- GeoJSON output: Automatically transforms to WGS84 (EPSG:4326)
- Other formats: Preserve original coordinate system
- Can be overridden with explicit target_epsg parameter

Performance Characteristics
---------------------------
- Time Complexity: O(N × M) where N = features, M = avg points per boundary
- Space Complexity: O(1) - streaming feature-by-feature
- Memory efficient for large polygon datasets

Dependencies
------------
- GDAL/OGR: Vector I/O and geometry operations
- OSR: Spatial reference system handling

See Also
--------
- get_vector_driver_from_extension: Automatic driver selection
- get_vector_format_from_extension: Format name extraction
- convert_vector_format: General format conversion with transformation
"""

import os
import logging
from typing import Optional
import osgeo
from osgeo import ogr, osr

from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_driver_from_extension,
    get_vector_format_from_extension,
)

# Configure logging
logger = logging.getLogger(__name__)


def convert_polygon_to_polyline_file(
    sFilename_polygon_in: str,
    sFilename_polyline_out: str,
    target_epsg: Optional[int] = None,
    iFlag_copy_attributes: int = 0,
) -> bool:
    """
    Convert polygon geometries to polyline geometries by extracting exterior rings.

    This function reads a polygon vector file and creates a polyline file by
    extracting the exterior ring (boundary) from each polygon. Multi-polygon
    features are decomposed into separate linestring features. The output format
    is automatically detected from the file extension.

    Coordinate transformation behavior:
    - GeoJSON output: Always transforms to WGS84 (EPSG:4326) per RFC 7946
    - Other formats: Preserve original CRS unless target_epsg is specified
    - Explicit override: Setting target_epsg forces transformation for any format

    Parameters
    ----------
    sFilename_polygon_in : str
        Absolute path to input polygon vector file. Format automatically detected
        from extension. Supports: Shapefile (.shp), GeoJSON (.geojson), GeoPackage
        (.gpkg), and all other GDAL-supported polygon formats.
    sFilename_polyline_out : str
        Absolute path for output polyline vector file. Format automatically detected
        from extension. Existing files will be overwritten. Common extensions:
        - .geojson or .json for GeoJSON (auto-transforms to WGS84)
        - .shp for Shapefile (preserves original CRS)
        - .gpkg for GeoPackage (preserves original CRS)
    target_epsg : Optional[int], optional
        EPSG code for output coordinate system. If None (default), behavior depends
        on output format:
        - GeoJSON output: Uses EPSG:4326 (WGS84)
        - Other formats: Preserve original input CRS
        If specified, forces transformation to this EPSG code regardless of format.
    iFlag_copy_attributes : int, optional
        Whether to copy attribute fields from input polygons to output polylines.
        - 0 (default): Only copy feature ID
        - 1: Copy all attribute fields from input layer
        Note: Copying attributes may increase processing time for large datasets.

    Returns
    -------
    bool
        True if conversion completed successfully, False if errors occurred.

    Raises
    ------
    FileNotFoundError
        If input file does not exist.
    ValueError
        If input or output format is not supported.
    RuntimeError
        If file cannot be opened, geometry extraction fails, or conversion fails.

    Notes
    -----
    1. **Geometry Type Handling**: Only processes Polygon and MultiPolygon geometries.
       For MultiPolygon features, each constituent polygon's exterior ring is
       extracted as a separate linestring feature. Other geometry types are skipped
       with a warning.

    2. **Ring Extraction**: Only the exterior ring (boundary) is extracted from
       each polygon. Interior rings (holes) are not included in the output. To
       preserve holes, use a different conversion approach.

    3. **Coordinate Transformation**: Follows the same logic as convert_vector_format:
       - GeoJSON output always uses WGS84 for web compatibility
       - Other formats preserve projection for GIS workflows
       - Explicit target_epsg overrides default behavior

    4. **Feature ID Preservation**: The original feature ID (FID) from polygon
       features is preserved in the output polyline features via an 'id' field.
       This allows tracking the source polygon for each polyline.

    5. **Attribute Handling**: When iFlag_copy_attributes=1, all attribute fields
       from the input layer are copied to the output layer with their original
       names, types, and values. This is useful when polygon attributes are
       meaningful for the boundary lines.

    6. **GDAL 3+ Compatibility**: For GDAL version 3 and above, axis mapping is
       set to traditional GIS order (longitude, latitude) to maintain consistency
       with earlier versions.

    7. **Memory Efficiency**: Features are processed one at a time in streaming
       fashion, making the function suitable for large polygon datasets without
       excessive memory usage.

    8. **File Cleanup**: Existing output files are automatically removed before
       creating new output. For Shapefiles, all component files (.shp, .shx, .dbf,
       .prj) are removed.

    9. **Multi-Format Support**: Both input and output formats are auto-detected
       using get_vector_driver_from_extension(). No need to specify drivers manually.

    10. **Progress Tracking**: For datasets with more than 100 features, progress
        is logged every 1000 features to provide feedback during long operations.

    Examples
    --------
    Convert polygon shapefile to polyline GeoJSON (automatic WGS84):

    >>> success = convert_polygon_to_polyline_file(
    ...     sFilename_polygon_in='/data/watersheds.shp',
    ...     sFilename_polyline_out='/output/boundaries.geojson'
    ... )
    >>> print(f"Success: {success}")
    # Output: boundaries.geojson with watershed boundary lines in WGS84

    Convert polygon GeoJSON to polyline shapefile (preserve CRS):

    >>> success = convert_polygon_to_polyline_file(
    ...     sFilename_polygon_in='/data/parcels.geojson',
    ...     sFilename_polyline_out='/output/parcel_lines.shp'
    ... )
    # Output: parcel_lines.shp with original CRS preserved

    Convert with attribute preservation and explicit projection:

    >>> success = convert_polygon_to_polyline_file(
    ...     sFilename_polygon_in='/data/admin_boundaries.gpkg',
    ...     sFilename_polyline_out='/output/admin_lines.gpkg',
    ...     target_epsg=3857,  # Web Mercator
    ...     iFlag_copy_attributes=1
    ... )
    # Output: admin_lines.gpkg in Web Mercator with all attributes copied

    Expected processing log:
    ```
    INFO: Input file: /data/watersheds.shp (2,345 features)
    INFO: Converting from ESRI Shapefile to GeoJSON
    INFO: Target CRS: EPSG:4326 (GeoJSON format requirement)
    INFO: Coordinate transformation required
    INFO: Processing 2,345 polygon features...
    INFO: Processed 1,000 features
    INFO: Processed 2,000 features
    INFO: Extracted 2,345 polylines from 2,345 polygons
    INFO: Conversion completed successfully
    INFO: Output: /output/boundaries.geojson (2,345 features)
    ```

    See Also
    --------
    get_vector_driver_from_extension : Automatic driver selection from filename
    get_vector_format_from_extension : Format name extraction from filename
    convert_vector_format : General vector format conversion
    """
    # Validate input file exists
    if not os.path.exists(sFilename_polygon_in):
        logger.error(f"Input file not found: {sFilename_polygon_in}")
        return False

    logger.info(f"Input file: {sFilename_polygon_in}")

    # Get input and output drivers based on file extensions
    try:
        driver_in = get_vector_driver_from_extension(sFilename_polygon_in)
        driver_out = get_vector_driver_from_extension(sFilename_polyline_out)

        format_in = get_vector_format_from_extension(sFilename_polygon_in)
        format_out = get_vector_format_from_extension(sFilename_polyline_out)

        logger.info(f"Converting from {format_in} to {format_out}")
    except ValueError as e:
        logger.error(f"Unsupported file format: {e}")
        return False

    if driver_out is None:
        logger.error(f"Output driver not available for format: {format_out}")
        return False

    # Remove existing output file
    if os.path.exists(sFilename_polyline_out):
        if format_out == "ESRI Shapefile":
            # Remove all shapefile component files
            base_name = os.path.splitext(sFilename_polyline_out)[0]
            for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg", ".shp.xml"]:
                component_file = base_name + ext
                if os.path.exists(component_file):
                    os.remove(component_file)
            logger.info(f"Removed existing shapefile: {sFilename_polyline_out}")
        else:
            os.remove(sFilename_polyline_out)
            logger.info(f"Removed existing file: {sFilename_polyline_out}")

    # Open input polygon dataset
    pDataset_in = ogr.Open(sFilename_polygon_in, 0)
    if pDataset_in is None:
        logger.error(f"Could not open input file: {sFilename_polygon_in}")
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
    pDataset_out = driver_out.CreateDataSource(sFilename_polyline_out)
    if pDataset_out is None:
        logger.error(f"Could not create output file: {sFilename_polyline_out}")
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
        if pSrs_out is None:
            pSrs_out = osr.SpatialReference()
            pSrs_out.ImportFromEPSG(4326)
            logger.warning("Using WGS84 as default for output")
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

    # Create output layer for polylines
    pLayer_out = pDataset_out.CreateLayer(
        "polyline", pSrs_out, geom_type=ogr.wkbLineString
    )
    if pLayer_out is None:
        logger.error("Could not create output layer")
        pDataset_in = None
        pDataset_out = None
        return False

    # Add ID field to output layer
    pLayer_out.CreateField(ogr.FieldDefn("id", ogr.OFTInteger64))

    # Copy field definitions from input if requested
    if iFlag_copy_attributes == 1:
        pFeatureDefn_in = pLayer_in.GetLayerDefn()
        nField_count = pFeatureDefn_in.GetFieldCount()
        for i in range(nField_count):
            pFieldDefn = pFeatureDefn_in.GetFieldDefn(i)
            pLayer_out.CreateField(pFieldDefn)
        logger.info(f"Copied {nField_count} attribute fields to output")

    # Get layer definition for creating features
    pLayerDefn_out = pLayer_out.GetLayerDefn()

    logger.info(f"Processing {nFeature_count} polygon features...")

    # Process polygon features
    pLayer_in.ResetReading()
    nProcessed = 0
    nPolylines_created = 0
    nSkipped = 0

    for pFeature_in in pLayer_in:
        fid = pFeature_in.GetFID()
        pGeometry_in = pFeature_in.GetGeometryRef()

        if pGeometry_in is None:
            logger.warning(f"Feature {fid} has no geometry, skipping")
            nSkipped += 1
            continue

        geom_type = pGeometry_in.GetGeometryType()

        # Handle Polygon geometry
        if geom_type == ogr.wkbPolygon or geom_type == ogr.wkbPolygon25D:
            # Get the exterior ring
            pRing = pGeometry_in.GetGeometryRef(0)
            if pRing is None:
                logger.warning(f"Feature {fid}: Could not extract exterior ring")
                nSkipped += 1
                continue

            # Create linestring from exterior ring
            pLineString = ogr.Geometry(ogr.wkbLineString)
            for i in range(pRing.GetPointCount()):
                point = pRing.GetPoint(i)
                pLineString.AddPoint(point[0], point[1])

            # Transform if needed
            if iFlag_transform == 1:
                try:
                    pLineString.Transform(transform)
                except Exception as e:
                    logger.warning(f"Failed to transform feature {fid}: {e}")
                    nSkipped += 1
                    continue

            # Create output feature
            pFeature_out = ogr.Feature(pLayerDefn_out)
            pFeature_out.SetGeometry(pLineString)
            pFeature_out.SetField("id", fid)

            # Copy attributes if requested
            if iFlag_copy_attributes == 1:
                pFeatureDefn_in = pLayer_in.GetLayerDefn()
                for i in range(pFeatureDefn_in.GetFieldCount()):
                    field_name = pFeatureDefn_in.GetFieldDefn(i).GetName()
                    pFeature_out.SetField(field_name, pFeature_in.GetField(i))

            # Add feature to output layer
            pLayer_out.CreateFeature(pFeature_out)
            pFeature_out = None
            nPolylines_created += 1

        # Handle MultiPolygon geometry
        elif geom_type == ogr.wkbMultiPolygon or geom_type == ogr.wkbMultiPolygon25D:
            nPolygon_count = pGeometry_in.GetGeometryCount()

            for j in range(nPolygon_count):
                pPolygon = pGeometry_in.GetGeometryRef(j)
                pRing = pPolygon.GetGeometryRef(0)

                if pRing is None:
                    logger.warning(
                        f"Feature {fid}, polygon {j}: Could not extract exterior ring"
                    )
                    continue

                # Create linestring from exterior ring
                pLineString = ogr.Geometry(ogr.wkbLineString)
                for i in range(pRing.GetPointCount()):
                    point = pRing.GetPoint(i)
                    pLineString.AddPoint(point[0], point[1])

                # Transform if needed
                if iFlag_transform == 1:
                    try:
                        pLineString.Transform(transform)
                    except Exception as e:
                        logger.warning(
                            f"Failed to transform feature {fid}, polygon {j}: {e}"
                        )
                        continue

                # Create output feature
                pFeature_out = ogr.Feature(pLayerDefn_out)
                pFeature_out.SetGeometry(pLineString)
                pFeature_out.SetField("id", fid)

                # Copy attributes if requested
                if iFlag_copy_attributes == 1:
                    pFeatureDefn_in = pLayer_in.GetLayerDefn()
                    for i in range(pFeatureDefn_in.GetFieldCount()):
                        field_name = pFeatureDefn_in.GetFieldDefn(i).GetName()
                        pFeature_out.SetField(field_name, pFeature_in.GetField(i))

                # Add feature to output layer
                pLayer_out.CreateFeature(pFeature_out)
                pFeature_out = None
                nPolylines_created += 1
        else:
            # Skip non-polygon geometries
            geom_name = ogr.GeometryTypeToName(geom_type)
            logger.warning(
                f"Feature {fid}: Skipping non-polygon geometry ({geom_name})"
            )
            nSkipped += 1
            continue

        nProcessed += 1
        if nProcessed % 1000 == 0 and nFeature_count > 100:
            logger.info(f"Processed {nProcessed}/{nFeature_count} features")

    # Cleanup
    pDataset_in = None
    pLayer_out = None
    pDataset_out = None

    logger.info(f"Extracted {nPolylines_created} polylines from {nProcessed} polygons")
    if nSkipped > 0:
        logger.warning(f"Skipped {nSkipped} features (non-polygon or invalid geometry)")

    logger.info("Conversion completed successfully")
    logger.info(f"Output: {sFilename_polyline_out} ({nPolylines_created} features)")

    return True
