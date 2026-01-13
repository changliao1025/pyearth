"""
Calculate Geometric Difference Between Polygon Datasets
=========================================================

This module provides functionality to compute the geometric difference between two
polygon datasets, identifying areas in the base dataset that are not covered by
a new dataset. The operation uses spatial indexing for efficient processing of
large datasets.

The difference operation is commonly used in geospatial workflows for:

* **Change Detection**: Identify areas that have been removed or changed between
  time periods (e.g., deforestation, urban development reversal)
* **Gap Analysis**: Find coverage gaps when comparing datasets from different sources
* **Data Quality**: Validate dataset completeness by finding missing areas
* **Update Tracking**: Monitor what portions of a base dataset are not covered by
  updates or new data

Key Features
------------
- Supports multiple vector formats (GeoJSON, Shapefile, GeoPackage, KML, etc.)
- Automatic format detection based on file extension
- Efficient spatial indexing with rtree
- Geometry validation and automatic repair
- Multi-part geometry handling (MULTIPOLYGON, GEOMETRYCOLLECTION)
- Progress logging for long-running operations
- Comprehensive error handling and reporting
- Tracks source features in output (base_fid, new_fid)

Technical Details
-----------------
**Spatial Indexing**: The function uses `setup_spatial_index()` to select rtree for spatial indexing.

**Difference Operation**: For each feature in the new dataset:
1. Uses spatial index to find potentially intersecting base features (bounding box test)
2. Performs exact geometry intersection test
3. Calculates geometric difference (base geometry NOT in new geometry)
4. Outputs resulting difference polygons

**Geometry Handling**:
- Invalid geometries are automatically repaired using Buffer(0)
- Multi-part results are decomposed into individual polygons
- Empty results are filtered out

**Performance**: Spatial indexing provides O(log n) lookup time for intersection
candidates, making this suitable for large datasets with millions of features.

See Also
--------
calculate_polygon_intersection : Find overlapping areas between datasets
difference_polyline_with_polygon : Calculate differences for line geometries
spatial_join : Join features based on spatial relationships

Notes
-----
The difference operation is NOT symmetric: difference(A, B) ≠ difference(B, A).
This function returns areas in the base dataset that are NOT in the new dataset.

Examples
--------
Detect deforestation (areas in old forest layer not in new forest layer):

    >>> calculate_polygon_difference(
    ...     'forest_2020.geojson',
    ...     'forest_2024.geojson',
    ...     'deforestation.gpkg',
    ...     spatial_ref_epsg=4326
    ... )

Find coverage gaps between datasets:

    >>> calculate_polygon_difference(
    ...     'complete_coverage.shp',
    ...     'partial_coverage.geojson',
    ...     'coverage_gaps.parquet'
    ... )
"""

import os
import sys
import logging
from typing import Optional, Union, Tuple, Dict
from osgeo import ogr, osr
from datetime import datetime
from rtree.index import Index as RTreeindex

from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_driver_from_extension,
    get_vector_format_from_extension,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def difference_polygon_with_polygon_file(
    sFilename_base: str,
    sFilename_new: str,
    sFilename_difference_out: str,
    spatial_ref_epsg: int = 4326,
) -> None:
    """
    Calculate the geometric difference between base polygons and new polygons.

    This function identifies areas in the base polygon dataset that are NOT covered
    by the new polygon dataset. The operation uses spatial indexing for efficient
    processing and supports automatic format detection for input/output files.

    The difference is computed as: base_geometry - new_geometry, which returns the
    portion of the base geometry that does not overlap with the new geometry. This
    is useful for change detection, gap analysis, and data quality validation.

    Parameters
    ----------
    sFilename_base : str
        Path to the base polygon file (the dataset to subtract FROM).
        Supports multiple formats: .geojson, .shp, .gpkg, .parquet, .kml, etc.
        Format is automatically detected from file extension.

    sFilename_new : str
        Path to the new polygon file (the dataset to subtract).
        Supports multiple formats with automatic detection.
        Areas covered by this dataset will be removed from base features.

    sFilename_difference_out : str
        Path for the output difference file.
        Output format is determined by file extension.
        Contains polygons representing base areas not covered by new features.

    spatial_ref_epsg : int, optional
        EPSG code for output spatial reference system.
        Default is 4326 (WGS84 lat/lon).
        Common values: 4326 (WGS84), 3857 (Web Mercator), 32633 (UTM Zone 33N), etc.

    Returns
    -------
    None
        Results are written directly to the output file.
        The function logs progress and summary statistics.

    Raises
    ------
    FileNotFoundError
        If input files don't exist or cannot be accessed.
    ImportError
        If required spatial indexing library (rtree) is not available.
    RuntimeError
        If GDAL/OGR operations fail (invalid datasets, layer access errors, etc.).
    ValueError
        If file format is not supported by GDAL/OGR.

    Notes
    -----
     1. **Spatial Indexing**: Uses rtree (libspatialindex backend) for efficient candidate lookup.
         See `setup_spatial_index()` for details.

    2. **Geometry Validation**: Invalid geometries are automatically repaired using
       the Buffer(0) technique. Features with empty or null geometries are skipped
       with warnings logged.

    3. **Multi-part Handling**: MULTIPOLYGON and GEOMETRYCOLLECTION results are
       decomposed into individual POLYGON features for cleaner output.

    4. **Output Attributes**:
       - id: Sequential polygon identifier
       - base_fid: Feature ID from base dataset (for traceability)
       - new_fid: Feature ID from new dataset that intersected

    5. **Performance**: Uses spatial indexing for O(log n) candidate lookup.
       Processing time scales well with dataset size. Progress is logged every
       100 features for long-running operations.

    6. **Memory**: Base features are cached in memory for repeated access.
       For very large datasets (millions of features), consider splitting into
       tiles or using on-disk spatial indexes.

    7. **Coordinate Systems**: Input files can be in different coordinate systems,
       but reprojection is NOT performed automatically. For accurate results,
       ensure inputs are in the same CRS or reproject before calling this function.

    8. **Empty Results**: If no differences are found (new dataset completely
       covers base dataset), an empty output file with correct schema is created.

    9. **Shapefile Cleanup**: When overwriting shapefiles, all associated files
       (.shp, .shx, .dbf, .prj, .cpg) are automatically removed.

    10. **Logging**: Progress and errors are logged using Python's logging module.
        Set logging level to DEBUG for detailed operation information.

    Examples
    --------
    Detect deforestation (areas in 2020 forest not in 2024 forest):

        >>> calculate_polygon_difference(
        ...     sFilename_base='forest_2020.geojson',
        ...     sFilename_new='forest_2024.geojson',
        ...     sFilename_difference_out='deforestation_2020_2024.gpkg',
        ...     spatial_ref_epsg=4326
        ... )
        INFO:root:Starting polygon difference calculation at 2024-10-14 12:00:00
    INFO:root:Using rtree for spatial indexing
        INFO:root:Base file contains 15420 features
        INFO:root:New file contains 14856 features
        INFO:root:Indexed 15420 valid base features
        INFO:root:Generated 2341 difference polygons
        INFO:root:Processing completed in 45.23 seconds (0.75 minutes)

    Find coverage gaps between complete and partial datasets:

        >>> calculate_polygon_difference(
        ...     sFilename_base='complete_coverage.shp',
        ...     sFilename_new='partial_coverage.geojson',
        ...     sFilename_difference_out='coverage_gaps.parquet'
        ... )

    Track urban development reversals (areas developed in 2020 but not in 2024):

        >>> calculate_polygon_difference(
        ...     sFilename_base='urban_2020.gpkg',
        ...     sFilename_new='urban_2024.gpkg',
        ...     sFilename_difference_out='urban_reversals.geojson',
        ...     spatial_ref_epsg=3857  # Web Mercator
        ... )

    Quality control - find areas missing from updated dataset:

        >>> calculate_polygon_difference(
        ...     sFilename_base='original_parcels.shp',
        ...     sFilename_new='updated_parcels.shp',
        ...     sFilename_difference_out='missing_parcels.gpkg'
        ... )

    See Also
    --------
    calculate_polygon_intersection : Find overlapping areas between datasets
    ogr.Geometry.Difference : Low-level OGR difference operation
    setup_spatial_index : Spatial indexing library selection
    difference_polyline_with_polygon : Difference operation for line geometries
    """

    # Setup spatial indexing library (rtree only)
    logger.info(f"Using rtree for spatial indexing")

    # Validate that input files exist before attempting to open them
    # This provides clearer error messages than GDAL's generic failures
    if not os.path.exists(sFilename_base):
        raise FileNotFoundError(f"Base file not found: {sFilename_base}")
    if not os.path.exists(sFilename_new):
        raise FileNotFoundError(f"New file not found: {sFilename_new}")

    # Ensure output directory exists, create if necessary
    # This prevents failures when writing to non-existent directories
    output_dir = os.path.dirname(sFilename_difference_out)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Created output directory: {output_dir}")

    # Record start time for performance reporting
    start_time = datetime.now()
    logger.info(f"Starting polygon difference calculation at {start_time}")

    try:
        # Determine appropriate drivers for all files based on their extensions
        # Uses centralized utility functions for consistent format handling
        try:
            driver_base = get_vector_driver_from_extension(sFilename_base)
            driver_new = get_vector_driver_from_extension(sFilename_new)
            driver_out = get_vector_driver_from_extension(sFilename_difference_out)

            format_base = get_vector_format_from_extension(sFilename_base)
            format_new = get_vector_format_from_extension(sFilename_new)
            format_out = get_vector_format_from_extension(sFilename_difference_out)

            logger.info(
                f"Input formats - Base: {format_base}, New: {format_new}, Output: {format_out}"
            )
        except ValueError as e:
            raise RuntimeError(f"Unsupported file format: {e}")

        # Remove existing output file to ensure clean write
        # For shapefiles, must remove all component files (.shp, .shx, .dbf, etc.)
        if os.path.exists(sFilename_difference_out):
            if format_out == "ESRI Shapefile":
                # Shapefiles consist of multiple files with different extensions
                base_name = os.path.splitext(sFilename_difference_out)[0]
                for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg"]:
                    if os.path.exists(base_name + ext):
                        os.remove(base_name + ext)
            else:
                # Single-file formats (GeoJSON, GeoPackage, Parquet, etc.)
                os.remove(sFilename_difference_out)
            logger.info(f"Removed existing output file: {sFilename_difference_out}")

        # Create output dataset using the format-specific driver
        pDataset_out = driver_out.CreateDataSource(sFilename_difference_out)
        if pDataset_out is None:
            raise RuntimeError(
                f"Could not create output file: {sFilename_difference_out}"
            )

        # Configure spatial reference system for output layer
        # All output features will be in this coordinate system
        pSpatial_reference_gcs = osr.SpatialReference()
        pSpatial_reference_gcs.ImportFromEPSG(spatial_ref_epsg)
        logger.info(f"Using spatial reference EPSG:{spatial_ref_epsg}")

        # Create output layer with polygon geometry type
        # Layer name 'diff' is arbitrary but descriptive
        pLayerOut = pDataset_out.CreateLayer(
            "diff", pSpatial_reference_gcs, ogr.wkbPolygon
        )
        if pLayerOut is None:
            raise RuntimeError("Could not create output layer")

        # Define output attribute schema
        # id: Sequential identifier for output polygons
        # base_fid: Original feature ID from base dataset (for traceability)
        # new_fid: Feature ID from new dataset that intersected (for analysis)
        pLayerOut.CreateField(ogr.FieldDefn("id", ogr.OFTInteger64))
        pLayerOut.CreateField(ogr.FieldDefn("base_fid", ogr.OFTInteger64))
        pLayerOut.CreateField(ogr.FieldDefn("new_fid", ogr.OFTInteger64))

        pLayerDefn = pLayerOut.GetLayerDefn()

        # Open base dataset in read-only mode (0)
        pDataset_base = driver_base.Open(sFilename_base, 0)
        if pDataset_base is None:
            raise RuntimeError(f"Could not open base file: {sFilename_base}")

        # Access the first layer from base dataset
        pLayer_base = pDataset_base.GetLayer()
        if pLayer_base is None:
            raise RuntimeError("Could not access layer in base file")

        nFeature_base = pLayer_base.GetFeatureCount()
        if nFeature_base == 0:
            logger.warning("Base file contains no features")
            return

        logger.info(f"Base file contains {nFeature_base} features")

        # Open and validate new file
        pDataset_new = driver_new.Open(sFilename_new, 0)
        if pDataset_new is None:
            raise RuntimeError(f"Could not open new file: {sFilename_new}")

        pLayer_new = pDataset_new.GetLayer()
        if pLayer_new is None:
            raise RuntimeError("Could not access layer in new file")

        nFeature_new = pLayer_new.GetFeatureCount()
        if nFeature_new == 0:
            logger.warning("New file contains no features")
            return

        logger.info(f"New file contains {nFeature_new} features")

        # Build spatial index for efficient intersection queries
        # The index stores bounding boxes of base features for fast lookup
        logger.info("Building spatial index for base features...")

        # rtree: Use default configuration
        index_base = RTreeindex()

        # Cache base features in memory to avoid repeated disk access
        # Key: feature ID, Value: cloned feature object
        base_features = {}
        indexed_count = 0

        # Iterate through all base features to build the spatial index
        for idx, pFeature_base in enumerate(pLayer_base):
            fid = pFeature_base.GetFID()
            pGeometry_base = pFeature_base.GetGeometryRef()

            # Skip features without geometry
            if pGeometry_base is None:
                logger.warning(f"Base feature {fid} has no geometry, skipping")
                continue

            # Validate and repair invalid geometries using Buffer(0) technique
            # Buffer(0) rebuilds the geometry, often fixing self-intersections and other issues
            if not pGeometry_base.IsValid():
                logger.warning(
                    f"Base feature {fid} has invalid geometry, attempting to fix"
                )
                pGeometry_base = pGeometry_base.Buffer(0)

            # Skip features that are still invalid or empty after repair
            if pGeometry_base is None or pGeometry_base.IsEmpty():
                logger.warning(
                    f"Base feature {fid} has empty geometry after validation, skipping"
                )
                continue

            try:
                # Get the bounding box (envelope) of the geometry
                # Returns: (minX, maxX, minY, maxY)
                envelope = pGeometry_base.GetEnvelope()
                left, right, bottom, top = envelope

                # Insert bounding box into spatial index
                # rtree use coordinate order
                pBound = (left, bottom, right, top)
                index_base.insert(fid, pBound)

                # Clone and cache the feature for later retrieval
                # Clone() creates a deep copy to prevent issues with feature reuse
                base_features[fid] = pFeature_base.Clone()
                indexed_count += 1
            except Exception as e:
                logger.error(f"Error indexing base feature {fid}: {e}")
                continue

        logger.info(f"Indexed {indexed_count} valid base features")

        # Exit early if no valid features were indexed
        if not base_features:
            logger.warning("No valid base features found for indexing")
            return

        # Process new features and calculate geometric differences
        # For each new feature, find overlapping base features and compute difference
        logger.info("Processing differences...")
        lID_polygon = 1  # Sequential ID for output polygons
        processed_features = 0

        # Iterate through all new features
        for idx, pFeature_new in enumerate(pLayer_new):
            processed_features += 1
            # Log progress every 100 features for long-running operations
            if processed_features % 100 == 0:
                logger.info(
                    f"Processed {processed_features}/{nFeature_new} new features"
                )

            new_fid = pFeature_new.GetFID()
            pGeometry_new = pFeature_new.GetGeometryRef()

            # Skip features without geometry
            if pGeometry_new is None:
                logger.warning(f"New feature {new_fid} has no geometry, skipping")
                continue

            # Validate and repair invalid geometries
            if not pGeometry_new.IsValid():
                logger.warning(
                    f"New feature {new_fid} has invalid geometry, attempting to fix"
                )
                pGeometry_new = pGeometry_new.Buffer(0)

            # Skip features that are still invalid or empty after repair
            if pGeometry_new is None or pGeometry_new.IsEmpty():
                logger.warning(
                    f"New feature {new_fid} has empty geometry after validation, skipping"
                )
                continue

            try:
                # Get bounding box of new feature for spatial index query
                envelope = pGeometry_new.GetEnvelope()
                left, right, bottom, top = envelope

                # Query spatial index to find base features with overlapping bounding boxes
                # This is a fast filter that returns candidates (may include false positives)
                aIntersect = list(index_base.intersection((left, bottom, right, top)))

                # Process each candidate base feature
                for base_fid in aIntersect:
                    # Verify feature is in cache (should always be true)
                    if base_fid not in base_features:
                        continue

                    pFeature_base = base_features[base_fid]
                    pGeometry_base = pFeature_base.GetGeometryRef()

                    try:
                        # Perform exact geometric intersection test
                        # Bounding box overlap doesn't guarantee geometry overlap
                        if not pGeometry_new.Intersects(pGeometry_base):
                            continue

                        # Calculate geometric difference: base - new
                        # Returns the portion of base geometry NOT covered by new geometry
                        pGeometry_diff = pGeometry_base.Difference(pGeometry_new)

                        # Skip if difference is empty (new completely covers base)
                        if pGeometry_diff is None or pGeometry_diff.IsEmpty():
                            continue

                        # Handle output based on geometry type
                        pGeometrytype_intersect = pGeometry_diff.GetGeometryName()

                        if pGeometrytype_intersect == "POLYGON":
                            # Simple polygon result - create single output feature
                            pFeatureOut = ogr.Feature(pLayerDefn)
                            pFeatureOut.SetGeometry(pGeometry_diff)
                            pFeatureOut.SetField("id", lID_polygon)
                            pFeatureOut.SetField("base_fid", base_fid)
                            pFeatureOut.SetField("new_fid", new_fid)
                            pLayerOut.CreateFeature(pFeatureOut)
                            lID_polygon += 1

                        elif pGeometrytype_intersect in [
                            "MULTIPOLYGON",
                            "GEOMETRYCOLLECTION",
                        ]:
                            # Multi-part geometry result - decompose into individual polygons
                            # This happens when new feature splits base feature into multiple parts
                            for i in range(pGeometry_diff.GetGeometryCount()):
                                pGeometry_part = pGeometry_diff.GetGeometryRef(i)
                                if (
                                    pGeometry_part is not None
                                    and not pGeometry_part.IsEmpty()
                                ):
                                    part_type = pGeometry_part.GetGeometryName()
                                    # Only extract polygon parts (ignore points, lines from collection)
                                    if part_type == "POLYGON":
                                        pFeatureOut = ogr.Feature(pLayerDefn)
                                        pFeatureOut.SetGeometry(pGeometry_part)
                                        pFeatureOut.SetField("id", lID_polygon)
                                        pFeatureOut.SetField("base_fid", base_fid)
                                        pFeatureOut.SetField("new_fid", new_fid)
                                        pLayerOut.CreateFeature(pFeatureOut)
                                        lID_polygon += 1

                    except Exception as e:
                        logger.error(
                            f"Error processing intersection between base {base_fid} and new {new_fid}: {e}"
                        )
                        continue

            except Exception as e:
                logger.error(f"Error processing new feature {new_fid}: {e}")
                continue

        logger.info(f"Generated {lID_polygon - 1} difference polygons")

    except Exception as e:
        logger.error(f"Error during processing: {e}")
        raise RuntimeError(f"Failed to calculate polygon difference: {e}")

    finally:
        # Clean up resources
        if "pDataset_base" in locals():
            pDataset_base = None
        if "pDataset_new" in locals():
            pDataset_new = None
        if "pDataset_out" in locals():
            pDataset_out = None

    # Report timing
    end_time = datetime.now()
    delta = end_time - start_time
    seconds = delta.total_seconds()
    minutes = seconds / 60

    logger.info(
        f"Processing completed in {seconds:.2f} seconds ({minutes:.2f} minutes)"
    )

    return None
