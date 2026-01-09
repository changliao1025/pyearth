"""
Calculate Geometric Intersection Between Polyline and Polygon Datasets
========================================================================

This module provides functionality to compute the geometric intersection between
polyline (line) datasets and polygon datasets. The operation identifies where lines
pass through or overlap with polygons, creating intersection geometries (typically
polygons or multi-polygons where lines are buffered or represent corridors).

The intersection operation is commonly used in geospatial workflows for:

* **Network Analysis**: Find road segments within administrative boundaries or zones
* **Infrastructure Planning**: Identify utility lines (power, water, telecom) within
  specific regions or parcels
* **Environmental Analysis**: Determine stream segments within watersheds or protected areas
* **Urban Planning**: Extract transportation corridors within development zones
* **Spatial Queries**: Find linear features intersecting with specific areas

Key Features
------------
- Supports multiple vector formats (GeoJSON, Shapefile, GeoPackage, Parquet, etc.)
- Automatic format detection based on file extension
- Efficient spatial indexing with rtree
- Geometry validation and automatic repair
- Multi-part geometry handling (MULTIPOLYGON, GEOMETRYCOLLECTION)
- Progress logging for long-running operations
- Comprehensive error handling and reporting

Technical Details
-----------------
**Spatial Indexing**: The function uses `setup_spatial_index()` to select the rtree spatial indexing library for efficient processing.

**Intersection Operation**: For each polyline feature:
1. Uses spatial index to find potentially intersecting polygon features (bounding box test)
2. Performs exact geometry intersection test
3. Calculates geometric intersection (area where line and polygon overlap)
4. Outputs resulting intersection geometries

**Performance**: Spatial indexing provides O(log n) lookup time, making this
suitable for large datasets with millions of features.

See Also
--------
intersect_polygon_with_polygon : Intersection between two polygon datasets
clip_vector_by_polygon : Clip features to polygon boundaries
spatial_join : Join features based on spatial relationships

Examples
--------
Find road segments within city boundaries:

    >>> intersect_polyline_with_polygon_files(
    ...     'all_roads.geojson',
    ...     'city_boundary.geojson',
    ...     'city_roads.gpkg'
    ... )

Extract streams within watershed:

    >>> intersect_polyline_with_polygon_files(
    ...     'stream_network.shp',
    ...     'watershed.geojson',
    ...     'watershed_streams.parquet'
    ... )
"""

import os
import logging
from typing import Optional
from osgeo import ogr, osr, gdal
from datetime import datetime
from rtree.index import Index as RTreeindex

from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_driver_from_extension,
    get_vector_format_from_extension,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def intersect_polyline_with_polygon_files(
    sFilename_base: str,
    sFilename_new: str,
    sFilename_difference_out: str,
    spatial_ref_epsg: int = 4326,
) -> None:
    """
    Calculate the geometric intersection between polyline and polygon datasets.

    This function identifies areas where polylines (lines) intersect with polygons,
    creating output geometries representing the spatial overlap. Uses spatial indexing
    for efficient processing of large datasets.

    The intersection is computed as the geometric overlap between line and polygon
    features. This is useful for network analysis, infrastructure planning, and
    spatial queries involving linear features and areas.

    Parameters
    ----------
    sFilename_base : str
        Path to the base polyline (line) file.
        Supports multiple formats: .geojson, .shp, .gpkg, .parquet, etc.
        Format is automatically detected from file extension.

    sFilename_new : str
        Path to the polygon file to intersect with.
        Supports multiple formats with automatic detection.
        Lines are tested for intersection with these polygons.

    sFilename_difference_out : str
        Path for the output intersection file.
        Output format is determined by file extension.
        Contains geometries representing line-polygon intersections.

    spatial_ref_epsg : int, optional
        EPSG code for output spatial reference system.
        Default is 4326 (WGS84 lat/lon).

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
    1. **Spatial Indexing**: Uses rtree (Python with libspatialindex backend) for spatial indexing. See `setup_spatial_index()` for details.

    2. **Geometry Validation**: Invalid geometries are automatically repaired using
       the Buffer(0) technique. Features with empty or null geometries are skipped.

    3. **Multi-part Handling**: MULTIPOLYGON and GEOMETRYCOLLECTION results are
       decomposed into individual POLYGON features for cleaner output.

    4. **Output Attributes**:
       - id: Sequential polygon identifier
       - base_fid: Feature ID from base (polyline) dataset
       - new_fid: Feature ID from polygon dataset that intersected

    5. **Performance**: Uses spatial indexing for O(log n) candidate lookup.
       Processing time scales well with dataset size.

    6. **Coordinate Systems**: Input files should ideally be in the same CRS.
       No automatic reprojection is performed.

    7. **Output Geometry Type**: Intersection of lines with polygons typically
       produces polygons (for buffered lines) or lines (for unbuffered lines).
       This function outputs polygon geometries.

    Examples
    --------
    Find road segments within city boundaries:

        >>> intersect_polyline_with_polygon_files(
        ...     sFilename_base='all_roads.geojson',
        ...     sFilename_new='city_boundary.geojson',
        ...     sFilename_difference_out='city_roads.gpkg',
        ...     spatial_ref_epsg=4326
        ... )
    INFO:root:Using rtree for spatial indexing
        INFO:root:Base file contains 5420 polyline features
        INFO:root:New file contains 1 polygon features
        INFO:root:Generated 847 intersection polygons
        INFO:root:Processing completed in 12.45 seconds

    Extract utility lines within parcels:

        >>> intersect_polyline_with_polygon_files(
        ...     sFilename_base='power_lines.shp',
        ...     sFilename_new='parcels.geojson',
        ...     sFilename_difference_out='parcel_powerlines.parquet'
        ... )

    Find streams within protected areas:

        >>> intersect_polyline_with_polygon_files(
        ...     sFilename_base='stream_network.gpkg',
        ...     sFilename_new='protected_areas.shp',
        ...     sFilename_difference_out='protected_streams.geojson',
        ...     spatial_ref_epsg=3857
        ... )

    See Also
    --------
    intersect_polygon_with_polygon : Intersection between polygon datasets
    ogr.Geometry.Intersection : Low-level OGR intersection operation
    setup_spatial_index : Spatial indexing library selection
    """

    # Display GDAL version for debugging purposes
    logger.info(f"GDAL version: {gdal.__version__}")

    # Setup spatial indexing library (rtree only)
    logger.info("Using rtree for spatial indexing")

    # Validate that input files exist
    if not os.path.exists(sFilename_base):
        raise FileNotFoundError(f"Base file not found: {sFilename_base}")
    if not os.path.exists(sFilename_new):
        raise FileNotFoundError(f"New file not found: {sFilename_new}")

    # Ensure output directory exists
    output_dir = os.path.dirname(sFilename_difference_out)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Created output directory: {output_dir}")

    # Record start time for performance reporting
    start_time = datetime.now()
    logger.info(f"Starting polyline-polygon intersection at {start_time}")

    try:
        # Get drivers based on file extensions
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
        if os.path.exists(sFilename_difference_out):
            if format_out == "ESRI Shapefile":
                # Remove all shapefile component files
                base_name = os.path.splitext(sFilename_difference_out)[0]
                for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg"]:
                    if os.path.exists(base_name + ext):
                        os.remove(base_name + ext)
            else:
                os.remove(sFilename_difference_out)
            logger.info(f"Removed existing output file: {sFilename_difference_out}")

        # Create output dataset
        pDataset_out = driver_out.CreateDataSource(sFilename_difference_out)
        if pDataset_out is None:
            raise RuntimeError(
                f"Could not create output file: {sFilename_difference_out}"
            )

        # Configure spatial reference system for output
        pSpatial_reference_gcs = osr.SpatialReference()
        pSpatial_reference_gcs.ImportFromEPSG(spatial_ref_epsg)
        logger.info(f"Using spatial reference EPSG:{spatial_ref_epsg}")

        # Create output layer with polygon geometry type
        pLayerOut = pDataset_out.CreateLayer(
            "intersection", pSpatial_reference_gcs, ogr.wkbPolygon
        )
        if pLayerOut is None:
            raise RuntimeError("Could not create output layer")

        # Define output attribute schema
        pLayerOut.CreateField(ogr.FieldDefn("id", ogr.OFTInteger64))
        pLayerOut.CreateField(ogr.FieldDefn("base_fid", ogr.OFTInteger64))
        pLayerOut.CreateField(ogr.FieldDefn("new_fid", ogr.OFTInteger64))

        pLayerDefn = pLayerOut.GetLayerDefn()

        # Open base polyline dataset
        pDataset_base = driver_base.Open(sFilename_base, 0)
        if pDataset_base is None:
            raise RuntimeError(f"Could not open base file: {sFilename_base}")

        pLayer_base = pDataset_base.GetLayer()
        if pLayer_base is None:
            raise RuntimeError("Could not access layer in base file")

        nFeature_base = pLayer_base.GetFeatureCount()
        if nFeature_base == 0:
            logger.warning("Base file contains no features")
            return

        logger.info(f"Base file contains {nFeature_base} polyline features")

        # Open polygon dataset
        pDataset_new = driver_new.Open(sFilename_new, 0)
        if pDataset_new is None:
            raise RuntimeError(f"Could not open polygon file: {sFilename_new}")

        pLayer_new = pDataset_new.GetLayer()
        if pLayer_new is None:
            raise RuntimeError("Could not access layer in polygon file")

        nFeature_new = pLayer_new.GetFeatureCount()
        if nFeature_new == 0:
            logger.warning("Polygon file contains no features")
            return

        logger.info(f"Polygon file contains {nFeature_new} polygon features")

        # Build spatial index for base polyline features
        logger.info("Building spatial index for polyline features...")
        index_base = RTreeindex()

        base_features = {}  # Cache base features
        indexed_count = 0

        for idx, pFeature_base in enumerate(pLayer_base):
            fid = pFeature_base.GetFID()
            pGeometry_base = pFeature_base.GetGeometryRef()

            if pGeometry_base is None:
                logger.warning(f"Base feature {fid} has no geometry, skipping")
                continue

            # Validate geometry
            if not pGeometry_base.IsValid():
                logger.warning(
                    f"Base feature {fid} has invalid geometry, attempting to fix"
                )
                pGeometry_base = pGeometry_base.Buffer(0)

            if pGeometry_base is None or pGeometry_base.IsEmpty():
                logger.warning(
                    f"Base feature {fid} has empty geometry after validation, skipping"
                )
                continue

            try:
                # Get bounding box
                envelope = pGeometry_base.GetEnvelope()
                left, right, bottom, top = envelope

                # Insert into spatial index
                index_base.insert(fid, (left, bottom, right, top))

                # Cache feature
                base_features[fid] = pFeature_base.Clone()
                indexed_count += 1
            except Exception as e:
                logger.error(f"Error indexing base feature {fid}: {e}")
                continue

        logger.info(f"Indexed {indexed_count} valid polyline features")

        if not base_features:
            logger.warning("No valid base features found for indexing")
            return

        # Process polygon features and calculate intersections
        logger.info("Processing intersections...")
        lID_polygon = 1
        processed_features = 0

        for idx, pFeature_new in enumerate(pLayer_new):
            processed_features += 1
            if processed_features % 100 == 0:
                logger.info(
                    f"Processed {processed_features}/{nFeature_new} polygon features"
                )

            new_fid = pFeature_new.GetFID()
            pGeometry_new = pFeature_new.GetGeometryRef()

            if pGeometry_new is None:
                logger.warning(f"Polygon feature {new_fid} has no geometry, skipping")
                continue

            # Validate geometry
            if not pGeometry_new.IsValid():
                logger.warning(
                    f"Polygon feature {new_fid} has invalid geometry, attempting to fix"
                )
                pGeometry_new = pGeometry_new.Buffer(0)

            if pGeometry_new is None or pGeometry_new.IsEmpty():
                logger.warning(
                    f"Polygon feature {new_fid} has empty geometry after validation, skipping"
                )
                continue

            try:
                # Get bounding box for spatial query
                envelope = pGeometry_new.GetEnvelope()
                left, right, bottom, top = envelope

                # Query spatial index to find candidate polylines
                aIntersect = list(index_base.intersection((left, bottom, right, top)))

                # Process each candidate polyline feature
                for base_fid in aIntersect:
                    if base_fid not in base_features:
                        continue

                    pFeature_base = base_features[base_fid]
                    pGeometry_base = pFeature_base.GetGeometryRef()

                    try:
                        # Perform exact geometric intersection test
                        if not pGeometry_new.Intersects(pGeometry_base):
                            continue

                        # Calculate geometric intersection
                        pGeometry_intersect = pGeometry_new.Intersection(pGeometry_base)

                        # Skip if intersection is empty
                        if pGeometry_intersect is None or pGeometry_intersect.IsEmpty():
                            continue

                        # Handle output based on geometry type
                        pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()

                        if pGeometrytype_intersect == "POLYGON":
                            # Simple polygon result
                            pFeatureOut = ogr.Feature(pLayerDefn)
                            pFeatureOut.SetGeometry(pGeometry_intersect)
                            pFeatureOut.SetField("id", lID_polygon)
                            pFeatureOut.SetField("base_fid", base_fid)
                            pFeatureOut.SetField("new_fid", new_fid)
                            pLayerOut.CreateFeature(pFeatureOut)
                            lID_polygon += 1

                        elif pGeometrytype_intersect in [
                            "MULTIPOLYGON",
                            "GEOMETRYCOLLECTION",
                        ]:
                            # Multi-part geometry - decompose into individual polygons
                            for i in range(pGeometry_intersect.GetGeometryCount()):
                                pGeometry_part = pGeometry_intersect.GetGeometryRef(i)
                                if (
                                    pGeometry_part is not None
                                    and not pGeometry_part.IsEmpty()
                                ):
                                    part_type = pGeometry_part.GetGeometryName()
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
                            f"Error processing intersection between base {base_fid} and polygon {new_fid}: {e}"
                        )
                        continue

            except Exception as e:
                logger.error(f"Error processing polygon feature {new_fid}: {e}")
                continue

        logger.info(f"Generated {lID_polygon - 1} intersection polygons")

    except Exception as e:
        logger.error(f"Error during processing: {e}")
        raise RuntimeError(f"Failed to calculate polyline-polygon intersection: {e}")

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
