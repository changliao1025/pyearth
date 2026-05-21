"""
Union overlapping polygons in a polygon file.

This module provides functionality to dissolve overlapping polygons while keeping
non-overlapping polygons as individual features.
"""

import os
import logging
from osgeo import ogr
from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_format_from_filename,
    get_vector_driver_from_filename,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def union_polygon_file(sFilename_in, sFilename_out, bFill_holes=False):
    """
    Union overlapping polygons while keeping non-overlapping polygons separate.

    This function processes a polygon file and:
    - Groups overlapping polygons together
    - Unions/dissolves each group of overlapping polygons
    - Keeps non-overlapping polygons as individual features

    Args:
        sFilename_in: Path to input polygon file
        sFilename_out: Path to output polygon file
        bFill_holes: If True, fills holes/voids between overlapping polygons.
                     If False (default), preserves holes using Shapely's unary_union.

    Returns:
        None

    Raises:
        RuntimeError: If GDAL/OGR operations fail
        FileNotFoundError: If input file doesn't exist

    Example:
        >>> # Preserve holes between polygons
        >>> union_polygon_file('input.geojson', 'output.geojson', bFill_holes=False)
        >>> # Fill all holes/voids
        >>> union_polygon_file('input.geojson', 'output.geojson', bFill_holes=True)
    """
    logger.info("=" * 80)
    logger.info("Starting polygon union operation")
    logger.info(f"Input file: {sFilename_in}")
    logger.info(f"Output file: {sFilename_out}")
    logger.info(f"Fill holes: {bFill_holes}")

    # Validate input file
    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"Input file not found: {sFilename_in}")

    # Remove existing output file
    if os.path.exists(sFilename_out):
        os.remove(sFilename_out)
        logger.info(f"Removed existing output file: {sFilename_out}")

    # Open the input file
    pDataset = ogr.Open(sFilename_in)
    if pDataset is None:
        raise RuntimeError(f"Could not open input file: {sFilename_in}")

    pLayer = pDataset.GetLayer()
    if pLayer is None:
        raise RuntimeError("Could not access layer in input file")

    srs = pLayer.GetSpatialRef()
    feature_count = pLayer.GetFeatureCount()
    logger.info(f"Input file contains {feature_count} polygon(s)")

    # Read all geometries into a list
    geometries = []
    pLayer.ResetReading()
    for feature in pLayer:
        geometry = feature.GetGeometryRef()
        if geometry is not None:
            geometries.append(geometry.Clone())

    logger.info(f"Loaded {len(geometries)} valid geometries")

    # Group overlapping polygons using a union-find approach
    # Each polygon starts in its own group
    n = len(geometries)
    groups = list(range(n))  # groups[i] = group id for polygon i

    def find_group(i):
        """Find the root group for polygon i."""
        if groups[i] != i:
            groups[i] = find_group(groups[i])  # Path compression
        return groups[i]

    def union_groups(i, j):
        """Union the groups containing polygons i and j."""
        root_i = find_group(i)
        root_j = find_group(j)
        if root_i != root_j:
            groups[root_j] = root_i

    # Check all pairs of polygons for overlap
    logger.info("Checking for overlapping polygons...")
    overlap_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if geometries[i].Intersects(geometries[j]):
                union_groups(i, j)
                overlap_count += 1

    logger.info(f"Found {overlap_count} overlapping polygon pair(s)")

    # Group polygons by their root group
    polygon_groups = {}
    for i in range(n):
        root = find_group(i)
        if root not in polygon_groups:
            polygon_groups[root] = []
        polygon_groups[root].append(i)

    logger.info(f"Organized into {len(polygon_groups)} group(s)")

    # Get driver from output file name
    driver = get_vector_driver_from_filename(sFilename_out)
    if driver is None:
        raise RuntimeError(f"Could not get driver for output file: {sFilename_out}")

    # Create the output dataset
    pDataset_out = driver.CreateDataSource(sFilename_out)
    if pDataset_out is None:
        raise RuntimeError(f"Could not create output file: {sFilename_out}")

    # Create the output layer with the same spatial reference as the input
    pLayer_out = pDataset_out.CreateLayer("union", srs, ogr.wkbPolygon)
    if pLayer_out is None:
        raise RuntimeError(f"Could not create layer in output file: {sFilename_out}")

    # Get feature definition for output
    feature_defn = pLayer_out.GetLayerDefn()

    # Process each group
    logger.info("Processing polygon groups...")
    single_polygon_count = 0
    union_polygon_count = 0

    for group_id, polygon_indices in polygon_groups.items():
        if len(polygon_indices) == 1:
            # Single polygon - no overlap, keep as is
            idx = polygon_indices[0]
            feature_out = ogr.Feature(feature_defn)
            feature_out.SetGeometry(geometries[idx])
            pLayer_out.CreateFeature(feature_out)
            feature_out = None
            single_polygon_count += 1
        else:
            # Multiple overlapping polygons - union them
            if bFill_holes:
                # Use GDAL Union - fills all holes/voids
                union_geometry = geometries[polygon_indices[0]].Clone()
                for idx in polygon_indices[1:]:
                    union_geometry = union_geometry.Union(geometries[idx])

                logger.info(f"  Unioned {len(polygon_indices)} overlapping polygons (holes filled)")
            else:
                # Use Shapely unary_union - preserves holes/voids
                try:
                    from shapely.geometry import shape, mapping
                    from shapely.ops import unary_union

                    # Convert OGR geometries to Shapely
                    shapely_polygons = []
                    for idx in polygon_indices:
                        geom_json = geometries[idx].ExportToJson()
                        shapely_geom = shape(eval(geom_json))
                        shapely_polygons.append(shapely_geom)

                    # Perform unary_union (preserves holes)
                    union_result = unary_union(shapely_polygons)

                    # Convert back to OGR
                    geom_json = mapping(union_result)
                    union_geometry = ogr.CreateGeometryFromJson(str(geom_json))

                    logger.info(f"  Unioned {len(polygon_indices)} overlapping polygons (holes preserved)")
                except ImportError:
                    logger.warning("Shapely not available, falling back to GDAL Union (will fill holes)")
                    union_geometry = geometries[polygon_indices[0]].Clone()
                    for idx in polygon_indices[1:]:
                        union_geometry = union_geometry.Union(geometries[idx])
                    logger.info(f"  Unioned {len(polygon_indices)} overlapping polygons (holes filled)")

            # Handle potential MultiPolygon result
            if union_geometry is not None and not union_geometry.IsEmpty():
                geom_type = union_geometry.GetGeometryType()

                if geom_type == ogr.wkbMultiPolygon:
                    # If result is MultiPolygon, save each polygon separately
                    for i in range(union_geometry.GetGeometryCount()):
                        sub_geom = union_geometry.GetGeometryRef(i)
                        feature_out = ogr.Feature(feature_defn)
                        feature_out.SetGeometry(sub_geom.Clone())
                        pLayer_out.CreateFeature(feature_out)
                        feature_out = None
                        union_polygon_count += 1
                else:
                    # Single polygon result
                    feature_out = ogr.Feature(feature_defn)
                    feature_out.SetGeometry(union_geometry)
                    pLayer_out.CreateFeature(feature_out)
                    feature_out = None
                    union_polygon_count += 1

    # Clean up
    pDataset = None
    pDataset_out = None

    logger.info("Union operation summary:")
    logger.info(f"  - Non-overlapping polygons: {single_polygon_count}")
    logger.info(f"  - Unioned polygon groups: {union_polygon_count}")
    logger.info(f"  - Total output features: {single_polygon_count + union_polygon_count}")
    logger.info(f"Union completed successfully: {sFilename_out}")

    return
