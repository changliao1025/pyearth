import os, sys
import numpy as np
import time
from typing import List, Tuple, Optional, Union, Dict, Any
from osgeo import gdal, osr, ogr
from rtree.index import Index as RTreeindex
from pyearth.system.filename import get_extension_from_path

from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_format_from_filename,
    print_supported_vector_formats,
    get_vector_driver_from_format,
    get_vector_driver_from_filename,
)


def create_spatial_index(
    geometries: List[Optional[Any]], verbose: bool = False
) -> Tuple[
    Union[Any, Dict[int, Tuple[float, float, float, float]]],
    Optional[Dict[int, Tuple[float, float, float, float]]],
    bool,
]:
    """
    Create a spatial index for geometries to speed up intersection queries.

    Uses RTree (rtree) for spatial indexing.

    Parameters
    ----------
    geometries : List[Optional[Any]]
        List of OGR geometry objects to index. None values are allowed and will be skipped.
    verbose : bool, optional
        Whether to print progress information. Default is False.

    Returns
    -------
    Tuple[Any, Dict[int, Tuple[float, float, float, float]], bool]
        A tuple containing:
        - spatial_index: RTree object
        - geometry_map: Dictionary mapping geometry indices to envelopes (minX, maxX, minY, maxY)

    Notes
    -----
    The spatial index enables efficient spatial queries for finding nearby geometries.
    RTree provides O(log n) query performance.
    """
    # Setup spatial indexing (rtree only)
    spatial_index = RTreeindex()
    geometry_map = {}

    for i, geom in enumerate(geometries):
        if geom is not None:
            envelope = geom.GetEnvelope()  # (minX, maxX, minY, maxY)
            # RTree expects (minX, minY, maxX, maxY)
            bbox = (envelope[0], envelope[2], envelope[1], envelope[3])
            spatial_index.insert(i, bbox)
            geometry_map[i] = envelope

    return spatial_index, geometry_map


def geometries_bbox_overlap(
    bbox1: Tuple[float, float, float, float],
    bbox2: Tuple[float, float, float, float],
    tolerance: float = 1e-10,
) -> bool:
    """
    Check if two bounding boxes overlap with optional tolerance.

    Parameters
    ----------
    bbox1 : Tuple[float, float, float, float]
        Bounding box of first geometry as (minX, maxX, minY, maxY).
    bbox2 : Tuple[float, float, float, float]
        Bounding box of second geometry as (minX, maxX, minY, maxY).
    tolerance : float, optional
        Small buffer for near-touching geometries. Default is 1e-10.

    Returns
    -------
    bool
        True if bounding boxes overlap (including touching), False otherwise.

    Notes
    -----
    Uses axis-aligned bounding box (AABB) overlap detection. Two boxes overlap if
    they overlap on both X and Y axes. Touching edges are considered overlapping
    when tolerance is applied.
    """
    return not (
        bbox1[1] + tolerance < bbox2[0]  # bbox1.maxX < bbox2.minX
        or bbox2[1] + tolerance < bbox1[0]  # bbox2.maxX < bbox1.minX
        or bbox1[3] + tolerance < bbox2[2]  # bbox1.maxY < bbox2.minY
        or bbox2[3] + tolerance < bbox1[2]
    )  # bbox2.maxY < bbox1.minY


def find_connectivity_groups_optimized(
    aGeometries: List[Optional[Any]], verbose: bool = False
) -> List[List[Any]]:
    """
    Find groups of connected geometries using optimized spatial indexing.

    Uses RTree for efficient spatial queries. Geometries are considered connected if they
    touch or intersect spatially.

    Parameters
    ----------
    aGeometries : List[Optional[Any]]
        List of OGR geometry objects to analyze for connectivity. None values are allowed and will be skipped.
    verbose : bool, optional
        Whether to print progress information. Default is False.

    Returns
    -------
    List[List[Any]]
        List of geometry groups, where each group is a list of connected geometries.
        Geometries in the same group are spatially connected (touching or intersecting).

    Notes
    -----
    Uses breadth-first search (BFS) algorithm to find connected components:
    1. Create spatial index for efficient candidate finding
    2. For each unprocessed geometry, start a new group
    3. Use BFS to find all connected geometries in the group
    4. Spatial relationships are checked using OGR's Touches() and Intersects() methods

    Performance:
    - RTree: O(n log n) for indexing + O(c) for connectivity where c is number of connections

    Examples
    --------
    >>> geometries = [point1, point2, polygon1, polygon2]  # OGR geometry objects
    >>> groups = find_connectivity_groups_optimized(geometries, verbose=True)
    >>> print(f"Found {len(groups)} connected groups")
    """
    if verbose:
        print(f"Finding connectivity groups for {len(aGeometries)} geometries...")

    start_time = time.time()

    # Create spatial index (rtree only)
    spatial_index, geometry_map = create_spatial_index(aGeometries, verbose)

    # Track processed geometries
    processed = set()
    groups = []

    for i, geom1 in enumerate(aGeometries):
        if i in processed or geom1 is None:
            continue

        # Start new group
        current_group = [geom1]
        current_indices = {i}
        processed.add(i)

        # Use queue for breadth-first search
        queue = [i]

        while queue:
            current_idx = queue.pop(0)
            current_geom = aGeometries[current_idx]

            if current_geom is None:
                continue

            # Get candidates using rtree spatial index
            current_envelope = current_geom.GetEnvelope()
            # RTree expects (minX, minY, maxX, maxY)
            query_bbox = (
                current_envelope[0],
                current_envelope[2],
                current_envelope[1],
                current_envelope[3],
            )

            # Add small tolerance for near-touching geometries
            tolerance = 1e-8
            expanded_bbox = (
                query_bbox[0] - tolerance,
                query_bbox[1] - tolerance,
                query_bbox[2] + tolerance,
                query_bbox[3] + tolerance,
            )

            candidate_indices = list(spatial_index.intersection(expanded_bbox))

            # Check spatial relationships for candidates
            for j in candidate_indices:
                if j in processed or j in current_indices:
                    continue

                geom2 = aGeometries[j]
                if geom2 is None:
                    continue

                # Expensive spatial relationship check only for candidates
                try:
                    if geom2.Touches(current_geom) or geom2.Intersects(current_geom):
                        current_group.append(geom2)
                        current_indices.add(j)
                        processed.add(j)
                        queue.append(j)
                except Exception as e:
                    if verbose:
                        print(
                            f"Warning: Failed to check geometry relationship: {str(e)}"
                        )
                    continue

        groups.append(current_group)

        if verbose and len(groups) % 100 == 0:
            elapsed = time.time() - start_time
            print(
                f"Processed {len(processed)}/{len(aGeometries)} geometries, found {len(groups)} groups (elapsed: {elapsed:.1f}s)"
            )

    if verbose:
        elapsed = time.time() - start_time
        print(
            f"Connectivity grouping completed in {elapsed:.1f}s using RTree, found {len(groups)} groups"
        )

    return groups


def cascaded_union(geometries: List[Any], verbose: bool = False) -> Optional[Any]:
    """
    Perform cascaded union for better performance with many geometries.

    Uses a divide-and-conquer approach for large geometry collections to avoid
    performance degradation that can occur with sequential union operations.

    Parameters
    ----------
    geometries : List[Any]
        List of OGR geometry objects to union together.
    verbose : bool, optional
        Whether to print progress information. Default is False.

    Returns
    -------
    Optional[Any]
        The unioned geometry result, or None if union fails or input is empty.

    Notes
    -----
    Union strategy depends on input size:
    - Empty list: Returns None
    - Single geometry: Returns the geometry unchanged
    - Two geometries: Direct union using geometry.Union()
    - 3-10 geometries: Sequential union starting with first geometry
    - 11+ geometries: Divide-and-conquer approach splitting into halves

    The cascaded approach helps avoid the O(n²) complexity that can occur
    with large sequential unions by balancing the union tree.

    Examples
    --------
    >>> polygons = [poly1, poly2, poly3, poly4, poly5]  # OGR polygon objects
    >>> merged = cascaded_union(polygons, verbose=True)
    >>> print(f"Union resulted in {merged.GetGeometryType()} geometry")
    """
    if not geometries:
        return None

    if len(geometries) == 1:
        return geometries[0]

    if len(geometries) == 2:
        return geometries[0].Union(geometries[1])

    # For larger groups, use divide-and-conquer approach
    if len(geometries) > 10:
        mid = len(geometries) // 2
        left_union = cascaded_union(geometries[:mid], verbose)
        right_union = cascaded_union(geometries[mid:], verbose)

        if left_union and right_union:
            return left_union.Union(right_union)
        elif left_union:
            return left_union
        elif right_union:
            return right_union
        else:
            return None

    # For smaller groups, use sequential union
    result = geometries[0].Clone()
    for i in range(1, len(geometries)):
        try:
            result = result.Union(geometries[i])
        except Exception as e:
            if verbose:
                print(f"Warning: Failed to union geometry {i}: {str(e)}")
            continue

    return result


def merge_features(
    sFilename_in: str,
    sFilename_out: str,
    verbose: bool = True,
    iFlag_force: bool = False,
) -> None:
    """
    Merge features in a vector file.
    By default, it merges based on connectivity (touching or intersecting geometries).
    If iFlag_force is True, all features are merged as one, ignoring connectivity.

    Features that touch or intersect will be merged into single features.
    Features that don't touch or intersect will remain as separate features.

    Parameters
    ----------
    sFilename_in : str
        Absolute path to input vector file. Supports all GDAL/OGR readable formats.
    sFilename_out : str
        Absolute path for output vector file. Format is auto-detected from file extension.
        Existing files will be overwritten.
    verbose : bool, optional
        Whether to print progress information and performance statistics. Default is True.

    Returns
    -------
    None
        Results are written to the output file. No return value.

    Raises
    ------
    This function handles errors internally and prints messages rather than raising exceptions.

    Notes
    -----
    **Algorithm Overview:**
    1. Load all geometries from input file
    2. Group geometries by spatial connectivity using BFS and spatial indexing
    3. Union geometries within each connected group
    4. Write merged features to output file

    **Connectivity Definition:**
    Geometries are considered connected if they touch or intersect spatially.
    Uses OGR's Touches() and Intersects() methods for precise spatial relationships.

    **Performance Characteristics:**
    - Time: O(n log n) with RTree, O(n²) worst case without spatial indexing
    - Memory: O(n) for geometry storage during processing
    - Scales well for thousands of features with proper spatial indexing

    **Output Schema:**
    - Inherits geometry type from first input feature
    - Copies all field definitions from input
    - Preserves spatial reference system
    - Each output feature represents one connected group

    **Supported Geometry Types:**
    - Point, LineString, Polygon (including Multi- variants)
    - Mixed geometry types in input are filtered to match the dominant type

    Examples
    --------
    Basic usage with verbose output:

    >>> merge_features(
    ...     sFilename_in='/data/parcels.shp',
    ...     sFilename_out='/output/merged_parcels.gpkg',
    ...     verbose=True
    ... )
    === Starting merge_features operation ===
    Input file: /data/parcels.shp
    Output file: /output/merged_parcels.gpkg
    Input format: ESRI Shapefile
    Output format: GPKG
    Number of features found in the input file: 1250
    Geometry type: POLYGON
    Starting geometry collection for 1250 features...
    Collected 1250 geometries, now grouping by connectivity...
    Finding connectivity groups for 1250 geometries...
    Connectivity grouping completed in 2.34s using RTree, found 45 groups
    Found 45 connected groups of geometries
    === Merge operation completed ===
    Total processing time: 3.12 seconds
    Successfully processed 1250 input features
    Created 45 output features based on connectivity
    Average processing rate: 400.6 features/second

    Silent processing:

    >>> merge_features('/input/data.geojson', '/output/merged.geojson', verbose=False)
    # No output, results written to file
    """
    start_time = time.time()

    if verbose:
        print("=== Starting merge_features operation ===")
        print(f"Input file: {sFilename_in}")
        print(f"Output file: {sFilename_out}")

    if not os.path.exists(sFilename_in):
        print(f"Error: Input file {sFilename_in} does not exist!")
        return

    input_format = get_vector_format_from_filename(sFilename_in)
    output_format = get_vector_format_from_filename(sFilename_out)
    if verbose:
        print(f"Input format: {input_format}")
        print(f"Output format: {output_format}")

    pDataset_in = ogr.Open(sFilename_in)
    if pDataset_in is None:
        print(f"Error: Could not open input file {sFilename_in}")
        return

    # Get the first layer in the file
    pLayer_in = pDataset_in.GetLayer(0)
    if pLayer_in is None:
        print("Error: No layer found in input file")
        pDataset_in = None
        return

    # Count the number of features (polygons)
    nFeature = pLayer_in.GetFeatureCount()
    # Get the spatial reference of the layer
    pSpatial_reference = pLayer_in.GetSpatialRef()

    if nFeature == 0:
        print("No features found in the input file")
        pDataset_in = None
        return
    else:
        if verbose:
            print(f"Number of features found in the input file: {nFeature}")

    # Create a new dataset using the output filename
    pDriver = get_vector_driver_from_format(output_format)
    if pDriver is None:
        print(f"Error: Driver {output_format} not available!")
        if verbose:
            print_supported_vector_formats()
        pDataset_in = None
        return

    if os.path.exists(sFilename_out):
        if verbose:
            print(f"Removing existing output file: {sFilename_out}")
        try:
            pDriver.DeleteDataSource(sFilename_out)
            print(f"Successfully removed existing output file: {sFilename_out}")
        except Exception as e:
            os.remove(sFilename_out)

    pDataset_out = pDriver.CreateDataSource(sFilename_out)
    if pDataset_out is None:
        print(f"Error: Could not create output file {sFilename_out}")
        pDataset_in = None
        return

    # obtain the geotype of first layer and
    iGeomType = pLayer_in.GetGeomType()
    # obtain the geotype of first geometry
    pLayer_in.ResetReading()
    # Obtain the first feature
    pFeature_first = pLayer_in.GetNextFeature()
    if pFeature_first is None:
        print("Error: No features found in input file")
        pDataset_in = None
        pDataset_out = None
        return

    pGeometry = pFeature_first.GetGeometryRef()
    if pGeometry is None:
        print("Error: No geometry found in first feature")
        pDataset_in = None
        pDataset_out = None
        return

    pGeometry.FlattenTo2D()
    # get the geometry type
    iGeomType = pGeometry.GetGeometryType()
    # get geometry type name
    sGeomType = ogr.GeometryTypeToName(iGeomType)

    if verbose:
        print(f"Geometry type: {sGeomType}")

    # Store all geometries for analysis
    aGeometries = [pGeometry.Clone()]

    # check whether it is a multi-geometry
    if (
        iGeomType == ogr.wkbMultiPoint
        or iGeomType == ogr.wkbMultiLineString
        or iGeomType == ogr.wkbMultiPolygon
    ):
        # get the number of geometries
        nGeom = pGeometry.GetGeometryCount()
        # get the first geometry
        pGeometry_single = pGeometry.GetGeometryRef(0)
        iGeomType = pGeometry_single.GetGeometryType()

    layer_name = "merged"

    if iGeomType == ogr.wkbPoint:
        pLayer_out = pDataset_out.CreateLayer(
            layer_name, pSpatial_reference, geom_type=ogr.wkbPoint
        )
    elif iGeomType == ogr.wkbLineString:
        pLayer_out = pDataset_out.CreateLayer(
            layer_name, pSpatial_reference, geom_type=ogr.wkbLineString
        )
    elif iGeomType == ogr.wkbPolygon:
        pLayer_out = pDataset_out.CreateLayer(
            layer_name, pSpatial_reference, geom_type=ogr.wkbPolygon
        )
    else:
        print(f"Error: Geometry type {sGeomType} not supported")
        pDataset_in = None
        pDataset_out = None
        return

    if pLayer_out is None:
        print("Error: Could not create output layer")
        pDataset_in = None
        pDataset_out = None
        return
    # Create a new layer in the output shapefile

    # Copy field definitions from input layer to output layer
    pLayerDefn_in = pLayer_in.GetLayerDefn()
    nFields = pLayerDefn_in.GetFieldCount()

    for i in range(nFields):
        pFieldDefn = pLayerDefn_in.GetFieldDefn(i)
        result = pLayer_out.CreateField(pFieldDefn)
        if result != ogr.OGRERR_NONE:
            print(
                f"Warning: Failed to create field {pFieldDefn.GetName()}, error code: {result}"
            )

    # Store first feature's attributes
    first_feature_attributes = {}
    for i in range(nFields):
        field_name = pLayerDefn_in.GetFieldDefn(i).GetName()
        field_value = pFeature_first.GetField(i)
        first_feature_attributes[field_name] = field_value

    if verbose:
        print(f"Starting geometry collection for {nFeature} features...")

    # Collect all geometries first
    collection_start = time.time()
    pLayer_in.ResetReading()  # Reset reading to start from the first feature again
    pFeature = pLayer_in.GetNextFeature()  # Get first feature again
    nProcessed = 0

    while pFeature:
        nProcessed += 1
        if verbose and nProcessed % 1000 == 0:
            print(f"Collected {nProcessed}/{nFeature} features")

        pGeometry = pFeature.GetGeometryRef()
        if pGeometry is not None:
            pGeometry.FlattenTo2D()
            # check geotype again
            iGeomType_new = pGeometry.GetGeometryType()
            if iGeomType_new == iGeomType:
                # Clone geometry before adding to collection
                geom_to_add = pGeometry.Clone()
                aGeometries.append(geom_to_add)
            else:
                # check whether the geometry type is a multi-geometry
                if (
                    iGeomType_new == ogr.wkbMultiPoint
                    or iGeomType_new == ogr.wkbMultiLineString
                    or iGeomType_new == ogr.wkbMultiPolygon
                ):
                    # get the number of geometries
                    nGeom = pGeometry.GetGeometryCount()
                    for j in range(nGeom):
                        pGeometry_single = pGeometry.GetGeometryRef(j)
                        # check again its geometry type
                        iGeomType_single = pGeometry_single.GetGeometryType()
                        if iGeomType_single == iGeomType:
                            # Clone geometry before adding to collection
                            geom_to_add = pGeometry_single.Clone()
                            aGeometries.append(geom_to_add)
                        else:
                            if verbose:
                                print(
                                    f"Warning: Geometry type {ogr.GeometryTypeToName(iGeomType_single)} not supported in feature {nProcessed}"
                                )
                else:
                    if verbose:
                        print(
                            f"Warning: Geometry type {ogr.GeometryTypeToName(iGeomType_new)} not supported in feature {nProcessed}"
                        )

        pFeature = pLayer_in.GetNextFeature()

    if verbose:
        print(
            f"Collected {len(aGeometries)} geometries, now grouping by connectivity..."
        )

    if iFlag_force:
        if verbose:
            print("iFlag_force is True, merging all features into a single group.")
        aGeometry_groups = [aGeometries]  # All geometries in one group
    else:
        # Group geometries by connectivity using optimized algorithm
        aGeometry_groups = find_connectivity_groups_optimized(aGeometries, verbose)

    if verbose:
        print(f"Found {len(aGeometry_groups)} connected groups of geometries")

    # Create output features for each group
    nOutput_features = 0
    for group_idx, geometry_group in enumerate(aGeometry_groups):
        if len(geometry_group) == 1:
            # Single geometry, no union needed
            pGeometry_result = geometry_group[0]
        else:
            # Multiple geometries, union them using optimized cascaded union
            if verbose and len(geometry_group) > 5:
                print(
                    f"Performing cascaded union for group {group_idx} with {len(geometry_group)} geometries"
                )

            pGeometry_result = cascaded_union(geometry_group, verbose)

            if pGeometry_result is None:
                if verbose:
                    print(
                        f"Warning: Failed to union geometries in group {group_idx}, skipping"
                    )
                continue

        # Skip if geometry is None after union
        if pGeometry_result is None:
            continue

        # Create a new feature for this group
        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        pFeature_out.SetGeometry(pGeometry_result)

        # Set all attributes from the first feature
        for field_name, field_value in first_feature_attributes.items():
            try:
                pFeature_out.SetField(field_name, field_value)
            except Exception as e:
                if verbose:
                    print(f"Warning: Could not set field {field_name}: {str(e)}")

        result = pLayer_out.CreateFeature(pFeature_out)
        if result != ogr.OGRERR_NONE:
            print(
                f"Error: Failed to create output feature {group_idx}, error code: {result}"
            )
        else:
            nOutput_features += 1

        pFeature_out = None

    # Cleanup
    pDataset_out = None
    pDataset_in = None

    # Performance summary
    total_time = time.time() - start_time

    if verbose:
        print("=== Merge operation completed ===")
        print(f"Total processing time: {total_time:.2f} seconds")
        print(f"Successfully processed {nProcessed} input features")
        print(f"Created {nOutput_features} output features based on connectivity")
        print(f"Found {len(aGeometry_groups)} connected groups")
        print(f"Average processing rate: {nProcessed/total_time:.1f} features/second")
        print(f"Output saved to: {sFilename_out}")

    return
