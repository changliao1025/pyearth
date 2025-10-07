import os, sys
import numpy as np
import time
from osgeo import gdal, osr, ogr, gdalconst
import importlib.util
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from tinyr import RTree
    iFlag_use_rtree = 1
else:
    iFlag_use_rtree = 0
from pyearth.gis.gdal.read.vector.get_supported_formats import get_supported_formats, print_supported_formats

def create_spatial_index(geometries, verbose=False):
    """
    Create a spatial index for geometries to speed up intersection queries
    Uses RTree if available, otherwise falls back to simple dictionary

    :param geometries: list of OGR geometry objects
    :param verbose: whether to print progress information
    :return: spatial index (RTree or dictionary) and geometry mapping
    """
    if verbose:
        index_type = "RTree" if iFlag_use_rtree else "dictionary"
        print(f"Creating spatial index using {index_type}...")

    if iFlag_use_rtree:
        # Use RTree for efficient spatial indexing
        rtree = RTree()
        geometry_map = {}

        for i, geom in enumerate(geometries):
            if geom is not None:
                envelope = geom.GetEnvelope()  # (minX, maxX, minY, maxY)
                # RTree expects (minX, minY, maxX, maxY)
                bbox = (envelope[0], envelope[2], envelope[1], envelope[3])
                rtree.insert(i, bbox)
                geometry_map[i] = envelope

        return rtree, geometry_map
    else:
        # Fallback to simple dictionary-based index
        spatial_index = {}
        for i, geom in enumerate(geometries):
            if geom is not None:
                envelope = geom.GetEnvelope()  # (minX, maxX, minY, maxY)
                spatial_index[i] = envelope

        return spatial_index, None

def geometries_bbox_overlap(bbox1, bbox2, tolerance=1e-10):
    """
    Check if two bounding boxes overlap with optional tolerance

    :param bbox1: (minX, maxX, minY, maxY) for first geometry
    :param bbox2: (minX, maxX, minY, maxY) for second geometry
    :param tolerance: small buffer for near-touching geometries
    :return: True if bounding boxes overlap
    """
    return not (bbox1[1] + tolerance < bbox2[0] or  # bbox1.maxX < bbox2.minX
                bbox2[1] + tolerance < bbox1[0] or  # bbox2.maxX < bbox1.minX
                bbox1[3] + tolerance < bbox2[2] or  # bbox1.maxY < bbox2.minY
                bbox2[3] + tolerance < bbox1[2])    # bbox2.maxY < bbox1.minY

def find_connectivity_groups_optimized(aGeometries, verbose=False):
    """
    Find groups of connected geometries using optimized spatial indexing
    Uses RTree for efficient spatial queries when available

    :param aGeometries: list of OGR geometry objects
    :param verbose: whether to print progress information
    :return: list of geometry groups (each group is a list of connected geometries)
    """
    if verbose:
        print(f'Finding connectivity groups for {len(aGeometries)} geometries...')

    start_time = time.time()

    # Create spatial index (RTree or dictionary)
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

            # Get candidates using spatial index
            if iFlag_use_rtree:
                # Use RTree for efficient spatial query
                current_envelope = current_geom.GetEnvelope()
                # RTree expects (minX, minY, maxX, maxY)
                query_bbox = (current_envelope[0], current_envelope[2], current_envelope[1], current_envelope[3])

                # Add small tolerance for near-touching geometries
                tolerance = 1e-8
                expanded_bbox = (
                    query_bbox[0] - tolerance,
                    query_bbox[1] - tolerance,
                    query_bbox[2] + tolerance,
                    query_bbox[3] + tolerance
                )

                candidate_indices = list(spatial_index.search(expanded_bbox))
            else:
                # Fallback to checking all geometries with bbox overlap
                current_bbox = spatial_index.get(current_idx)
                if current_bbox is None:
                    continue

                candidate_indices = []
                for j, geom2 in enumerate(aGeometries):
                    if j in processed or j in current_indices or geom2 is None:
                        continue

                    bbox2 = spatial_index.get(j)
                    if bbox2 is None:
                        continue

                    # Quick bounding box check first
                    if geometries_bbox_overlap(current_bbox, bbox2, tolerance=1e-8):
                        candidate_indices.append(j)

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
                        print(f'Warning: Failed to check geometry relationship: {str(e)}')
                    continue

        groups.append(current_group)

        if verbose and len(groups) % 100 == 0:
            elapsed = time.time() - start_time
            print(f'Processed {len(processed)}/{len(aGeometries)} geometries, found {len(groups)} groups (elapsed: {elapsed:.1f}s)')

    if verbose:
        elapsed = time.time() - start_time
        index_type = "RTree" if iFlag_use_rtree else "dictionary"
        print(f'Connectivity grouping completed in {elapsed:.1f}s using {index_type}, found {len(groups)} groups')

    return groups

def cascaded_union(geometries, verbose=False):
    """
    Perform cascaded union for better performance with many geometries

    :param geometries: list of OGR geometry objects to union
    :param verbose: whether to print progress information
    :return: resulting union geometry
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
                print(f'Warning: Failed to union geometry {i}: {str(e)}')
            continue

    return result

def merge_features(sFilename_in, sFilename_out, sFormat=None, verbose=True):
    """
    Merge features in a vector file based on connectivity (touching or intersecting geometries)

    Features that touch or intersect will be merged into single features.
    Features that don't touch or intersect will remain as separate features.

    :param sFilename_in: input vector file
    :param sFilename_out: output vector file
    :param sFormat: output format (auto-detected from extension if None)
    :param verbose: whether to print progress information
    :return: None
    """

    start_time = time.time()

    if verbose:
        print(f"=== Starting merge_features operation ===")
        print(f"Input file: {sFilename_in}")
        print(f"Output file: {sFilename_out}")

    if not os.path.exists(sFilename_in):
        print(f'Error: Input file {sFilename_in} does not exist!')
        return

    # Auto-detect format from file extension if not specified
    if sFormat is None:
        output_ext = os.path.splitext(sFilename_out)[1].lower()
        driver_map = get_supported_formats()
        sFormat = driver_map.get(output_ext, 'GeoJSON')

    if verbose:
        input_ext = os.path.splitext(sFilename_in)[1].lower()
        print(f'Input format: {input_ext}')
        print(f'Output format: {sFormat}')

    pDataset_in = ogr.Open(sFilename_in)
    if pDataset_in is None:
        print(f'Error: Could not open input file {sFilename_in}')
        return

    # Get the first layer in the file
    pLayer_in = pDataset_in.GetLayer(0)
    if pLayer_in is None:
        print('Error: No layer found in input file')
        pDataset_in = None
        return

    # Count the number of features (polygons)
    nFeature = pLayer_in.GetFeatureCount()
    # Get the spatial reference of the layer
    pSpatial_reference = pLayer_in.GetSpatialRef()

    if nFeature == 0:
        print('No features found in the input file')
        pDataset_in = None
        return
    else:
        if verbose:
            print(f'Number of features found in the input file: {nFeature}')

    # Create a new dataset using the output filename
    pDriver = ogr.GetDriverByName(sFormat)
    if pDriver is None:
        print(f'Error: Driver {sFormat} not available!')
        if verbose:
            print_supported_formats()
        pDataset_in = None
        return

    if os.path.exists(sFilename_out):
        if verbose:
            print(f'Removing existing output file: {sFilename_out}')
        try:
            pDriver.DeleteDataSource(sFilename_out)
            print(f'Successfully removed existing output file: {sFilename_out}')
        except Exception as e:
            os.remove(sFilename_out)


    pDataset_out = pDriver.CreateDataSource(sFilename_out)
    if pDataset_out is None:
        print(f'Error: Could not create output file {sFilename_out}')
        pDataset_in = None
        return

    #obtain the geotype of first layer and
    iGeomType = pLayer_in.GetGeomType()
    #obtain the geotype of first geometry
    pLayer_in.ResetReading()
    # Obtain the first feature
    pFeature_first = pLayer_in.GetNextFeature()
    if pFeature_first is None:
        print('Error: No features found in input file')
        pDataset_in = None
        pDataset_out = None
        return

    pGeometry = pFeature_first.GetGeometryRef()
    if pGeometry is None:
        print('Error: No geometry found in first feature')
        pDataset_in = None
        pDataset_out = None
        return

    pGeometry.FlattenTo2D()
    #get the geometry type
    iGeomType = pGeometry.GetGeometryType()
    #get geometry type name
    sGeomType = ogr.GeometryTypeToName(iGeomType)

    if verbose:
        print(f'Geometry type: {sGeomType}')

    # Store all geometries for analysis
    aGeometries = [pGeometry.Clone()]

    #check whether it is a multi-geometry
    if iGeomType == ogr.wkbMultiPoint or iGeomType == ogr.wkbMultiLineString or iGeomType == ogr.wkbMultiPolygon:
        #get the number of geometries
        nGeom = pGeometry.GetGeometryCount()
        #get the first geometry
        pGeometry_single = pGeometry.GetGeometryRef(0)
        iGeomType = pGeometry_single.GetGeometryType()

    # Create output layer with appropriate geometry type and naming
    if sFormat == 'ESRI Shapefile':
        # Shapefile requires simpler layer name
        layer_name = os.path.splitext(os.path.basename(sFilename_out))[0]
    else:
        layer_name = 'merged'

    if iGeomType == ogr.wkbPoint:
        pLayer_out = pDataset_out.CreateLayer(layer_name, pSpatial_reference, geom_type=ogr.wkbPoint)
    elif iGeomType == ogr.wkbLineString:
        pLayer_out = pDataset_out.CreateLayer(layer_name, pSpatial_reference, geom_type=ogr.wkbLineString)
    elif iGeomType == ogr.wkbPolygon:
        pLayer_out = pDataset_out.CreateLayer(layer_name, pSpatial_reference, geom_type=ogr.wkbPolygon)
    else:
        print(f'Error: Geometry type {sGeomType} not supported')
        pDataset_in = None
        pDataset_out = None
        return

    if pLayer_out is None:
        print('Error: Could not create output layer')
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
            print(f"Warning: Failed to create field {pFieldDefn.GetName()}, error code: {result}")

    # Store first feature's attributes
    first_feature_attributes = {}
    for i in range(nFields):
        field_name = pLayerDefn_in.GetFieldDefn(i).GetName()
        field_value = pFeature_first.GetField(i)
        first_feature_attributes[field_name] = field_value

    if verbose:
        print(f'Starting geometry collection for {nFeature} features...')

    # Collect all geometries first
    collection_start = time.time()
    pLayer_in.ResetReading()  # Reset reading to start from the first feature again
    pFeature = pLayer_in.GetNextFeature()  # Get first feature again
    nProcessed = 0

    while pFeature:
        nProcessed += 1
        if verbose and nProcessed % 1000 == 0:
            print(f'Collected {nProcessed}/{nFeature} features')

        pGeometry = pFeature.GetGeometryRef()
        if pGeometry is not None:
            pGeometry.FlattenTo2D()
            #check geotype again
            iGeomType_new = pGeometry.GetGeometryType()
            if iGeomType_new == iGeomType:
                # Clone geometry before adding to collection
                geom_to_add = pGeometry.Clone()
                aGeometries.append(geom_to_add)
            else:
                #check whether the geometry type is a multi-geometry
                if iGeomType_new == ogr.wkbMultiPoint or iGeomType_new == ogr.wkbMultiLineString or iGeomType_new == ogr.wkbMultiPolygon:
                    #get the number of geometries
                    nGeom = pGeometry.GetGeometryCount()
                    for j in range(nGeom):
                        pGeometry_single = pGeometry.GetGeometryRef(j)
                        #check again its geometry type
                        iGeomType_single = pGeometry_single.GetGeometryType()
                        if iGeomType_single == iGeomType:
                            # Clone geometry before adding to collection
                            geom_to_add = pGeometry_single.Clone()
                            aGeometries.append(geom_to_add)
                        else:
                            if verbose:
                                print(f'Warning: Geometry type {ogr.GeometryTypeToName(iGeomType_single)} not supported in feature {nProcessed}')
                else:
                    if verbose:
                        print(f'Warning: Geometry type {ogr.GeometryTypeToName(iGeomType_new)} not supported in feature {nProcessed}')

        pFeature = pLayer_in.GetNextFeature()

    if verbose:
        print(f'Collected {len(aGeometries)} geometries, now grouping by connectivity...')

    # Group geometries by connectivity using optimized algorithm
    aGeometry_groups = find_connectivity_groups_optimized(aGeometries, verbose)

    if verbose:
        print(f'Found {len(aGeometry_groups)} connected groups of geometries')

    # Create output features for each group
    nOutput_features = 0
    for group_idx, geometry_group in enumerate(aGeometry_groups):
        if len(geometry_group) == 1:
            # Single geometry, no union needed
            pGeometry_result = geometry_group[0]
        else:
            # Multiple geometries, union them using optimized cascaded union
            if verbose and len(geometry_group) > 5:
                print(f'Performing cascaded union for group {group_idx} with {len(geometry_group)} geometries')

            pGeometry_result = cascaded_union(geometry_group, verbose)

            if pGeometry_result is None:
                if verbose:
                    print(f'Warning: Failed to union geometries in group {group_idx}, skipping')
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
                    print(f'Warning: Could not set field {field_name}: {str(e)}')

        result = pLayer_out.CreateFeature(pFeature_out)
        if result != ogr.OGRERR_NONE:
            print(f'Error: Failed to create output feature {group_idx}, error code: {result}')
        else:
            nOutput_features += 1

        pFeature_out = None

    # Cleanup
    pDataset_out = None
    pDataset_in = None

    # Performance summary
    total_time = time.time() - start_time

    if verbose:
        print(f'=== Merge operation completed ===')
        print(f'Total processing time: {total_time:.2f} seconds')
        print(f'Successfully processed {nProcessed} input features')
        print(f'Created {nOutput_features} output features based on connectivity')
        print(f'Found {len(aGeometry_groups)} connected groups')
        print(f'Average processing rate: {nProcessed/total_time:.1f} features/second')
        print(f'Output saved to: {sFilename_out}')

    return

def main():
    """
    Main function for command line usage
    """
    import argparse

    parser = argparse.ArgumentParser(description='Merge connected features in a vector file (touching or intersecting geometries)')
    parser.add_argument('input', help='Input vector file')
    parser.add_argument('output', help='Output vector file')
    parser.add_argument('--format', help='Output format (auto-detected from extension if not specified)')
    parser.add_argument('--quiet', action='store_true', help='Suppress progress output')
    parser.add_argument('--formats', action='store_true', help='Show supported formats and exit')

    args = parser.parse_args()

    if args.formats:
        print_supported_formats()
        return

    merge_features(
        args.input,
        args.output,
        sFormat=args.format,
        verbose=not args.quiet
    )

if __name__ == '__main__':
    main()