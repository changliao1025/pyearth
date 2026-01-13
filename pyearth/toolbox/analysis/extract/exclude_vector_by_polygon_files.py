"""
Exclude vector features by multiple polygon files with support for multiple vector formats.

This module provides functionality to exclude (remove) features from a vector file
that fall within or intersect a union of multiple polygon boundaries.
"""

import os
import logging
from typing import List
from osgeo import ogr

from pyearth.system.define_global_variables import *
from pyearth.toolbox.management.vector.reproject import reproject_vector
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_extension
from pyearth.gis.geometry.get_output_geometry_type import get_output_geometry_type

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def exclude_vector_by_polygon_files(
    sFilename_vector_in: str, aFilename_polygon_in: List[str], sFilename_vector_out: str
) -> None:
    """
    Exclude (remove) vector features that fall within a union of multiple polygon files.

    This function performs the inverse of clipping using multiple polygon files. It creates
    a union of all polygons from all input files, then keeps only the features or parts of
    features that are OUTSIDE this combined polygon boundary.

    Args:
        sFilename_vector_in: Path to input vector file (supports .geojson, .shp, .gpkg, etc.)
        aFilename_polygon_in: List of paths to exclusion polygon files
        sFilename_vector_out: Path to output vector file with excluded features

    Returns:
        None

    Raises:
        FileNotFoundError: If input files don't exist
        RuntimeError: If GDAL/OGR operations fail
        ValueError: If geometry type is not supported or polygon files are empty

    Note:
        - Supports multiple vector formats (GeoJSON, Shapefile, GeoPackage, KML, etc.)
        - Automatically handles projection differences (reprojects polygons if needed)
        - Merges multiple polygons within each file if necessary
        - Creates union of all polygons from all files
        - Uses geometric difference operation for partially overlapping features
        - Handles multi-part geometries (e.g., MultiLineString) properly

    Example:
        >>> exclude_vector_by_polygon_files(
        ...     'rivers.geojson',
        ...     ['lake1.shp', 'lake2.shp', 'reservoir.gpkg'],
        ...     'land_rivers.gpkg'  # Rivers excluding all water bodies
        ... )
    """
    logger.info("=" * 80)
    logger.info("Starting multi-polygon vector exclusion operation")
    logger.info(f"Input vector: {sFilename_vector_in}")
    logger.info(f"Number of polygon files: {len(aFilename_polygon_in)}")
    logger.info(f"Output vector: {sFilename_vector_out}")

    try:
        # Validate input files
        if not os.path.exists(sFilename_vector_in):
            raise FileNotFoundError(
                f"Input vector file not found: {sFilename_vector_in}"
            )

        for polygon_file in aFilename_polygon_in:
            if not os.path.exists(polygon_file):
                raise FileNotFoundError(f"Polygon file not found: {polygon_file}")

        # Remove existing output file
        if os.path.exists(sFilename_vector_out):
            os.remove(sFilename_vector_out)
            logger.info(f"Removed existing output file: {sFilename_vector_out}")

        sFolder_out = os.path.dirname(sFilename_vector_out)

        # Get output driver using multi-format support
        sExtension_out = os.path.splitext(sFilename_vector_out)[1]
        pDriver_out = get_vector_driver_from_extension(sExtension_out)
        if pDriver_out is None:
            raise RuntimeError(f"Could not get driver for extension: {sExtension_out}")

        logger.info(f"Using output driver: {pDriver_out.GetName()}")

        # Open input vector dataset
        pDataset_source = ogr.Open(sFilename_vector_in)
        if pDataset_source is None:
            raise RuntimeError(
                f"Could not open input vector file: {sFilename_vector_in}"
            )

        pLayer_source = pDataset_source.GetLayer()
        if pLayer_source is None:
            raise RuntimeError("Could not access layer in input vector file")

        feature_count = pLayer_source.GetFeatureCount()
        logger.info(f"Input vector contains {feature_count} features")

        # Get spatial reference from source
        pSpatial_reference_target = pLayer_source.GetSpatialRef()
        if pSpatial_reference_target is None:
            logger.warning("Input vector has no spatial reference")
            pProjection_target = None
        else:
            pProjection_target = pSpatial_reference_target.ExportToWkt()
            logger.info(
                f"Target projection: {pSpatial_reference_target.GetName() if pSpatial_reference_target else 'Unknown'}"
            )

        # Create output dataset
        pDataset_excluded = pDriver_out.CreateDataSource(sFilename_vector_out)
        if pDataset_excluded is None:
            raise RuntimeError(f"Could not create output file: {sFilename_vector_out}")

        # Determine output geometry type
        iGeomType = pLayer_source.GetGeomType()

        try:
            output_geom_type = get_output_geometry_type(iGeomType)
            geom_name = ogr.GeometryTypeToName(output_geom_type)
            logger.info(f"Output geometry type: {geom_name}")
        except ValueError as e:
            raise ValueError(str(e))

        # Create output layer
        pLayer_excluded = pDataset_excluded.CreateLayer(
            "layer", pSpatial_reference_target, geom_type=output_geom_type
        )
        if pLayer_excluded is None:
            raise RuntimeError("Could not create output layer")

        # Copy field definitions
        pFeatureDefn_source = pLayer_source.GetLayerDefn()
        for i in range(pFeatureDefn_source.GetFieldCount()):
            pFieldDefn_source = pFeatureDefn_source.GetFieldDefn(i)
            pLayer_excluded.CreateField(pFieldDefn_source)

        logger.info("Processing polygon files for union...")
        aFilename_polygon_processed = []

        # Process each polygon file
        for idx, sFilename_polygon_in in enumerate(aFilename_polygon_in, 1):
            logger.info(
                f"Processing polygon file {idx}/{len(aFilename_polygon_in)}: {sFilename_polygon_in}"
            )

            # Open polygon dataset
            pDataset_clip = ogr.Open(sFilename_polygon_in)
            if pDataset_clip is None:
                logger.warning(
                    f"Could not open polygon file, skipping: {sFilename_polygon_in}"
                )
                continue

            pLayer_clip = pDataset_clip.GetLayer()
            nPolygon = pLayer_clip.GetFeatureCount()
            logger.info(f"Polygon file contains {nPolygon} feature(s)")

            if nPolygon == 0:
                logger.warning("Polygon file contains no features, skipping")
                pDataset_clip = None
                continue

            sExtension_clip = os.path.splitext(sFilename_polygon_in)[1]

            # Handle multiple polygons by merging
            if nPolygon > 1:
                logger.info("Multiple polygons detected, merging into single polygon")
                sFilename_clip_merged = sFilename_polygon_in.replace(
                    sExtension_clip, "_merged" + sExtension_clip
                )

                # Close current dataset before merging
                pDataset_clip = None
                pLayer_clip = None

                # Merge features
                merge_features(sFilename_polygon_in, sFilename_clip_merged)
                sFilename_polygon_in = sFilename_clip_merged

                # Reopen merged file
                pDataset_clip = ogr.Open(sFilename_polygon_in)
                pLayer_clip = pDataset_clip.GetLayer()

            # Check projection compatibility
            pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
            if pSpatial_reference_clip is not None:
                pProjection_clip = pSpatial_reference_clip.ExportToWkt()
            else:
                pProjection_clip = None
                logger.warning(f"Polygon file {idx} has no spatial reference")

            sFilename_clip_final = sFilename_polygon_in
            if (
                pProjection_target is not None
                and pProjection_clip is not None
                and pProjection_target != pProjection_clip
            ):
                logger.info("Projection mismatch detected, reprojecting polygon file")

                # Close dataset before reprojection
                pDataset_clip = None
                pLayer_clip = None

                # Reproject to target projection
                sFolder = os.path.dirname(sFilename_polygon_in)
                sName = os.path.basename(sFilename_polygon_in)
                sName_no_ext = os.path.splitext(sName)[0]
                sFilename_clip_final = os.path.join(
                    sFolder, f"{sName_no_ext}_transformed{sExtension_clip}"
                )

                reproject_vector(
                    sFilename_polygon_in, sFilename_clip_final, pProjection_target
                )

                # Reopen reprojected file
                pDataset_clip = ogr.Open(sFilename_clip_final)
                pLayer_clip = pDataset_clip.GetLayer()

            # Add to processed list
            aFilename_polygon_processed.append(sFilename_clip_final)
            pDataset_clip = None
            pLayer_clip = None

        if len(aFilename_polygon_processed) == 0:
            raise ValueError("No valid polygon files could be processed!")

        logger.info(
            f"Successfully processed {len(aFilename_polygon_processed)} polygon file(s)"
        )

        # Create union of all polygons
        logger.info("Creating union of all polygon geometries...")
        pGeometry_union = None

        if len(aFilename_polygon_processed) == 1:
            # Single polygon file - just use its geometry
            pDataset_clip = ogr.Open(aFilename_polygon_processed[0])
            pLayer_clip = pDataset_clip.GetLayer()
            pFeature_clip = pLayer_clip.GetNextFeature()
            pGeometry_union = pFeature_clip.GetGeometryRef().Clone()
            pDataset_clip = None
            logger.info("Using single polygon geometry")
        else:
            # Multiple polygon files - create union
            logger.info(
                f"Creating union of {len(aFilename_polygon_processed)} polygon geometries..."
            )

            for idx, sFilename_polygon in enumerate(aFilename_polygon_processed, 1):
                pDataset_clip = ogr.Open(sFilename_polygon)
                pLayer_clip = pDataset_clip.GetLayer()

                for pFeature_clip in pLayer_clip:
                    pGeometry_current = pFeature_clip.GetGeometryRef()

                    if pGeometry_current is None:
                        continue

                    if pGeometry_union is None:
                        # First geometry - initialize union
                        pGeometry_union = pGeometry_current.Clone()
                    else:
                        # Subsequent geometries - union with existing
                        pGeometry_union = pGeometry_union.Union(pGeometry_current)

                pDataset_clip = None
                logger.info(f"Merged polygon {idx}/{len(aFilename_polygon_processed)}")

        if pGeometry_union is None:
            raise RuntimeError("Failed to create union of polygon geometries")

        # Get envelope for logging
        minx, maxx, miny, maxy = pGeometry_union.GetEnvelope()
        logger.info(
            f"Union polygon envelope: ({minx:.2f}, {miny:.2f}) to ({maxx:.2f}, {maxy:.2f})"
        )

        # Use the union as the exclusion polygon
        pPolygon_exclude = pGeometry_union

        logger.info("Processing features for exclusion...")
        excluded_count = 0
        kept_count = 0
        trimmed_count = 0
        trimmed_multipart_count = 0

        pLayer_source.ResetReading()

        # Process each feature
        for pFeature in pLayer_source:
            pGeometry_source = pFeature.GetGeometryRef()
            if pGeometry_source is None:
                logger.debug("Feature has no geometry, skipping")
                continue

            # Check if feature is entirely within the exclusion polygon
            if pPolygon_exclude.Contains(pGeometry_source):
                # Feature is entirely inside - EXCLUDE it
                excluded_count += 1
                continue

            # Check if feature intersects the exclusion polygon
            if pGeometry_source.Intersects(pPolygon_exclude):
                # Feature partially overlaps - need to trim it
                pGeometry_difference = pGeometry_source.Difference(pPolygon_exclude)

                if pGeometry_difference is None or pGeometry_difference.IsEmpty():
                    logger.debug(
                        "Difference resulted in empty geometry, excluding feature"
                    )
                    excluded_count += 1
                    continue

                iGeomType_difference = pGeometry_difference.GetGeometryType()

                # Check if geometry type matches expected output
                if get_output_geometry_type(iGeomType_difference) == output_geom_type:
                    pFeature_excluded = ogr.Feature(pLayer_excluded.GetLayerDefn())
                    pFeature_excluded.SetGeometry(pGeometry_difference)

                    # Copy attributes
                    for i in range(pFeature.GetFieldCount()):
                        pFeature_excluded.SetField(
                            pFeature.GetFieldDefnRef(i).GetNameRef(),
                            pFeature.GetField(i),
                        )

                    pLayer_excluded.CreateFeature(pFeature_excluded)
                    trimmed_count += 1

                elif (
                    iGeomType_difference == ogr.wkbMultiLineString
                    and output_geom_type == ogr.wkbLineString
                ):
                    # Handle MultiLineString -> LineString conversion
                    nPart = pGeometry_difference.GetGeometryCount()
                    logger.debug(
                        f"Splitting MultiLineString into {nPart} LineString parts"
                    )

                    for iPart in range(nPart):
                        pGeometry_part = pGeometry_difference.GetGeometryRef(iPart)
                        pFeature_excluded = ogr.Feature(pLayer_excluded.GetLayerDefn())
                        pFeature_excluded.SetGeometry(pGeometry_part)

                        # Copy attributes
                        for i in range(pFeature.GetFieldCount()):
                            pFeature_excluded.SetField(
                                pFeature.GetFieldDefnRef(i).GetNameRef(),
                                pFeature.GetField(i),
                            )

                        pLayer_excluded.CreateFeature(pFeature_excluded)

                    trimmed_multipart_count += 1
                else:
                    sGeomType_diff = ogr.GeometryTypeToName(iGeomType_difference)
                    logger.debug(
                        f"Difference geometry type mismatch: {sGeomType_diff}, excluding"
                    )
                    excluded_count += 1
            else:
                # Feature is entirely outside the exclusion polygon - KEEP it
                pFeature_excluded = ogr.Feature(pLayer_excluded.GetLayerDefn())
                pFeature_excluded.SetGeometry(pGeometry_source)

                # Copy attributes
                for i in range(pFeature.GetFieldCount()):
                    pFeature_excluded.SetField(
                        pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i)
                    )

                pLayer_excluded.CreateFeature(pFeature_excluded)
                kept_count += 1

        logger.info(f"Exclusion summary:")
        logger.info(f"  - Features kept (entirely outside): {kept_count}")
        logger.info(f"  - Features trimmed (partially inside): {trimmed_count}")
        logger.info(f"  - Features split from multipart: {trimmed_multipart_count}")
        logger.info(f"  - Features excluded (entirely inside): {excluded_count}")
        logger.info(f"  - Total output features: {kept_count + trimmed_count}")
        logger.info(
            f"Multi-polygon exclusion completed successfully: {sFilename_vector_out}"
        )

    except Exception as e:
        logger.error(f"Error during multi-polygon vector exclusion: {e}")
        raise RuntimeError(f"Failed to exclude vector with multiple polygons: {e}")

    finally:
        # Clean up resources
        if "pDataset_source" in locals():
            pDataset_source = None
        if "pDataset_clip" in locals():
            pDataset_clip = None
        if "pDataset_excluded" in locals():
            pDataset_excluded = None
        if "pSpatial_reference_clip" in locals():
            pSpatial_reference_clip = None
        if "pSpatial_reference_target" in locals():
            pSpatial_reference_target = None
        if "pGeometry_union" in locals():
            pGeometry_union = None

    return
