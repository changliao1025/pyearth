"""
Exclude vector features by polygon with support for multiple vector formats.

This module provides functionality to exclude (remove) features from a vector file
that fall within or intersect a polygon boundary - the inverse of clipping.
"""

import os
import logging
from typing import Optional
from osgeo import ogr

from pyearth.system.define_global_variables import *
from pyearth.toolbox.management.vector.reproject import reproject_vector
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_extension
from pyearth.gis.geometry.get_output_geometry_type import get_output_geometry_type

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def exclude_vector_by_polygon_file(
    sFilename_vector_in: str, sFilename_polygon_in: str, sFilename_vector_out: str
) -> None:
    """
    Exclude (remove) vector features that fall within a polygon boundary.

    This function performs the inverse of clipping - it keeps only the features or
    parts of features that are OUTSIDE the polygon boundary. Features entirely within
    the polygon are removed, and features that partially overlap are trimmed.

    Args:
        sFilename_vector_in: Path to input vector file (supports .geojson, .shp, .gpkg, etc.)
        sFilename_polygon_in: Path to exclusion polygon file (supports multiple formats)
        sFilename_vector_out: Path to output vector file with excluded features

    Returns:
        None

    Raises:
        FileNotFoundError: If input files don't exist
        RuntimeError: If GDAL/OGR operations fail
        ValueError: If geometry type is not supported or polygon file is empty

    Note:
        - Supports multiple vector formats (GeoJSON, Shapefile, GeoPackage, KML, etc.)
        - Automatically handles projection differences (reprojects polygon if needed)
        - Merges multiple polygons in the exclusion file if necessary
        - Uses geometric difference operation for partially overlapping features

    Example:
        >>> exclude_vector_by_polygon_file(
        ...     'roads.geojson',
        ...     'urban_area.shp',
        ...     'rural_roads.gpkg'  # Only roads outside urban area
        ... )
    """
    logger.info("=" * 80)
    logger.info("Starting vector exclusion operation")
    logger.info(f"Input vector: {sFilename_vector_in}")
    logger.info(f"Exclusion polygon: {sFilename_polygon_in}")
    logger.info(f"Output vector: {sFilename_vector_out}")

    try:
        # Validate input files
        if not os.path.exists(sFilename_vector_in):
            raise FileNotFoundError(
                f"Input vector file not found: {sFilename_vector_in}"
            )

        if not os.path.exists(sFilename_polygon_in):
            raise FileNotFoundError(f"Polygon file not found: {sFilename_polygon_in}")

        # Remove existing output file
        if os.path.exists(sFilename_vector_out):
            os.remove(sFilename_vector_out)
            logger.info(f"Removed existing output file: {sFilename_vector_out}")

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

        # Open exclusion polygon file
        pDataset_clip = ogr.Open(sFilename_polygon_in)
        if pDataset_clip is None:
            raise RuntimeError(f"Could not open polygon file: {sFilename_polygon_in}")

        pLayer_clip = pDataset_clip.GetLayer()
        if pLayer_clip is None:
            raise RuntimeError("Could not access layer in polygon file")

        pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
        if pSpatial_reference_clip is not None:
            pProjection_clip = pSpatial_reference_clip.ExportToWkt()
            logger.info(
                f"Polygon projection: {pSpatial_reference_clip.GetName() if pSpatial_reference_clip else 'Unknown'}"
            )
        else:
            pProjection_clip = None
            logger.warning("Polygon has no spatial reference")

        # Check polygon count and merge if needed
        nPolygon = pLayer_clip.GetFeatureCount()
        if nPolygon == 0:
            raise ValueError("The polygon file does not contain any features!")

        logger.info(f"Exclusion polygon contains {nPolygon} feature(s)")

        sExtension_clip = os.path.splitext(sFilename_polygon_in)[1]

        if nPolygon > 1:
            pDataset_clip = None
            pLayer_clip = None
            logger.info("Multiple polygons detected, merging into single feature...")
            sFilename_clip_new = sFilename_polygon_in.replace(
                sExtension_clip, "_merged" + sExtension_clip
            )
            merge_features(sFilename_polygon_in, sFilename_clip_new)
            sFilename_polygon_in = sFilename_clip_new

            # Reopen merged file
            pDataset_clip = ogr.Open(sFilename_polygon_in)
            pLayer_clip = pDataset_clip.GetLayer(0)
            pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
            if pSpatial_reference_clip is not None:
                pProjection_clip = pSpatial_reference_clip.ExportToWkt()

        # Handle projection mismatch
        if (
            pProjection_target is not None
            and pProjection_clip is not None
            and pProjection_target != pProjection_clip
        ):
            logger.info("Projection mismatch detected, reprojecting polygon...")
            pDataset_clip = None
            pLayer_clip = None

            sFolder = os.path.dirname(sFilename_polygon_in)
            sName = os.path.basename(sFilename_polygon_in)
            sName_no_extension = os.path.splitext(sName)[0]
            sFilename_clip_out = os.path.join(
                sFolder, sName_no_extension + "_transformed" + sExtension_clip
            )

            reproject_vector(
                sFilename_polygon_in, sFilename_clip_out, pProjection_target
            )
            pDataset_clip = ogr.Open(sFilename_clip_out)
            pLayer_clip = pDataset_clip.GetLayer(0)
            logger.info("Polygon reprojected")

        # Get exclusion polygon geometry
        pLayer_clip.ResetReading()
        pFeature_clip = pLayer_clip.GetNextFeature()
        if pFeature_clip is None:
            raise RuntimeError("Could not read polygon feature")

        pPolygon_clip = pFeature_clip.GetGeometryRef()
        if pPolygon_clip is None:
            raise RuntimeError("Polygon feature has no geometry")

        # Get output driver using multi-format support
        sExtension_out = os.path.splitext(sFilename_vector_out)[1]
        pDriver_out = get_vector_driver_from_extension(sExtension_out)
        if pDriver_out is None:
            raise RuntimeError(f"Could not get driver for extension: {sExtension_out}")

        logger.info(f"Using driver: {pDriver_out.GetName()}")

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

        logger.info("Processing features for exclusion...")
        excluded_count = 0
        kept_count = 0
        trimmed_count = 0

        pLayer_source.ResetReading()

        # Process each feature
        for pFeature in pLayer_source:
            pGeometry_source = pFeature.GetGeometryRef()
            if pGeometry_source is None:
                logger.debug("Feature has no geometry, skipping")
                continue

            # Check if feature is entirely within the exclusion polygon
            if pPolygon_clip.Contains(pGeometry_source):
                # Feature is entirely inside - EXCLUDE it (don't add to output)
                excluded_count += 1
                continue

            # Check if feature intersects the exclusion polygon
            if pGeometry_source.Intersects(pPolygon_clip):
                # Feature partially overlaps - need to trim it
                # Use Difference operation to keep only parts outside the polygon
                pGeometry_difference = pGeometry_source.Difference(pPolygon_clip)

                if pGeometry_difference is None or pGeometry_difference.IsEmpty():
                    logger.debug(
                        "Difference resulted in empty geometry, excluding feature"
                    )
                    excluded_count += 1
                    continue

                iGeomType_difference = pGeometry_difference.GetGeometryType()

                # Only keep if geometry type matches expected output
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
        logger.info(f"  - Features excluded (entirely inside): {excluded_count}")
        logger.info(f"  - Total output features: {kept_count + trimmed_count}")
        logger.info(f"Exclusion completed successfully: {sFilename_vector_out}")

    except Exception as e:
        logger.error(f"Error during vector exclusion: {e}")
        raise RuntimeError(f"Failed to exclude vector: {e}")

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

    return
