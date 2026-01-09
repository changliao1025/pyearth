"""
Clip vector files by polygon with support for multiple vector formats.
"""

import os
import logging
from typing import Optional, List, Dict
from osgeo import ogr, osr

from pyearth.system.define_global_variables import *
from pyearth.toolbox.management.vector.reproject import reproject_vector
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_extension
from pyearth.gis.geometry.get_output_geometry_type import get_output_geometry_type

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def clip_vector_by_polygon_file(
    sFilename_vector_in: str,
    sFilename_polygon_in: str,
    sFilename_vector_out: str,
    iFlag_endorheic: int = 0,
) -> None:
    """
    Clip a vector file by a polygon file with support for multiple vector formats.

    This function clips an input vector file (points, lines, or polygons) to the extent
    of a polygon boundary. Supports automatic reprojection and merging of multiple polygons.
    Works with any OGR-supported vector format.

    Args:
        sFilename_vector_in: Path to input vector file (supports .geojson, .shp, .gpkg, etc.)
        sFilename_polygon_in: Path to clipping polygon file (supports multiple formats)
        sFilename_vector_out: Path to output clipped vector file
        iFlag_endorheic: If 1, skip merging for endorheic basins (default: 0)

    Returns:
        None

    Raises:
        FileNotFoundError: If input files don't exist
        RuntimeError: If GDAL/OGR operations fail
        ValueError: If geometry type is not supported or polygon file is empty

    Example:
        >>> clip_vector_by_polygon_file(
        ...     'rivers.geojson',
        ...     'study_area.shp',
        ...     'clipped_rivers.gpkg'
        ... )
    """

    logger.info("Starting vector clipping operation")
    logger.info(f"Input vector: {sFilename_vector_in}")
    logger.info(f"Clip polygon: {sFilename_polygon_in}")
    logger.info(f"Output vector: {sFilename_vector_out}")

    # Validate input files
    if not os.path.exists(sFilename_vector_in):
        raise FileNotFoundError(f"Input vector file not found: {sFilename_vector_in}")

    if not os.path.exists(sFilename_polygon_in):
        raise FileNotFoundError(f"Polygon file not found: {sFilename_polygon_in}")

    # Remove existing output file
    if os.path.exists(sFilename_vector_out):
        os.remove(sFilename_vector_out)
        logger.info(f"Removed existing output file: {sFilename_vector_out}")

    try:
        # Get drivers using multi-format support
        try:
            pDriver_vector_in = get_vector_driver_from_extension(sFilename_vector_in)
            pDriver_vector_out = get_vector_driver_from_extension(sFilename_vector_out)
            logger.info("Vector drivers initialized")
        except ValueError as e:
            raise RuntimeError(f"Unsupported file format: {e}")

        # Open input vector file
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
        else:
            pProjection_target = pSpatial_reference_target.ExportToWkt()
            logger.info(
                f"Target projection: {pSpatial_reference_target.GetName() if pSpatial_reference_target else 'Unknown'}"
            )

        # Open clip polygon file
        pDataset_clip = ogr.Open(sFilename_polygon_in)
        if pDataset_clip is None:
            raise RuntimeError(
                f"Could not open clip polygon file: {sFilename_polygon_in}"
            )

        pLayer_clip = pDataset_clip.GetLayer()
        if pLayer_clip is None:
            raise RuntimeError("Could not access layer in polygon file")

        pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
        if pSpatial_reference_clip is not None:
            pProjection_clip = pSpatial_reference_clip.ExportToWkt()
            logger.info(
                f"Clip polygon projection: {pSpatial_reference_clip.GetName() if pSpatial_reference_clip else 'Unknown'}"
            )
        else:
            pProjection_clip = None
            logger.warning("Clip polygon has no spatial reference")

        # Check polygon count and merge if needed
        nPolygon = pLayer_clip.GetFeatureCount()
        if nPolygon == 0:
            raise ValueError("The polygon file does not contain any features!")

        logger.info(f"Clip polygon contains {nPolygon} feature(s)")

        sExtension_clip = os.path.splitext(sFilename_polygon_in)[1]

        if nPolygon > 1:
            if iFlag_endorheic == 1:
                logger.info(
                    "Multiple polygons detected, but endorheic flag is set - skipping merge"
                )
            else:
                pDataset_clip = None
                pLayer_clip = None
                logger.info(
                    "Multiple polygons detected, merging into single feature..."
                )
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
            logger.info("Projection mismatch detected, reprojecting clip polygon...")
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
            logger.info("Clip polygon reprojected")

        # Apply spatial filter using polygon extent
        pEnvelope = pLayer_clip.GetExtent()
        minx, maxx, miny, maxy = pEnvelope
        pLayer_source.SetSpatialFilterRect(minx, miny, maxx, maxy)
        logger.info(f"Spatial filter applied: ({minx}, {miny}) to ({maxx}, {maxy})")

        # Get clip polygon geometry
        pLayer_clip.ResetReading()
        pFeature_clip = pLayer_clip.GetNextFeature()
        if pFeature_clip is None:
            raise RuntimeError("Could not read clip polygon feature")

        pPolygon_clip = pFeature_clip.GetGeometryRef()
        if pPolygon_clip is None:
            raise RuntimeError("Clip polygon feature has no geometry")

        # Create output dataset
        pDataset_clipped = pDriver_vector_out.CreateDataSource(sFilename_vector_out)
        if pDataset_clipped is None:
            raise RuntimeError(f"Could not create output file: {sFilename_vector_out}")

        # Determine output geometry type
        pLayer_source.ResetReading()
        pFeature_in = pLayer_source.GetNextFeature()
        if pFeature_in is None:
            raise RuntimeError("Input vector has no features")

        pGeometry_in = pFeature_in.GetGeometryRef()
        if pGeometry_in is None:
            raise RuntimeError("Input feature has no geometry")

        iGeomType = pGeometry_in.GetGeometryType()

        try:
            output_geom_type = get_output_geometry_type(iGeomType)
            geom_name = ogr.GeometryTypeToName(output_geom_type)
            logger.info(f"Output geometry type: {geom_name}")
        except ValueError as e:
            raise ValueError(str(e))

        # Create output layer
        pLayer_clipped = pDataset_clipped.CreateLayer(
            "layer", pSpatial_reference_target, geom_type=output_geom_type
        )
        if pLayer_clipped is None:
            raise RuntimeError("Could not create output layer")

        # Copy field definitions
        pFeatureDefn_source = pLayer_source.GetLayerDefn()
        for i in range(pFeatureDefn_source.GetFieldCount()):
            pFieldDefn_source = pFeatureDefn_source.GetFieldDefn(i)
            pLayer_clipped.CreateField(pFieldDefn_source)

        logger.info("Processing features for clipping...")
        clipped_count = 0
        pLayer_source.ResetReading()  # Reset to process all features

        # Process each feature
        for pFeature in pLayer_source:
            pGeometry_source = pFeature.GetGeometryRef()
            if pGeometry_source is None:
                continue

            # Check if feature is entirely within clip polygon
            if pPolygon_clip.Contains(pGeometry_source):
                pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                pFeature_clipped.SetGeometry(pGeometry_source)

                # Copy attributes
                for i in range(pFeature.GetFieldCount()):
                    pFeature_clipped.SetField(
                        pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i)
                    )

                pLayer_clipped.CreateFeature(pFeature_clipped)
                clipped_count += 1

            elif pGeometry_source.Intersects(pPolygon_clip):
                # Perform intersection for partially overlapping features
                pGeometry_intersect = pGeometry_source.Intersection(pPolygon_clip)

                if pGeometry_intersect is None or pGeometry_intersect.IsEmpty():
                    logger.debug("Intersection resulted in empty geometry")
                    continue

                iGeomType_intersect = pGeometry_intersect.GetGeometryType()
                sGeomType_intersect = ogr.GeometryTypeToName(iGeomType_intersect)
                logger.debug(f"Intersection geometry type: {sGeomType_intersect}")

                # Create clipped feature
                pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                pFeature_clipped.SetGeometry(pGeometry_intersect)

                # Copy attributes
                for i in range(pFeature.GetFieldCount()):
                    pFeature_clipped.SetField(
                        pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i)
                    )

                pLayer_clipped.CreateFeature(pFeature_clipped)
                clipped_count += 1

        logger.info(f"Clipped {clipped_count} features to output file")
        logger.info(f"Clipping completed successfully: {sFilename_vector_out}")

    except Exception as e:
        logger.error(f"Error during vector clipping: {e}")
        raise RuntimeError(f"Failed to clip vector: {e}")

    finally:
        # Clean up resources
        if "pDataset_source" in locals():
            pDataset_source = None
        if "pDataset_clip" in locals():
            pDataset_clip = None
        if "pDataset_clipped" in locals():
            pDataset_clipped = None
        if "pSpatial_reference_clip" in locals():
            pSpatial_reference_clip = None
        if "pSpatial_reference_target" in locals():
            pSpatial_reference_target = None

    return


def clip_vector_by_polygon_files(
    sFilename_vector_in: str, aFilename_polygon_in: List[str], sFilename_vector_out: str
) -> None:
    """
    Clip a vector file by multiple polygon files (union of all polygons).

    This function clips features from an input vector file using the union of
    multiple polygon files. All polygons from all input files are combined to
    create a single clipping region.

    Args:
        sFilename_vector_in: Path to input vector file (various formats supported)
        aFilename_polygon_in: List of paths to polygon files for clipping
        sFilename_vector_out: Path to output clipped vector file

    Raises:
        FileNotFoundError: If input or polygon files don't exist
        RuntimeError: If vector clipping fails

    Note:
        - Supports multiple vector formats (GeoJSON, Shapefile, GeoPackage, etc.)
        - Automatically handles projection differences
        - Merges multiple polygons within each file if necessary
        - Uses spatial filtering for performance optimization
    """
    logger.info("=" * 80)
    logger.info(f"Starting multi-polygon vector clipping")
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

        # Open input vector dataset
        pDataset_source = ogr.Open(sFilename_vector_in)
        if pDataset_source is None:
            raise RuntimeError(
                f"Could not open input vector file: {sFilename_vector_in}"
            )

        logger.info(
            f"Opened input vector dataset with {pDataset_source.GetLayerCount()} layer(s)"
        )

        # Remove output file if it exists
        if os.path.exists(sFilename_vector_out):
            logger.warning(f"Output file exists, removing: {sFilename_vector_out}")
            os.remove(sFilename_vector_out)

        # Get output driver from extension
        sExtension_out = os.path.splitext(sFilename_vector_out)[1]
        pDriver_out = get_vector_driver_from_extension(sExtension_out)
        if pDriver_out is None:
            raise RuntimeError(f"Could not get driver for extension: {sExtension_out}")

        # Get source layer and spatial reference
        pLayer_source = pDataset_source.GetLayer()
        pSpatial_reference_target = pLayer_source.GetSpatialRef()
        pProjection_target = pSpatial_reference_target.ExportToWkt()
        logger.info(f"Source projection: {pProjection_target[:100]}...")

        # Create output dataset
        pDataset_clipped = pDriver_out.CreateDataSource(sFilename_vector_out)
        if pDataset_clipped is None:
            raise RuntimeError(f"Could not create output file: {sFilename_vector_out}")

        # Determine output geometry type
        pLayer_source.ResetReading()
        pFeature_sample = pLayer_source.GetNextFeature()
        if pFeature_sample is None:
            raise RuntimeError("Source layer contains no features")

        pGeometry_sample = pFeature_sample.GetGeometryRef()
        iGeomType_input = pGeometry_sample.GetGeometryType()
        iGeomType_output = get_output_geometry_type(iGeomType_input)

        if iGeomType_output is None:
            sGeomType = ogr.GeometryTypeToName(iGeomType_input)
            raise RuntimeError(f"Unsupported geometry type: {sGeomType}")

        sGeomType_output = ogr.GeometryTypeToName(iGeomType_output)
        logger.info(f"Output geometry type: {sGeomType_output}")

        # Create output layer
        pLayer_clipped = pDataset_clipped.CreateLayer(
            "layer", pSpatial_reference_target, geom_type=iGeomType_output
        )
        if pLayer_clipped is None:
            raise RuntimeError("Could not create output layer")

        # Copy field definitions
        pFeatureDefn_source = pLayer_source.GetLayerDefn()
        for i in range(pFeatureDefn_source.GetFieldCount()):
            pFieldDefn_source = pFeatureDefn_source.GetFieldDefn(i)
            pLayer_clipped.CreateField(pFieldDefn_source)

        # Process each polygon file
        logger.info("Processing polygon files for clipping...")
        total_clipped = 0

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

            # Handle multiple polygons by merging
            if nPolygon > 1:
                logger.info("Multiple polygons detected, merging into single polygon")
                sExtension_clip = os.path.splitext(sFilename_polygon_in)[1]
                sFilename_clip_merged = sFilename_polygon_in.replace(
                    sExtension_clip, "_merged" + sExtension_clip
                )

                # Close current dataset before merging
                pDataset_clip = None
                pLayer_clip = None

                # Merge features
                from pyearth.toolbox.management.vector.merge_features import (
                    merge_features,
                )

                merge_features(sFilename_polygon_in, sFilename_clip_merged)

                # Reopen merged file
                pDataset_clip = ogr.Open(sFilename_clip_merged)
                pLayer_clip = pDataset_clip.GetLayer()

            # Check projection compatibility
            pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
            pProjection_clip = pSpatial_reference_clip.ExportToWkt()

            sFilename_clip_final = sFilename_polygon_in
            if pProjection_target != pProjection_clip:
                logger.info("Projection mismatch detected, reprojecting polygon file")

                # Close dataset before reprojection
                pDataset_clip = None
                pLayer_clip = None

                # Reproject to target projection
                sFolder = os.path.dirname(sFilename_polygon_in)
                sName = os.path.basename(sFilename_polygon_in)
                sName_no_ext = os.path.splitext(sName)[0]
                sExtension_clip = os.path.splitext(sFilename_polygon_in)[1]
                sFilename_clip_final = os.path.join(
                    sFolder, f"{sName_no_ext}_transformed{sExtension_clip}"
                )

                from pyearth.toolbox.management.vector.reproject import reproject_vector

                reproject_vector(
                    sFilename_polygon_in, sFilename_clip_final, pProjection_target
                )

                # Reopen reprojected file
                pDataset_clip = ogr.Open(sFilename_clip_final)
                pLayer_clip = pDataset_clip.GetLayer()

            # Get clip polygon geometry
            pLayer_clip.ResetReading()
            pFeature_clip = pLayer_clip.GetNextFeature()
            if pFeature_clip is None:
                logger.warning("Could not read clip polygon feature, skipping")
                pDataset_clip = None
                continue

            pPolygon_clip = pFeature_clip.GetGeometryRef()
            if pPolygon_clip is None:
                logger.warning("Clip polygon has no geometry, skipping")
                pDataset_clip = None
                continue

            # Set spatial filter for performance
            pEnvelope = pLayer_clip.GetExtent()
            minx, maxx, miny, maxy = pEnvelope
            pLayer_source.SetSpatialFilterRect(minx, miny, maxx, maxy)
            logger.info(
                f"Applied spatial filter: [{minx:.2f}, {miny:.2f}, {maxx:.2f}, {maxy:.2f}]"
            )

            # Process features
            clipped_count = 0
            pLayer_source.ResetReading()

            for pFeature in pLayer_source:
                pGeometry_source = pFeature.GetGeometryRef()
                if pGeometry_source is None:
                    continue

                # Check if feature is entirely within clip polygon
                if pPolygon_clip.Contains(pGeometry_source):
                    pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                    pFeature_clipped.SetGeometry(pGeometry_source)

                    # Copy attributes
                    for i in range(pFeature.GetFieldCount()):
                        pFeature_clipped.SetField(
                            pFeature.GetFieldDefnRef(i).GetNameRef(),
                            pFeature.GetField(i),
                        )

                    pLayer_clipped.CreateFeature(pFeature_clipped)
                    clipped_count += 1

                elif pGeometry_source.Intersects(pPolygon_clip):
                    # Perform intersection for partially overlapping features
                    pGeometry_intersect = pGeometry_source.Intersection(pPolygon_clip)

                    if pGeometry_intersect is None or pGeometry_intersect.IsEmpty():
                        logger.debug("Intersection resulted in empty geometry")
                        continue

                    # Create clipped feature
                    pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                    pFeature_clipped.SetGeometry(pGeometry_intersect)

                    # Copy attributes
                    for i in range(pFeature.GetFieldCount()):
                        pFeature_clipped.SetField(
                            pFeature.GetFieldDefnRef(i).GetNameRef(),
                            pFeature.GetField(i),
                        )

                    pLayer_clipped.CreateFeature(pFeature_clipped)
                    clipped_count += 1

            logger.info(f"Clipped {clipped_count} features from polygon file {idx}")
            total_clipped += clipped_count

            # Clean up for this iteration
            pDataset_clip = None
            pLayer_clip = None

        logger.info(f"Total features clipped: {total_clipped}")
        logger.info(
            f"Multi-polygon clipping completed successfully: {sFilename_vector_out}"
        )

    except Exception as e:
        logger.error(f"Error during multi-polygon vector clipping: {e}")
        raise RuntimeError(f"Failed to clip vector with multiple polygons: {e}")

    finally:
        # Clean up resources
        if "pDataset_source" in locals():
            pDataset_source = None
        if "pDataset_clip" in locals():
            pDataset_clip = None
        if "pDataset_clipped" in locals():
            pDataset_clipped = None
        if "pSpatial_reference_clip" in locals():
            pSpatial_reference_clip = None
        if "pSpatial_reference_target" in locals():
            pSpatial_reference_target = None

    return
