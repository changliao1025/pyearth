"""
Clip vector files by a bounding box with support for multiple vector formats.
"""

import logging
import os
from typing import Dict, List, Optional

from osgeo import ogr, osr

from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.gis.gdal.read.vector.gdal_get_vector_boundary import (
    gdal_get_vector_boundary,
)
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.gis.geometry.get_output_geometry_type import get_output_geometry_type
from pyearth.system.define_global_variables import *
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.management.vector.reproject import reproject_vector

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def clip_vector_by_bounding_box(
    sFilename_vector_in: str,
    sFilename_polygon_in: str,
    sFilename_vector_out: str,
    dBuffer_in=None,
) -> None:
    """
    Select features from a vector file using the bounding box of a polygon file.

    Instead of clipping against the individual polygon shape, this function derives the
    outer bounding box (extent) of the polygon file and uses that rectangle as the
    selection region. Any input feature (points, lines, or polygons) that touches, is
    within, or intersects the bounding box is written to the output as-is; the original
    geometry is preserved and NOT trimmed to the box edges.

    Parameters:
    sFilename_vector_in (str): The input vector file to be clipped.
    sFilename_polygon_in (str): The polygon file whose bounding box defines the clip extent.
    sFilename_vector_out (str): The output vector file after clipping.
    iFlag_endorheic (int, optional): A flag indicating whether to clip the vector file
        to the endorheic (1) or exorheic (0) area. Default is 0.

    Returns:
    None
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
            pDriver_vector_in = get_vector_driver_from_filename(sFilename_vector_in)
            pDriver_vector_out = get_vector_driver_from_filename(sFilename_vector_out)
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

        # Use the outer bounding box extent of the clip polygon(s) as the clip region.
        # The extent is read from the (possibly reprojected) clip layer so that it is
        # always expressed in the source projection.
        pEnvelope = pLayer_clip.GetExtent()
        minx, maxx, miny, maxy = pEnvelope

        # add buffer if user provided buffer
        if dBuffer_in is not None:
            minx = minx - dBuffer_in
            maxx = maxx + dBuffer_in
            miny = miny - dBuffer_in
            maxy = maxy + dBuffer_in

        # Build a rectangular polygon from the bounding box extent
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(minx, miny)
        ring.AddPoint(maxx, miny)
        ring.AddPoint(maxx, maxy)
        ring.AddPoint(minx, maxy)
        ring.AddPoint(minx, miny)  # close the ring
        pPolygon_clip = ogr.Geometry(ogr.wkbPolygon)
        pPolygon_clip.AddGeometry(ring)

        # Apply spatial filter using the bounding box extent for performance
        pLayer_source.SetSpatialFilterRect(minx, miny, maxx, maxy)
        logger.info(f"Spatial filter applied: ({minx}, {miny}) to ({maxx}, {maxy})")

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

        # Process each feature.
        # Any feature that touches, is within, or intersects the bounding box is kept
        # as-is (the original geometry is saved without trimming to the box edges).
        for pFeature in pLayer_source:
            pGeometry_source = pFeature.GetGeometryRef()
            if pGeometry_source is None:
                continue

            # Intersects covers features that are within, touching, or overlapping the box
            if pGeometry_source.Intersects(pPolygon_clip):
                pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                pFeature_clipped.SetGeometry(pGeometry_source)

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
