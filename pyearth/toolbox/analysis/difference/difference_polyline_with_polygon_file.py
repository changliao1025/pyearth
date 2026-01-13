"""
Geometric difference analysis between polylines and polygons using spatial indexing.
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


def difference_polyline_with_polygon_file(
    sFilename_base: str,
    sFilename_new: str,
    sFilename_difference_out: str,
    spatial_ref_epsg: int = 4326,
) -> None:
    """
    Calculate the geometric difference between base polygons and new polygons.

    This function finds areas in the base polygons that are not covered by
    the new polygons and outputs them as difference polygons. Supports multiple
    vector formats including GeoJSON, Shapefile, GeoPackage, KML, and more.

    Args:
        sFilename_base: Path to the base polygon file (supports multiple formats: .geojson, .shp, .gpkg, etc.)
        sFilename_new: Path to the new polygon file (supports multiple formats)
        sFilename_difference_out: Path for the output difference file (format determined by extension)
        spatial_ref_epsg: EPSG code for output spatial reference system (default: 4326 WGS84)

    Returns:
        None

    Raises:
        FileNotFoundError: If input files don't exist
        ImportError: If required spatial indexing libraries are not available
        RuntimeError: If GDAL/OGR operations fail
        ValueError: If file format is not supported

    Example:
        >>> calculate_polyline_polygon_difference(
        ...     'base_polygons.geojson',
        ...     'new_polygons.geojson',
        ...     'difference.shp'
        ... )
    """

    # Setup spatial indexing (rtree only)

    # Validate input files
    if not os.path.exists(sFilename_base):
        raise FileNotFoundError(f"Base file not found: {sFilename_base}")
    if not os.path.exists(sFilename_new):
        raise FileNotFoundError(f"New file not found: {sFilename_new}")

    # Validate output directory
    output_dir = os.path.dirname(sFilename_difference_out)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Created output directory: {output_dir}")

    start_time = datetime.now()
    logger.info(f"Starting polyline-polygon difference calculation at {start_time}")

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

        # Remove existing output file
        if os.path.exists(sFilename_difference_out):
            # For shapefiles, remove associated files too
            if format_out == "ESRI Shapefile":
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

        # Setup spatial reference system
        pSpatial_reference_gcs = osr.SpatialReference()
        pSpatial_reference_gcs.ImportFromEPSG(spatial_ref_epsg)
        logger.info(f"Using spatial reference EPSG:{spatial_ref_epsg}")

        # Create output layer
        pLayerOut = pDataset_out.CreateLayer(
            "diff", pSpatial_reference_gcs, ogr.wkbPolygon
        )
        if pLayerOut is None:
            raise RuntimeError("Could not create output layer")

        # Add attributes
        pLayerOut.CreateField(ogr.FieldDefn("id", ogr.OFTInteger64))
        pLayerOut.CreateField(
            ogr.FieldDefn("base_fid", ogr.OFTInteger64)
        )  # Track source feature
        pLayerOut.CreateField(
            ogr.FieldDefn("new_fid", ogr.OFTInteger64)
        )  # Track intersecting feature

        pLayerDefn = pLayerOut.GetLayerDefn()

        # Open and validate base file
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

        # Build spatial index for base features
        logger.info("Building spatial index for base features...")
        index_base = RTreeindex()

        base_features = {}  # Cache base features to avoid repeated queries
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
                envelope = pGeometry_base.GetEnvelope()
                left, right, bottom, top = envelope

                index_base.insert(fid, (left, bottom, right, top))

                # Cache feature for later use
                base_features[fid] = pFeature_base.Clone()
                indexed_count += 1
            except Exception as e:
                logger.error(f"Error indexing base feature {fid}: {e}")
                continue

        logger.info(f"Indexed {indexed_count} valid base features")

        if not base_features:
            logger.warning("No valid base features found for indexing")
            return

        # Process new features and calculate differences
        logger.info("Processing differences...")
        lID_polygon = 1
        processed_features = 0

        for idx, pFeature_new in enumerate(pLayer_new):
            processed_features += 1
            if processed_features % 100 == 0:
                logger.info(
                    f"Processed {processed_features}/{nFeature_new} new features"
                )

            new_fid = pFeature_new.GetFID()
            pGeometry_new = pFeature_new.GetGeometryRef()

            if pGeometry_new is None:
                logger.warning(f"New feature {new_fid} has no geometry, skipping")
                continue

            # Validate geometry
            if not pGeometry_new.IsValid():
                logger.warning(
                    f"New feature {new_fid} has invalid geometry, attempting to fix"
                )
                pGeometry_new = pGeometry_new.Buffer(0)

            if pGeometry_new is None or pGeometry_new.IsEmpty():
                logger.warning(
                    f"New feature {new_fid} has empty geometry after validation, skipping"
                )
                continue

            try:
                envelope = pGeometry_new.GetEnvelope()
                left, right, bottom, top = envelope

                # Find intersecting base features
                aIntersect = list(index_base.intersection((left, bottom, right, top)))

                for base_fid in aIntersect:
                    if base_fid not in base_features:
                        continue

                    pFeature_base = base_features[base_fid]
                    pGeometry_base = pFeature_base.GetGeometryRef()

                    try:
                        # Check for actual intersection
                        if not pGeometry_new.Intersects(pGeometry_base):
                            continue

                        # Calculate difference
                        pGeometry_diff = pGeometry_base.Difference(pGeometry_new)

                        if pGeometry_diff is None or pGeometry_diff.IsEmpty():
                            continue

                        # Create output feature(s)
                        pGeometrytype_intersect = pGeometry_diff.GetGeometryName()

                        if pGeometrytype_intersect == "POLYGON":
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
                            # Handle multi-part geometries
                            for i in range(pGeometry_diff.GetGeometryCount()):
                                pGeometry_part = pGeometry_diff.GetGeometryRef(i)
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
                            f"Error processing intersection between base {base_fid} and new {new_fid}: {e}"
                        )
                        continue

            except Exception as e:
                logger.error(f"Error processing new feature {new_fid}: {e}")
                continue

        logger.info(f"Generated {lID_polygon - 1} difference polygons")

    except Exception as e:
        logger.error(f"Error during processing: {e}")
        raise RuntimeError(f"Failed to calculate polyline-polygon difference: {e}")

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
