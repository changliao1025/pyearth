"""
Clip raster by polygon file with support for multiple vector formats.
"""

import os
import sys
import logging
from typing import Optional, Tuple
import numpy as np
from osgeo import gdal, ogr, osr

from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_format_from_filename

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def clip_raster_by_polygon_file(
    sFilename_raster_in: str,
    sFilename_polygon_in: str,
    sFilename_raster_out: str,
    iFlag_use_raster_extent: int = 0,
    sFormat_in: str = "GTiff",
    missing_value: float = -9999,
    resample_algorithm: str = "near",
) -> None:
    """
    Clip a raster by a polygon file with support for multiple vector formats.

    This function clips an input raster to the extent of a polygon. If the polygon
    file contains multiple features, they will be automatically merged. Supports
    automatic reprojection if the polygon and raster have different coordinate systems.

    Args:
        sFilename_raster_in: Path to input raster file
        sFilename_polygon_in: Path to input polygon file (supports .geojson, .shp, .gpkg, etc.)
        sFilename_raster_out: Path to output clipped raster file
        iFlag_use_raster_extent: If 1, use raster extent; if 0, crop to polygon extent (default: 0)
        sFormat_in: Input raster format driver name (default: 'GTiff')
        missing_value: NoData value for the output raster (default: -9999)
        resample_algorithm: Resampling algorithm: 'near', 'bilinear', 'cubic', etc. (default: 'near')

    Returns:
        None

    Raises:
        FileNotFoundError: If input files don't exist
        RuntimeError: If GDAL/OGR operations fail
        ValueError: If polygon is out of raster bounds

    Example:
        >>> clip_raster_by_polygon_file(
        ...     'input_raster.tif',
        ...     'clip_boundary.geojson',
        ...     'output_clipped.tif',
        ...     iFlag_use_raster_extent=0
        ... )
    """

    logger.info(f"Starting raster clipping operation")
    logger.info(f"Input raster: {sFilename_raster_in}")
    logger.info(f"Clip polygon: {sFilename_polygon_in}")
    logger.info(f"Output raster: {sFilename_raster_out}")

    # Validate input files
    if not os.path.exists(sFilename_raster_in):
        raise FileNotFoundError(f"Raster file not found: {sFilename_raster_in}")

    if not os.path.exists(sFilename_polygon_in):
        raise FileNotFoundError(f"Polygon file not found: {sFilename_polygon_in}")

    # Remove existing output file
    if os.path.exists(sFilename_raster_out):
        os.remove(sFilename_raster_out)
        logger.info(f"Removed existing output file: {sFilename_raster_out}")

    try:
        # Get raster driver
        if sFormat_in is not None:
            sDriverName = sFormat_in
        else:
            sDriverName = "GTiff"

        pDriver_raster = gdal.GetDriverByName(sDriverName)
        if pDriver_raster is None:
            raise RuntimeError(f"Raster driver '{sDriverName}' not available")

        # Get vector driver using multi-format support
        try:
            pDriver_vector = get_vector_format_from_filename(sFilename_polygon_in)
            logger.info(f"Using vector driver for polygon file")
        except ValueError as e:
            raise RuntimeError(f"Unsupported polygon file format: {e}")

        # Read raster data and metadata
        logger.info("Reading raster metadata...")
        dummy = gdal_read_geotiff_file(sFilename_raster_in)
        aData = dummy["dataOut"]
        eType = dummy["dataType"]
        dPixelWidth = dummy["pixelWidth"]
        dPixelHeight = dummy["pixelHeight"]
        dOriginX = dummy["originX"]
        dOriginY = dummy["originY"]
        nrow = dummy["nrow"]
        ncolumn = dummy["ncolumn"]
        dMissing_value = dummy["missingValue"]

        if dMissing_value is None:
            dMissing_value = missing_value

        eType_out = gdal.GDT_Int16  # Support missing value
        pProjection_target = dummy["projection"]
        min_x, max_x, min_y, max_y = gdal_get_raster_extent(sFilename_raster_in)

        aRaster_extent = [min_x, min_y, max_x, max_y]
        pSpatial_reference_target = osr.SpatialReference()
        pSpatial_reference_target.ImportFromWkt(pProjection_target)

        dX_left = dOriginX
        dX_right = dOriginX + ncolumn * dPixelWidth
        dY_top = dOriginY
        dY_bot = dOriginY + nrow * dPixelHeight

        logger.info(f"Raster extent: ({min_x}, {min_y}) to ({max_x}, {max_y})")
        logger.info(f"Raster size: {ncolumn} x {nrow} pixels")

        # Open and validate polygon file
        pDataset_clip = ogr.Open(sFilename_polygon_in)
        if pDataset_clip is None:
            raise RuntimeError(f"Could not open polygon file: {sFilename_polygon_in}")

        pLayer_clip = pDataset_clip.GetLayer(0)
        if pLayer_clip is None:
            raise RuntimeError("Could not access layer in polygon file")

        nPolygon = pLayer_clip.GetFeatureCount()
        if nPolygon == 0:
            raise ValueError("The polygon file does not contain any features!")

        logger.info(f"Polygon file contains {nPolygon} feature(s)")

        # Handle multiple polygons by merging
        sExtension_vector = os.path.splitext(sFilename_polygon_in)[1]
        if nPolygon > 1:
            pDataset_clip = None
            pLayer_clip = None
            logger.info("Multiple polygons detected, merging into single feature...")
            sFilename_clip_new = sFilename_polygon_in.replace(
                sExtension_vector, "_merged" + sExtension_vector
            )
            merge_features(sFilename_polygon_in, sFilename_clip_new)
            sFilename_polygon_in = sFilename_clip_new

            # Open merged file
            pDataset_clip = ogr.Open(sFilename_polygon_in)
            pLayer_clip = pDataset_clip.GetLayer(0)

        # Get spatial reference and geometry
        pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
        pProjection_clip = pSpatial_reference_clip.ExportToWkt()
        logger.info(
            f"Polygon spatial reference: {pSpatial_reference_clip.GetName() if pSpatial_reference_clip else 'Unknown'}"
        )

        pLayer_clip.ResetReading()
        pFeature_clip = pLayer_clip.GetNextFeature()
        if pFeature_clip is None:
            raise RuntimeError("Could not read polygon feature")

        pPolygon = pFeature_clip.GetGeometryRef()
        if pPolygon is None:
            raise RuntimeError("Polygon feature has no geometry")

        # Handle projection mismatch
        if pProjection_target != pProjection_clip:
            logger.info(
                "Projection mismatch detected, transforming polygon to raster projection..."
            )
            pSRS_source = pSpatial_reference_clip
            pSRS_target = pSpatial_reference_target

            # Create coordinate transformation
            pCoordinateTransform = osr.CoordinateTransformation(
                pSRS_source, pSRS_target
            )

            # Clone and transform polygon
            pPolygon_transformed = pPolygon.Clone()
            pPolygon_transformed.Transform(pCoordinateTransform)
            pPolygon = pPolygon_transformed
            logger.info("Polygon transformed to raster projection")

        # Check if polygon intersects with raster
        minX, maxX, minY, maxY = pPolygon.GetEnvelope()
        logger.info(f"Polygon extent: ({minX}, {minY}) to ({maxX}, {maxY})")

        if minX > dX_right or maxX < dX_left or minY > dY_top or maxY < dY_bot:
            raise ValueError("Polygon is completely outside the raster extent!")

        # Prepare for clipping
        pPolygonWKT = pPolygon.ExportToWkt()
        create_options = ["COMPRESS=DEFLATE", "PREDICTOR=2"]

        # Set up warp options
        if iFlag_use_raster_extent == 1:
            logger.info("Clipping with original raster extent...")
            pWrapOption = gdal.WarpOptions(
                cropToCutline=False,
                cutlineWKT=pPolygonWKT,
                xRes=dPixelWidth,
                yRes=abs(dPixelHeight),
                outputBounds=aRaster_extent,
                dstSRS=pSpatial_reference_target,
                format="MEM",
                resampleAlg=resample_algorithm,
                dstNodata=dMissing_value,
                outputType=eType_out,
            )
        else:
            logger.info("Clipping to polygon extent...")
            pWrapOption = gdal.WarpOptions(
                cropToCutline=True,
                cutlineWKT=pPolygonWKT,
                xRes=dPixelWidth,
                yRes=abs(dPixelHeight),
                dstSRS=pSpatial_reference_target,
                format="MEM",
                resampleAlg=resample_algorithm,
                dstNodata=dMissing_value,
                outputType=eType_out,
            )

        # Perform clipping
        logger.info("Warping raster...")
        pDataset_clip_warped = gdal.Warp("", sFilename_raster_in, options=pWrapOption)
        if pDataset_clip_warped is None:
            raise RuntimeError("Raster warping failed")

        newGeoTransform = pDataset_clip_warped.GetGeoTransform()

        # Convert warped dataset to array
        aData_clip = pDataset_clip_warped.ReadAsArray()
        if aData_clip is None:
            raise RuntimeError("Failed to read clipped raster data")

        logger.info(f"Clipped data min value: {np.min(aData_clip)}")

        iNewWidth = aData_clip.shape[1]
        iNewHeigh = aData_clip.shape[0]
        logger.info(f"Clipped raster size: {iNewWidth} x {iNewHeigh} pixels")

        # Create output raster
        pDataset_clip = pDriver_raster.Create(
            sFilename_raster_out,
            iNewWidth,
            iNewHeigh,
            1,
            eType_out,
            options=create_options,
        )
        if pDataset_clip is None:
            raise RuntimeError(f"Could not create output file: {sFilename_raster_out}")

        pDataset_clip.SetGeoTransform(newGeoTransform)
        pDataset_clip.SetProjection(pProjection_target)

        # Change missing value to standard value
        aData_clip[aData_clip == dMissing_value] = missing_value

        # Write output
        pBand = pDataset_clip.GetRasterBand(1)
        pBand.SetNoDataValue(missing_value)
        pBand.WriteArray(aData_clip)
        pBand.FlushCache()
        pDataset_clip.FlushCache()

        logger.info(f"Successfully clipped raster saved to: {sFilename_raster_out}")

    except Exception as e:
        logger.error(f"Error during raster clipping: {e}")
        raise RuntimeError(f"Failed to clip raster: {e}")

    finally:
        # Clean up resources
        if "pDataset_clip" in locals():
            pDataset_clip = None
        if "pDataset_clip_warped" in locals():
            pDataset_clip_warped = None
        if "pSpatial_reference_target" in locals():
            pSpatial_reference_target = None
        if "pSpatial_reference_clip" in locals():
            pSpatial_reference_clip = None

    return
