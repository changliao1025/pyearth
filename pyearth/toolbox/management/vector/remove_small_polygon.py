import os
import logging
import time
from typing import Optional, Union
import numpy as np
from osgeo import ogr, osr, gdal

gdal.UseExceptions()

from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_format_from_filename,
    print_supported_vector_formats,
    get_vector_driver_from_format,
)


def remove_small_polygon(
    sFilename_vector_in: str,
    sFilename_vector_out: str,
    dThreshold_in: Union[float, int],
    iFlag_algorithm: int = 2,
    verbose: bool = True,
    progress_interval: int = 1000,
) -> None:
    """
    Remove small polygons from a vector file based on area threshold.

    This function filters polygon geometries from an input vector file, keeping only
    those with areas greater than the specified threshold. The area calculation uses
    geodesic methods for accurate results on large regions. Both POLYGON and
    MULTIPOLYGON geometries are supported, with holes (inner rings) preserved.

    Parameters
    ----------
    sFilename_vector_in : str
        Path to the input vector file containing polygon geometries.
    sFilename_vector_out : str
        Path where the filtered output vector file will be created.
    dThreshold_in : float or int
        Minimum area threshold in square kilometers. Polygons with areas
        less than or equal to this value will be removed.
    iFlag_algorithm : int, optional
        Algorithm flag for area calculation (default is 2 for geodesic).
        - 1: Planar area calculation
        - 2: Geodesic area calculation (recommended for large areas)
    verbose : bool, optional
        If True, print progress information and statistics (default is True).
    progress_interval : int, optional
        Number of features to process before printing progress updates
        (default is 1000).

    Returns
    -------
    None
        The function creates a new vector file at `sFilename_vector_out`
        containing only polygons that meet the area threshold.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If the threshold value is invalid or if vector format is unsupported.
    RuntimeError
        If GDAL operations fail (file opening, creation, or processing).

    Notes
    -----
    - The function assumes input geometries are in WGS84 (EPSG:4326) coordinate system
      for accurate geodesic area calculations.
    - Output file will contain additional fields: 'id' (sequential integer) and
      'area' (calculated area in square kilometers).
    - All original attribute fields are preserved except for 'id' and 'area' which
      are overwritten.
    - Degenerate polygons (less than 3 vertices) are automatically skipped.
    - For MULTIPOLYGON geometries, each individual polygon part is evaluated
      separately against the threshold.

    Examples
    --------
    Remove polygons smaller than 1 square kilometer:

    >>> remove_small_polygon(
    ...     'input_polygons.shp',
    ...     'filtered_polygons.shp',
    ...     1.0
    ... )

    Process with custom progress reporting:

    >>> remove_small_polygon(
    ...     'world_countries.geojson',
    ...     'large_countries.geojson',
    ...     1000.0,
    ...     verbose=True,
    ...     progress_interval=500
    ... )
    """
    # Input validation
    if not isinstance(sFilename_vector_in, str) or not sFilename_vector_in.strip():
        raise ValueError("Input filename must be a non-empty string")

    if not isinstance(sFilename_vector_out, str) or not sFilename_vector_out.strip():
        raise ValueError("Output filename must be a non-empty string")

    if not os.path.exists(sFilename_vector_in):
        raise FileNotFoundError(f"Input file does not exist: {sFilename_vector_in}")

    try:
        dThreshold = float(dThreshold_in)
    except (TypeError, ValueError) as e:
        raise ValueError(
            f"Threshold must be a numeric value, got {type(dThreshold_in)}: {e}"
        )

    if dThreshold <= 0:
        raise ValueError(f"Threshold must be positive, got {dThreshold}")

    if iFlag_algorithm not in [1, 2]:
        raise ValueError(f"Algorithm flag must be 1 or 2, got {iFlag_algorithm}")

    # Remove existing output file if it exists
    if os.path.exists(sFilename_vector_out):
        if verbose:
            logging.info(f"Removing existing output file: {sFilename_vector_out}")
        try:
            os.remove(sFilename_vector_out)
        except OSError as e:
            raise RuntimeError(f"Could not remove existing output file: {e}")

    # Determine formats and drivers
    try:
        sFormat_in = get_vector_format_from_filename(sFilename_vector_in)
        sFormat_out = get_vector_format_from_filename(sFilename_vector_out)
        pDriver_out = get_vector_driver_from_format(sFormat_out)
    except Exception as e:
        raise RuntimeError(f"Failed to determine vector format: {e}")

    if pDriver_out is None:
        available_formats = print_supported_vector_formats()
        raise ValueError(
            f"Output format '{sFormat_out}' is not supported. "
            f"Available formats: {available_formats}"
        )

    if verbose:
        logging.info(f"Input format: {sFormat_in}")
        logging.info(f"Output format: {sFormat_out}")
        logging.info(f"Area threshold: {dThreshold} m²")
        logging.info(
            f"Area calculation algorithm: {'Geodesic' if iFlag_algorithm == 2 else 'Planar'}"
        )

    # Set up spatial reference (WGS84 for geodesic calculations)
    try:
        pSrs = osr.SpatialReference()
        pSrs.ImportFromEPSG(4326)
    except Exception as e:
        raise RuntimeError(f"Failed to create spatial reference: {e}")

    # Open input datasource
    try:
        pDataSource_in = ogr.Open(sFilename_vector_in, 0)  # Read-only
        if pDataSource_in is None:
            raise RuntimeError(f"Could not open input file: {sFilename_vector_in}")
    except Exception as e:
        raise RuntimeError(f"GDAL error opening input file: {e}")

    try:
        pLayer_in = pDataSource_in.GetLayer()
        if pLayer_in is None:
            raise RuntimeError("Could not access input layer")

        # Get feature count for progress tracking
        nTotal_features = pLayer_in.GetFeatureCount()
        if nTotal_features < 0:
            logging.warning(
                "Could not determine total feature count, progress reporting may be inaccurate"
            )

        if verbose:
            logging.info(f"Processing {nTotal_features} features...")

        # Create output datasource
        try:
            pDataSource_out = pDriver_out.CreateDataSource(sFilename_vector_out)
            if pDataSource_out is None:
                raise RuntimeError(
                    f"Could not create output file: {sFilename_vector_out}"
                )
        except Exception as e:
            raise RuntimeError(f"GDAL error creating output file: {e}")

        try:
            # Create output layer
            pLayer_out = pDataSource_out.CreateLayer(
                "filtered_polygons", pSrs, ogr.wkbPolygon
            )
            if pLayer_out is None:
                raise RuntimeError("Could not create output layer")

            # Add standard fields
            field_id = ogr.FieldDefn("id", ogr.OFTInteger)
            field_area = ogr.FieldDefn("area", ogr.OFTReal)
            pLayer_out.CreateField(field_id)
            pLayer_out.CreateField(field_area)

            # Copy existing fields from input (excluding id and area if they exist)
            pLayerDefn_in = pLayer_in.GetLayerDefn()
            nFieldCount = pLayerDefn_in.GetFieldCount()

            for i in range(nFieldCount):
                pFieldDefn = pLayerDefn_in.GetFieldDefn(i)
                field_name = pFieldDefn.GetName().lower()
                if field_name not in ["id", "area"]:
                    try:
                        pLayer_out.CreateField(pFieldDefn)
                    except Exception as e:
                        logging.warning(f"Could not create field '{field_name}': {e}")

            pLayerDefn_out = pLayer_out.GetLayerDefn()

            # Process features
            start_time = time.time()
            lID = 1
            nProcessed = 0
            nKept = 0
            nSkipped = 0

            pLayer_in.ResetReading()  # Ensure we're at the beginning

            for pFeature_in in pLayer_in:
                nProcessed += 1

                if (
                    verbose
                    and progress_interval > 0
                    and nProcessed % progress_interval == 0
                ):
                    elapsed = time.time() - start_time
                    rate = nProcessed / elapsed if elapsed > 0 else 0
                    logging.info(
                        f"Processed {nProcessed}/{nTotal_features} features "
                        f"({rate:.1f} features/sec), kept {nKept}"
                    )

                pGeometry = pFeature_in.GetGeometryRef()
                if pGeometry is None:
                    nSkipped += 1
                    continue

                geometry_type = pGeometry.GetGeometryName()

                # Process POLYGON geometries
                if geometry_type == "POLYGON":
                    kept = _process_single_polygon(
                        pGeometry,
                        dThreshold,
                        iFlag_algorithm,
                        pSrs,
                        pLayerDefn_out,
                        pFeature_in,
                        pLayerDefn_in,
                        lID,
                        pLayer_out,
                    )
                    if kept:
                        lID += 1
                        nKept += 1

                # Process MULTIPOLYGON geometries
                elif geometry_type == "MULTIPOLYGON":
                    polygons_kept = _process_multipolygon(
                        pGeometry,
                        dThreshold,
                        iFlag_algorithm,
                        pSrs,
                        pLayerDefn_out,
                        pFeature_in,
                        pLayerDefn_in,
                        lID,
                        pLayer_out,
                    )
                    lID += polygons_kept
                    nKept += polygons_kept

                else:
                    nSkipped += 1
                    if verbose and nSkipped <= 5:  # Log first few non-polygon features
                        logging.debug(f"Skipping non-polygon geometry: {geometry_type}")

            # Final statistics
            total_time = time.time() - start_time

            if verbose:
                logging.info("Processing complete!")
                logging.info(f"Total features processed: {nProcessed}")
                logging.info(f"Features kept (area > {dThreshold} km²): {nKept}")
                logging.info(f"Features removed: {nProcessed - nKept - nSkipped}")
                logging.info(f"Features skipped (non-polygon): {nSkipped}")
                logging.info(f"Processing time: {total_time:.2f} seconds")
                logging.info(
                    f"Average processing rate: {nProcessed/total_time:.1f} features/sec"
                )
                logging.info(f"Output saved to: {sFilename_vector_out}")

        finally:
            # Clean up output datasource
            pDataSource_out = None

    finally:
        # Clean up input datasource
        pDataSource_in = None


def _process_single_polygon(
    pGeometry: ogr.Geometry,
    dThreshold: float,
    iFlag_algorithm: int,
    pSrs: osr.SpatialReference,
    pLayerDefn_out: ogr.FeatureDefn,
    pFeature_in: ogr.Feature,
    pLayerDefn_in: ogr.FeatureDefn,
    lID: int,
    pLayer_out: ogr.Layer,
) -> bool:
    """
    Process a single POLYGON geometry.

    Returns True if the polygon was kept, False otherwise.
    """
    try:
        # Extract outer ring coordinates
        pOuterRing = pGeometry.GetGeometryRef(0)
        if pOuterRing is None or pOuterRing.GetPointCount() < 3:
            return False

        aCoords_outer = np.array(
            [
                [pOuterRing.GetPoint(i)[0], pOuterRing.GetPoint(i)[1]]
                for i in range(pOuterRing.GetPointCount())
            ]
        )

        # Calculate area
        dArea = calculate_polygon_area(
            aCoords_outer[:, 0], aCoords_outer[:, 1], iFlag_algorithm
        )
        dAreakm = dArea / 1e6  # Convert to square kilometers

        if dArea <= dThreshold:
            return False

        # Create output polygon with all rings
        pGeometry_out = ogr.Geometry(ogr.wkbPolygon)

        # Add outer ring
        pRing_outer = ogr.Geometry(ogr.wkbLinearRing)
        for coord in aCoords_outer:
            pRing_outer.AddPoint(coord[0], coord[1])
        pRing_outer.CloseRings()
        pGeometry_out.AddGeometry(pRing_outer)

        # Add inner rings (holes)
        for iRing in range(1, pGeometry.GetGeometryCount()):
            pInnerRing_src = pGeometry.GetGeometryRef(iRing)
            if pInnerRing_src is None:
                continue

            pInnerRing_new = ogr.Geometry(ogr.wkbLinearRing)
            for iPoint in range(pInnerRing_src.GetPointCount()):
                x, y, z = pInnerRing_src.GetPoint(iPoint)
                pInnerRing_new.AddPoint(x, y)
            pInnerRing_new.CloseRings()
            pGeometry_out.AddGeometry(pInnerRing_new)

        pGeometry_out.AssignSpatialReference(pSrs)

        # Create and populate output feature
        pFeature_out = ogr.Feature(pLayerDefn_out)
        pFeature_out.SetGeometry(pGeometry_out)
        pFeature_out.SetField("id", lID)
        pFeature_out.SetField("area", dAreakm)

        # Copy other fields
        for i in range(pLayerDefn_in.GetFieldCount()):
            field_name = pLayerDefn_in.GetFieldDefn(i).GetName()
            if field_name.lower() not in ["id", "area"]:
                try:
                    pFeature_out.SetField(field_name, pFeature_in.GetField(field_name))
                except Exception:
                    # Skip fields that can't be copied
                    pass

        pLayer_out.CreateFeature(pFeature_out)
        pFeature_out = None

        return True

    except Exception as e:
        logging.warning(f"Error processing polygon: {e}")
        return False


def _process_multipolygon(
    pGeometry: ogr.Geometry,
    dThreshold: float,
    iFlag_algorithm: int,
    pSrs: osr.SpatialReference,
    pLayerDefn_out: ogr.FeatureDefn,
    pFeature_in: ogr.Feature,
    pLayerDefn_in: ogr.FeatureDefn,
    lID_start: int,
    pLayer_out: ogr.Layer,
) -> int:
    """
    Process a MULTIPOLYGON geometry, creating separate features for each polygon part.

    Returns the number of polygons that were kept.
    """
    polygons_kept = 0

    try:
        for iPoly in range(pGeometry.GetGeometryCount()):
            pPolygon = pGeometry.GetGeometryRef(iPoly)
            if pPolygon is None or pPolygon.GetGeometryName() != "POLYGON":
                continue

            kept = _process_single_polygon(
                pPolygon,
                dThreshold,
                iFlag_algorithm,
                pSrs,
                pLayerDefn_out,
                pFeature_in,
                pLayerDefn_in,
                lID_start + polygons_kept,
                pLayer_out,
            )

            if kept:
                polygons_kept += 1

    except Exception as e:
        logging.warning(f"Error processing multipolygon: {e}")

    return polygons_kept
