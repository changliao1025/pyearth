import logging
import os
import time
from typing import Optional

from osgeo import gdal, ogr, osr

gdal.UseExceptions()

from pyearth.gis.gdal.gdal_vector_format_support import (
    get_available_vector_formats,
    get_vector_driver_from_format,
    get_vector_format_from_filename,
)


def remove_internal_polygon(
    sFilename_vector_in: str,
    sFilename_vector_out: str,
    verbose: bool = True,
) -> None:
    """
    Remove internal rings (holes) from polygons in a vector file.

    This function processes polygon and multipolygon features from an input vector
    file and removes all interior rings (holes/islands/lakes), keeping only the
    exterior boundary ring for each polygon part. All attribute fields and spatial
    reference information are preserved.

    Parameters
    ----------
    sFilename_vector_in : str
        Path to the input vector file containing polygon or multipolygon geometries.
    sFilename_vector_out : str
        Path where the output vector file without interior rings will be created.
    verbose : bool, optional
        If True, print progress information and statistics (default is True).

    Returns
    -------
    None
        The function creates a new vector file at `sFilename_vector_out`.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If the input filenames are invalid or vector format is unsupported.
    RuntimeError
        If GDAL operations fail (file opening, creation, or processing).
    """
    # Input validation
    if not isinstance(sFilename_vector_in, str) or not sFilename_vector_in.strip():
        raise ValueError("Input filename must be a non-empty string")

    if not isinstance(sFilename_vector_out, str) or not sFilename_vector_out.strip():
        raise ValueError("Output filename must be a non-empty string")

    if not os.path.exists(sFilename_vector_in):
        raise FileNotFoundError(f"Input file does not exist: {sFilename_vector_in}")

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
        available_formats = get_available_vector_formats()
        raise ValueError(
            f"Output format '{sFormat_out}' is not supported. "
            f"Available formats: {available_formats}"
        )

    if verbose:
        logging.info(f"Input format: {sFormat_in}")
        logging.info(f"Output format: {sFormat_out}")

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

        # Get spatial reference from input layer, default to WGS84 if not defined
        pSrs = pLayer_in.GetSpatialRef()
        if pSrs is None:
            pSrs = osr.SpatialReference()
            pSrs.ImportFromEPSG(4326)

        iGeomType = pLayer_in.GetGeomType()
        # Ensure output layer supports multipolygons if present or general polygons
        if iGeomType not in (
            ogr.wkbMultiPolygon,
            ogr.wkbMultiPolygon25D,
            ogr.wkbPolygon,
            ogr.wkbPolygon25D,
        ):
            iGeomType = ogr.wkbMultiPolygon

        nTotal_features = pLayer_in.GetFeatureCount()
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
            pLayer_out = pDataSource_out.CreateLayer(
                "no_internal_rings", pSrs, iGeomType
            )
            if pLayer_out is None:
                raise RuntimeError("Could not create output layer")

            # Copy field definitions
            pLayerDefn_in = pLayer_in.GetLayerDefn()
            nFieldCount = pLayerDefn_in.GetFieldCount()
            for i in range(nFieldCount):
                pFieldDefn = pLayerDefn_in.GetFieldDefn(i)
                try:
                    pLayer_out.CreateField(pFieldDefn)
                except Exception as e:
                    logging.warning(
                        f"Could not create field '{pFieldDefn.GetName()}': {e}"
                    )

            pLayerDefn_out = pLayer_out.GetLayerDefn()

            # Process features
            start_time = time.time()
            nProcessed = 0
            nHolesRemoved = 0

            pLayer_in.ResetReading()

            for pFeature_in in pLayer_in:
                nProcessed += 1
                pGeometry_in = pFeature_in.GetGeometryRef()
                if pGeometry_in is None:
                    continue

                sGeomName = pGeometry_in.GetGeometryName().upper()
                pGeometry_out = None

                if sGeomName == "POLYGON":
                    n_rings = pGeometry_in.GetGeometryCount()
                    if n_rings > 1:
                        nHolesRemoved += n_rings - 1
                        pOuterRing = pGeometry_in.GetGeometryRef(0)
                        pGeometry_out = ogr.Geometry(ogr.wkbPolygon)
                        pGeometry_out.AddGeometry(pOuterRing.Clone())
                    else:
                        pGeometry_out = pGeometry_in.Clone()

                elif sGeomName == "MULTIPOLYGON":
                    pGeometry_out = ogr.Geometry(ogr.wkbMultiPolygon)
                    for j in range(pGeometry_in.GetGeometryCount()):
                        pPart = pGeometry_in.GetGeometryRef(j)
                        if (
                            pPart is not None
                            and pPart.GetGeometryName().upper() == "POLYGON"
                        ):
                            n_rings = pPart.GetGeometryCount()
                            if n_rings > 1:
                                nHolesRemoved += n_rings - 1
                                pOuterRing = pPart.GetGeometryRef(0)
                                pNewPoly = ogr.Geometry(ogr.wkbPolygon)
                                pNewPoly.AddGeometry(pOuterRing.Clone())
                                pGeometry_out.AddGeometry(pNewPoly)
                            else:
                                pGeometry_out.AddGeometry(pPart.Clone())
                        elif pPart is not None:
                            pGeometry_out.AddGeometry(pPart.Clone())

                else:
                    # Non-polygon geometry: pass through unchanged
                    pGeometry_out = pGeometry_in.Clone()

                if pGeometry_out is not None:
                    pGeometry_out.AssignSpatialReference(pSrs)
                    pFeature_out = ogr.Feature(pLayerDefn_out)
                    pFeature_out.SetGeometry(pGeometry_out)

                    for i in range(nFieldCount):
                        try:
                            pFeature_out.SetField(i, pFeature_in.GetField(i))
                        except Exception:
                            pass

                    pLayer_out.CreateFeature(pFeature_out)
                    pFeature_out = None

            total_time = time.time() - start_time
            if verbose:
                logging.info("Processing complete!")
                logging.info(f"Total features processed: {nProcessed}")
                logging.info(f"Total internal rings removed: {nHolesRemoved}")
                logging.info(f"Processing time: {total_time:.2f} seconds")
                logging.info(f"Output saved to: {sFilename_vector_out}")

        finally:
            pDataSource_out = None

    finally:
        pDataSource_in = None
