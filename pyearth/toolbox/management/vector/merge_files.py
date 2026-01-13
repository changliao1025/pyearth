import os, sys
import numpy as np
from typing import List
from osgeo import gdal, osr, ogr, gdalconst
from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_format_from_filename,
    print_supported_vector_formats,
    get_vector_driver_from_format,
    get_vector_driver_from_filename,
)


def merge_files(
    aFilename_in: List[str],
    sFilename_out: str,
    copy_attributes: bool = False,
    add_id_field: bool = True,
    verbose: bool = True,
) -> None:
    """
    Merge multiple vector files into a single output file.

    This function combines multiple vector files of the same geometry type into a single
    output file. All input files must have compatible geometry types (e.g., all points,
    all polygons, etc.). The output format is automatically determined from the output
    filename extension.

    Parameters
    ----------
    aFilename_in : List[str]
        List of input vector filenames to merge. All files must exist and be readable
        by GDAL/OGR. Files should have compatible geometry types.
    sFilename_out : str
        Output filename with extension. Format will be auto-detected from extension.
        Supported formats include GeoJSON (.geojson), Shapefile (.shp), GeoPackage (.gpkg), etc.
    copy_attributes : bool, optional
        If True, copy all attribute fields from input files to output. If False, only
        geometries are copied. Default is False.
    add_id_field : bool, optional
        If True, add a sequential 'id' field to the output features. Default is True.
    verbose : bool, optional
        If True, print progress messages and warnings. Default is True.

    Returns
    -------
    None
        Function modifies files on disk but returns no value.

    Raises
    ------
    FileNotFoundError
        If any input file does not exist (implicit via GDAL).
    RuntimeError
        If output file cannot be created or input files cannot be read (implicit via GDAL).

    Notes
    -----
    1. **Geometry Type Compatibility**: All input files must have the same base geometry
       type (Point, LineString, Polygon, etc.). Multi-geometry types are supported.

    2. **Spatial Reference System**: The SRS from the first input file is used for the
       output. No coordinate transformation is performed.

    3. **Attribute Handling**: When copy_attributes=True, field definitions are copied
       from the first input file. Field name conflicts are not resolved.

    4. **ID Field**: The sequential ID starts from 0 and increments for each feature
       across all input files.

    5. **Output Format Support**: Format availability depends on GDAL/OGR installation.
       Use print_supported_vector_formats() to see available formats.

    6. **File Overwrite**: Existing output files are automatically deleted before
       creating new ones. For Shapefiles, all associated files (.shp, .shx, .dbf, etc.)
       are removed.

    Examples
    --------
    Basic merge without attributes:

    >>> merge_files(
    ...     ['input1.shp', 'input2.shp'],
    ...     'merged_output.shp',
    ...     copy_attributes=False,
    ...     add_id_field=True
    ... )

    Merge with attributes preserved:

    >>> merge_files(
    ...     ['data/polygons1.geojson', 'data/polygons2.geojson'],
    ...     'combined_polygons.gpkg',
    ...     copy_attributes=True,
    ...     verbose=False
    ... )

    See Also
    --------
    merge_features : Merge connected features within a single file
    print_supported_vector_formats : List available vector formats
    """

    # Auto-detect format from output filename if not specified
    sFormat = get_vector_format_from_filename(sFilename_out)
    if sFormat is None:
        print("Could not determine output format from filename extension.")
        print("Supported formats are:")
        print_supported_vector_formats()
        return

    if verbose:
        print(f"=== Starting merge_files operation ===")
        print(f"Input files: {len(aFilename_in)} files")
        print(f"Output file: {sFilename_out}")
        print(f"Copy attributes: {copy_attributes}")
        print(f"Add ID field: {add_id_field}")

    # Register all drivers
    ogr.RegisterAll()

    # Get the driver based on format
    pDriver = ogr.GetDriverByName(sFormat)
    if pDriver is None:
        print(f"{sFormat} driver not available.")
        return  # Check if the output file exists and delete it if it does
    if os.path.exists(sFilename_out):
        # For shapefiles, we need to delete all associated files
        if sFormat == "ESRI Shapefile":
            # Get the base name without extension
            base_name = os.path.splitext(sFilename_out)[0]
            # Delete all shapefile components
            for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg", ".qpj", ".sbn", ".sbx"]:
                file_to_delete = base_name + ext
                if os.path.exists(file_to_delete):
                    os.remove(file_to_delete)
        else:
            pDriver.DeleteDataSource(sFilename_out)

    # Create the output data source
    pDataset_out = pDriver.CreateDataSource(sFilename_out)
    if pDataset_out is None:
        print("Dataset not created")
        return

    # Loop through each input file
    iFlag_first = 1
    lid = 0
    for sFilename_in in aFilename_in:
        # Open the input file
        pDataset_in = ogr.Open(sFilename_in)
        if pDataset_in is None:
            print(f"Failed to open file: {sFilename_in}")
            continue

        # Get the input layer
        pLayer_in = pDataset_in.GetLayer()
        if pLayer_in is None:
            print(f"Failed to get layer from file: {sFilename_in}")
            continue

        # Obtain the geometry type
        iGeomType = pLayer_in.GetGeomType()
        if verbose:
            print(f"Processing {sFilename_in} - Geometry type: {iGeomType}")

        # Create the output layer based on the geometry type
        if iFlag_first == 1:
            # Set spatial reference system from first input layer
            pSpatialRef = pLayer_in.GetSpatialRef()

            if iGeomType == ogr.wkbPoint:
                pLayer_out = pDataset_out.CreateLayer(
                    "layer", srs=pSpatialRef, geom_type=ogr.wkbPoint
                )
            elif iGeomType == ogr.wkbLineString:
                pLayer_out = pDataset_out.CreateLayer(
                    "layer", srs=pSpatialRef, geom_type=ogr.wkbLineString
                )
            elif iGeomType == ogr.wkbPolygon:
                pLayer_out = pDataset_out.CreateLayer(
                    "layer", srs=pSpatialRef, geom_type=ogr.wkbPolygon
                )
            elif iGeomType == ogr.wkbMultiPoint:
                pLayer_out = pDataset_out.CreateLayer(
                    "layer", srs=pSpatialRef, geom_type=ogr.wkbMultiPoint
                )
            elif iGeomType == ogr.wkbMultiLineString:
                pLayer_out = pDataset_out.CreateLayer(
                    "layer", srs=pSpatialRef, geom_type=ogr.wkbMultiLineString
                )
            elif iGeomType == ogr.wkbMultiPolygon:
                pLayer_out = pDataset_out.CreateLayer(
                    "layer", srs=pSpatialRef, geom_type=ogr.wkbMultiPolygon
                )
            else:
                print(f"Unsupported geometry type: {iGeomType}")
                sGeomType = ogr.GeometryTypeToName(iGeomType)
                print("Geometry type not supported:", sGeomType)
                continue

            # Copy field definitions from the first layer (if copy_attributes is True)
            if copy_attributes:
                pLayerDefn_in = pLayer_in.GetLayerDefn()
                for i in range(pLayerDefn_in.GetFieldCount()):
                    pFieldDefn = pLayerDefn_in.GetFieldDefn(i)
                    if pLayer_out.CreateField(pFieldDefn) != 0:
                        print(
                            f"Failed to create field {pFieldDefn.GetName()} in output layer"
                        )
                        return

            # Add an id field if requested and it doesn't exist
            if add_id_field:
                pLayerDefn_in = pLayer_in.GetLayerDefn()
                if not copy_attributes or pLayerDefn_in.GetFieldIndex("id") == -1:
                    pField = ogr.FieldDefn("id", ogr.OFTInteger)
                    if pLayer_out.CreateField(pField) != 0:
                        print("Failed to create id field in output layer")
                        return

            iFlag_first = 0
        else:
            pass

        # Copy features from the input layer to the output layer
        nFeatures = pLayer_in.GetFeatureCount()
        nProcessed = 0
        nSkipped = 0

        if verbose:
            print(f"Processing {nFeatures} features from {sFilename_in}")

        for feature in pLayer_in:
            nProcessed += 1
            if verbose and nProcessed % 1000 == 0:
                print(f"  Processed {nProcessed}/{nFeatures} features")

            # Copy the feature
            pFeature_new = ogr.Feature(pLayer_out.GetLayerDefn())

            # Get the geometry
            original_geom = feature.GetGeometryRef()
            if original_geom is not None:
                # Set geometry directly without validation
                pFeature_new.SetGeometry(original_geom)
            else:
                if verbose:
                    print(f"Warning: Feature {nProcessed} has no geometry, skipping")
                nSkipped += 1
                continue

            # Copy field values from the original feature (if copy_attributes is True)
            if copy_attributes:
                pLayerDefn_in = pLayer_in.GetLayerDefn()
                pLayerDefn_out = pLayer_out.GetLayerDefn()

                for i in range(pLayerDefn_in.GetFieldCount()):
                    pFieldDefn_in = pLayerDefn_in.GetFieldDefn(i)
                    sFieldName = pFieldDefn_in.GetName()

                    # Check if the field exists in the output layer
                    iFieldIndex_out = pLayerDefn_out.GetFieldIndex(sFieldName)
                    if iFieldIndex_out >= 0:
                        pFeature_new.SetField(sFieldName, feature.GetField(sFieldName))

            # Set id for each feature (if add_id_field is True)
            if add_id_field:
                pFeature_new.SetField("id", lid)
            pLayer_out.CreateFeature(pFeature_new)
            lid = lid + 1

        # Clean up
        pDataset_in = None

        # Print summary for this file
        if verbose:
            print(
                f"Completed {sFilename_in}: {nProcessed} processed, {nSkipped} skipped"
            )

    # Clean up
    pDataset_out = None
    if verbose:
        print("=== Merge operation completed ===")
        print(f"Total features processed: {lid}")
        print("Merge completed successfully.")

    return
