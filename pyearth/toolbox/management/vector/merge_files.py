import os, sys
import numpy as np
from osgeo import gdal, osr, ogr, gdalconst

def get_format_from_extension(filename):
    """
    Determine the OGR format string from file extension.

    Args:
        filename: Input filename with extension

    Returns:
        Format string for OGR driver
    """
    _, ext = os.path.splitext(filename.lower())

    format_map = {
        '.geojson': 'GeoJSON',
        '.json': 'GeoJSON',
        '.shp': 'ESRI Shapefile',
        '.gpkg': 'GPKG',
        '.kml': 'KML',
        '.gml': 'GML',
        '.sqlite': 'SQLite'
    }

    return format_map.get(ext, 'GeoJSON')



def merge_files(aFilename_in, sFilename_out, sFormat=None, copy_attributes=False, add_id_field=True, verbose=True):
    """
    Merge multiple vector files into a single output file.

    Args:
        aFilename_in: List of input filenames
        sFilename_out: Output filename
        sFormat: Output format ('GeoJSON', 'ESRI Shapefile', 'GPKG', etc.)
                If None, format will be determined from output file extension
        copy_attributes: If True, copy all attributes from input files. If False, only copy geometries
        add_id_field: If True, add an 'id' field with sequential numbering
        verbose: If True, print progress and warning messages
    """

    # Auto-detect format from output filename if not specified
    if sFormat is None:
        sFormat = get_format_from_extension(sFilename_out)
        if verbose:
            print(f'Auto-detected format: {sFormat}')

    if verbose:
        print(f'=== Starting merge_files operation ===')
        print(f'Input files: {len(aFilename_in)} files')
        print(f'Output file: {sFilename_out}')
        print(f'Copy attributes: {copy_attributes}')
        print(f'Add ID field: {add_id_field}')

    # Register all drivers
    ogr.RegisterAll()

    # Get the driver based on format
    pDriver = ogr.GetDriverByName(sFormat)
    if pDriver is None:
        print(f'{sFormat} driver not available.')
        return    # Check if the output file exists and delete it if it does
    if os.path.exists(sFilename_out):
        # For shapefiles, we need to delete all associated files
        if sFormat == 'ESRI Shapefile':
            # Get the base name without extension
            base_name = os.path.splitext(sFilename_out)[0]
            # Delete all shapefile components
            for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj', '.sbn', '.sbx']:
                file_to_delete = base_name + ext
                if os.path.exists(file_to_delete):
                    os.remove(file_to_delete)
        else:
            pDriver.DeleteDataSource(sFilename_out)

    # Create the output data source
    pDataset_out = pDriver.CreateDataSource(sFilename_out)
    if pDataset_out is None:
        print('Dataset not created')
        return

    # Loop through each input file
    iFlag_first = 1
    lid = 0
    for sFilename_in in aFilename_in:
        # Open the input file
        pDataset_in = ogr.Open(sFilename_in)
        if pDataset_in is None:
            print(f'Failed to open file: {sFilename_in}')
            continue

        # Get the input layer
        pLayer_in = pDataset_in.GetLayer()
        if pLayer_in is None:
            print(f'Failed to get layer from file: {sFilename_in}')
            continue

        # Obtain the geometry type
        iGeomType = pLayer_in.GetGeomType()
        if verbose:
            print(f'Processing {sFilename_in} - Geometry type: {iGeomType}')

        # Create the output layer based on the geometry type
        if iFlag_first == 1:
            # Set spatial reference system from first input layer
            pSpatialRef = pLayer_in.GetSpatialRef()

            if iGeomType == ogr.wkbPoint:
                pLayer_out = pDataset_out.CreateLayer('layer', srs=pSpatialRef, geom_type=ogr.wkbPoint)
            elif iGeomType == ogr.wkbLineString:
                pLayer_out = pDataset_out.CreateLayer('layer', srs=pSpatialRef, geom_type=ogr.wkbLineString)
            elif iGeomType == ogr.wkbPolygon:
                pLayer_out = pDataset_out.CreateLayer('layer', srs=pSpatialRef, geom_type=ogr.wkbPolygon)
            elif iGeomType == ogr.wkbMultiPoint:
                pLayer_out = pDataset_out.CreateLayer('layer', srs=pSpatialRef, geom_type=ogr.wkbMultiPoint)
            elif iGeomType == ogr.wkbMultiLineString:
                pLayer_out = pDataset_out.CreateLayer('layer', srs=pSpatialRef, geom_type=ogr.wkbMultiLineString)
            elif iGeomType == ogr.wkbMultiPolygon:
                pLayer_out = pDataset_out.CreateLayer('layer', srs=pSpatialRef, geom_type=ogr.wkbMultiPolygon)
            else:
                print(f'Unsupported geometry type: {iGeomType}')
                sGeomType = ogr.GeometryTypeToName(iGeomType)
                print('Geometry type not supported:', sGeomType)
                continue

            # Copy field definitions from the first layer (if copy_attributes is True)
            if copy_attributes:
                pLayerDefn_in = pLayer_in.GetLayerDefn()
                for i in range(pLayerDefn_in.GetFieldCount()):
                    pFieldDefn = pLayerDefn_in.GetFieldDefn(i)
                    if pLayer_out.CreateField(pFieldDefn) != 0:
                        print(f'Failed to create field {pFieldDefn.GetName()} in output layer')
                        return

            # Add an id field if requested and it doesn't exist
            if add_id_field:
                pLayerDefn_in = pLayer_in.GetLayerDefn()
                if not copy_attributes or pLayerDefn_in.GetFieldIndex('id') == -1:
                    pField = ogr.FieldDefn('id', ogr.OFTInteger)
                    if pLayer_out.CreateField(pField) != 0:
                        print('Failed to create id field in output layer')
                        return

            iFlag_first = 0
        else:
            pass

        # Copy features from the input layer to the output layer
        nFeatures = pLayer_in.GetFeatureCount()
        nProcessed = 0
        nSkipped = 0

        if verbose:
            print(f'Processing {nFeatures} features from {sFilename_in}')

        for feature in pLayer_in:
            nProcessed += 1
            if verbose and nProcessed % 1000 == 0:
                print(f'  Processed {nProcessed}/{nFeatures} features')

            # Copy the feature
            pFeature_new = ogr.Feature(pLayer_out.GetLayerDefn())

            # Get the geometry
            original_geom = feature.GetGeometryRef()
            if original_geom is not None:
                # Set geometry directly without validation
                pFeature_new.SetGeometry(original_geom)
            else:
                if verbose:
                    print(f'Warning: Feature {nProcessed} has no geometry, skipping')
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
                pFeature_new.SetField('id', lid)
            pLayer_out.CreateFeature(pFeature_new)
            lid = lid + 1

        # Clean up
        pDataset_in = None

        # Print summary for this file
        if verbose:
            print(f'Completed {sFilename_in}: {nProcessed} processed, {nSkipped} skipped')

    # Clean up
    pDataset_out = None
    if verbose:
        print('=== Merge operation completed ===')
        print(f'Total features processed: {lid}')
        print('Merge completed successfully.')

    return

def main():
    """
    Main function for command line usage
    """
    import argparse

    parser = argparse.ArgumentParser(description='Merge multiple vector files into a single output file')
    parser.add_argument('inputs', nargs='+', help='Input vector files to merge')
    parser.add_argument('output', help='Output vector file')
    parser.add_argument('--format', help='Output format (auto-detected from extension if not specified)')
    parser.add_argument('--no-attributes', action='store_true',
                       help='Do not copy attributes from input files (geometries only)')
    parser.add_argument('--no-id', action='store_true',
                       help='Do not add sequential ID field to output')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress output')

    args = parser.parse_args()

    # Convert input arguments
    copy_attributes = not args.no_attributes
    add_id_field = not args.no_id
    verbose = not args.quiet

    if verbose:
        print(f'Merging {len(args.inputs)} files into {args.output}')
        print(f'Copy attributes: {copy_attributes}')
        print(f'Add ID field: {add_id_field}')

    merge_files(
        args.inputs,
        args.output,
        sFormat=args.format,
        copy_attributes=copy_attributes,
        add_id_field=add_id_field,
        verbose=verbose
    )

if __name__ == '__main__':
    main()

# Example usage:
# merge_files(['file1.shp', 'file2.shp'], 'merged.shp')  # Auto-detect shapefile format, copy all attributes
# merge_files(['file1.geojson', 'file2.geojson'], 'merged.gpkg', 'GPKG', copy_attributes=False)  # No attributes
# merge_files(['file1.shp', 'file2.geojson'], 'merged.geojson', add_id_field=False)  # No ID field