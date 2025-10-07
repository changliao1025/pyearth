import os, sys
import numpy as np
from osgeo import ogr, osr, gdal
gdal.UseExceptions()
from pyearth.system.define_global_variables import *
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.douglas_peucker_geodetic import douglas_peucker_geodetic
from pyearth.gis.geometry.visvalingam_whyatt_geodetic import  visvalingam_whyatt_geodetic

from pyearth.gis.gdal.read.vector.get_supported_formats import get_supported_formats, print_supported_formats

def remove_small_polygon(sFilename_vector_in, sFilename_vector_out, dThreshold_in,
                        verbose=True, progress_interval=1000):
    """
    This function is used to remove small polygons from a vector file
    :param sFilename_vector_in: input vector file
    :param sFilename_vector_out: output vector file
    :param dThreshold_in: threshold for removing small polygons, in square kilo meters
    :param verbose: whether to print progress information
    :param progress_interval: how often to print progress updates
    :return: None
    """

    if not os.path.exists(sFilename_vector_in):
        print('Error: file %s does not exist!' % sFilename_vector_in)
        return

    if os.path.exists(sFilename_vector_out):
        if verbose:
            print(f'Removing existing output file: {sFilename_vector_out}')
        os.remove(sFilename_vector_out)

    dThreshold = float(dThreshold_in)

    # Automatically detect input format and set appropriate output driver
    input_ext = os.path.splitext(sFilename_vector_in)[1].lower()
    output_ext = os.path.splitext(sFilename_vector_out)[1].lower()

    # Map file extensions to OGR driver names
    driver_map = {
        '.geojson': 'GeoJSON',
        '.json': 'GeoJSON',
        '.shp': 'ESRI Shapefile',
        '.gpkg': 'GPKG',
        '.kml': 'KML',
        '.gml': 'GML',
        '.csv': 'CSV'
    }

    # Get input driver (for reading, we can use any driver that supports the format)
    input_driver_name = driver_map.get(input_ext, 'GeoJSON')  # default to GeoJSON
    output_driver_name = driver_map.get(output_ext, 'GeoJSON')  # default to GeoJSON

    if verbose:
        print(f'Input format: {input_ext} -> {input_driver_name}')
        print(f'Output format: {output_ext} -> {output_driver_name}')

    pDriver_out = ogr.GetDriverByName(output_driver_name)
    if pDriver_out is None:
        print(f'Error: Driver {output_driver_name} not available!')
        if verbose:
            print_supported_formats()
        return


    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)

    # Open input datasource (GDAL can auto-detect format)
    pDataSource_in = ogr.Open(sFilename_vector_in, 0)
    if pDataSource_in is None:
        print(f'Error: Could not open input file {sFilename_vector_in}')
        return

    pLayer_in = pDataSource_in.GetLayer()

    # Get total feature count for progress tracking
    nTotal_features = pLayer_in.GetFeatureCount()
    if verbose:
        print(f'Processing {nTotal_features} features with area threshold > {dThreshold} km²...')

    pDataSource_out = pDriver_out.CreateDataSource(sFilename_vector_out)
    if pDataSource_out is None:
        print(f'Error: Could not create output file {sFilename_vector_out}')
        pDataSource_in = None
        return

    # Handle layer creation based on output format
    if output_driver_name == 'ESRI Shapefile':
        # Shapefile requires simpler layer name
        layer_name = os.path.splitext(os.path.basename(sFilename_vector_out))[0]
        pLayer_out = pDataSource_out.CreateLayer(layer_name, pSrs, ogr.wkbPolygon)
    else:
        pLayer_out = pDataSource_out.CreateLayer('layer', pSrs, ogr.wkbPolygon)
    #create fields for id and area
    pFieldDefn = ogr.FieldDefn('id', ogr.OFTInteger)
    pLayer_out.CreateField(pFieldDefn)
    pFieldDefn = ogr.FieldDefn('area', ogr.OFTReal)
    pLayer_out.CreateField(pFieldDefn)
    pLayerDefn_in = pLayer_in.GetLayerDefn()
    nFieldCount = pLayerDefn_in.GetFieldCount()
    for i in range(nFieldCount):
        pFieldDefn = pLayerDefn_in.GetFieldDefn(i)
        if pFieldDefn.GetName() == 'id':
            continue
        if pFieldDefn.GetName() == 'area':
            continue
        pLayer_out.CreateField(pFieldDefn)

    pLayerDefn_out = pLayer_out.GetLayerDefn()

    # Pre-allocate arrays and optimize memory usage
    lID = 1
    nProcessed = 0
    nKept = 0

    # Use batch processing for better performance
    pLayer_in.SetNextByIndex(0)  # Reset to beginning

    for pFeature_in in pLayer_in:  # More pythonic iteration
        nProcessed += 1
        if verbose and nProcessed % progress_interval == 0:  # Progress indicator
            print(f'Processed {nProcessed}/{nTotal_features} features, kept {nKept}')

        pGeometry = pFeature_in.GetGeometryRef()
        if pGeometry is None:
            continue

        sGeometry_type = pGeometry.GetGeometryName()

        # Handle different geometry types more efficiently
        if sGeometry_type == 'POLYGON':
            # Get outer ring coordinates
            pOuterRing = pGeometry.GetGeometryRef(0)
            aCoords_outer = []
            for iPoint in range(pOuterRing.GetPointCount()):
                x, y, z = pOuterRing.GetPoint(iPoint)
                aCoords_outer.append([x, y])

            if len(aCoords_outer) < 3:  # Skip degenerate polygons
                continue

            aCoords_outer = np.array(aCoords_outer)

            # Calculate area of outer ring
            dArea = calculate_polygon_area(aCoords_outer[:,0], aCoords_outer[:,1], iFlag_algorithm=2)
            dAreakm = dArea * 1.0E-6

            if dAreakm > dThreshold:
                # Create new polygon with outer ring and all inner rings
                pGeometry_partial = ogr.Geometry(ogr.wkbPolygon)

                # Add outer ring
                pRing_outer = ogr.Geometry(ogr.wkbLinearRing)
                for j in range(len(aCoords_outer)):
                    pRing_outer.AddPoint(aCoords_outer[j, 0], aCoords_outer[j, 1])
                pRing_outer.CloseRings()
                pGeometry_partial.AddGeometry(pRing_outer)

                # Add inner rings (holes)
                nRings = pGeometry.GetGeometryCount()
                for iRing in range(1, nRings):  # Start from 1 to skip outer ring
                    pInnerRing_src = pGeometry.GetGeometryRef(iRing)
                    pInnerRing_new = ogr.Geometry(ogr.wkbLinearRing)
                    for iPoint in range(pInnerRing_src.GetPointCount()):
                        x, y, z = pInnerRing_src.GetPoint(iPoint)
                        pInnerRing_new.AddPoint(x, y)
                    pInnerRing_new.CloseRings()
                    pGeometry_partial.AddGeometry(pInnerRing_new)

                pGeometry_partial.AssignSpatialReference(pSrs)

                # Create output feature
                pFeature_out = ogr.Feature(pLayerDefn_out)
                pFeature_out.SetGeometry(pGeometry_partial)
                pFeature_out.SetField('id', lID)
                pFeature_out.SetField('area', dAreakm)

                # Copy other fields efficiently
                for k in range(nFieldCount):
                    field_name = pLayerDefn_in.GetFieldDefn(k).GetName()
                    if field_name not in ['id', 'area']:
                        pFeature_out.SetField(field_name, pFeature_in.GetField(field_name))

                pLayer_out.CreateFeature(pFeature_out)
                pFeature_out = None
                lID += 1
                nKept += 1

        elif sGeometry_type == 'MULTIPOLYGON':
            # For multipolygons, process each polygon part
            nPolygons = pGeometry.GetGeometryCount()
            for iPoly in range(nPolygons):
                pPolygon = pGeometry.GetGeometryRef(iPoly)

                # Get outer ring coordinates
                pOuterRing = pPolygon.GetGeometryRef(0)
                aCoords_outer = []
                for iPoint in range(pOuterRing.GetPointCount()):
                    x, y, z = pOuterRing.GetPoint(iPoint)
                    aCoords_outer.append([x, y])

                aCoords_outer = np.array(aCoords_outer)

                # Calculate area of outer ring
                dArea = calculate_polygon_area(aCoords_outer[:,0], aCoords_outer[:,1], iFlag_algorithm=2)
                dAreakm = dArea * 1.0E-6

                if dAreakm > dThreshold:
                    # Create new polygon with outer ring and all inner rings
                    pGeometry_partial = ogr.Geometry(ogr.wkbPolygon)

                    # Add outer ring
                    pRing_outer = ogr.Geometry(ogr.wkbLinearRing)
                    for j in range(len(aCoords_outer)):
                        pRing_outer.AddPoint(aCoords_outer[j, 0], aCoords_outer[j, 1])
                    pRing_outer.CloseRings()
                    pGeometry_partial.AddGeometry(pRing_outer)

                    # Add inner rings (holes)
                    nRings = pPolygon.GetGeometryCount()
                    for iRing in range(1, nRings):  # Start from 1 to skip outer ring
                        pInnerRing_src = pPolygon.GetGeometryRef(iRing)
                        pInnerRing_new = ogr.Geometry(ogr.wkbLinearRing)
                        for iPoint in range(pInnerRing_src.GetPointCount()):
                            x, y, z = pInnerRing_src.GetPoint(iPoint)
                            pInnerRing_new.AddPoint(x, y)
                        pInnerRing_new.CloseRings()
                        pGeometry_partial.AddGeometry(pInnerRing_new)

                    pGeometry_partial.AssignSpatialReference(pSrs)

                    # Create output feature
                    pFeature_out = ogr.Feature(pLayerDefn_out)
                    pFeature_out.SetGeometry(pGeometry_partial)
                    pFeature_out.SetField('id', lID)
                    pFeature_out.SetField('area', dAreakm)

                    # Copy other fields efficiently
                    for k in range(nFieldCount):
                        field_name = pLayerDefn_in.GetFieldDefn(k).GetName()
                        if field_name not in ['id', 'area']:
                            pFeature_out.SetField(field_name, pFeature_in.GetField(field_name))

                    pLayer_out.CreateFeature(pFeature_out)
                    pFeature_out = None
                    lID += 1
                    nKept += 1
        else:
            continue  # Skip non-polygon geometries

    # Cleanup and final statistics
    pDataSource_in = None  # Proper cleanup instead of Destroy()
    pDataSource_out = None

    if verbose:
        print(f'Processing complete!')
        print(f'Total features processed: {nProcessed}')
        print(f'Features kept (area > {dThreshold} km²): {nKept}')
        print(f'Features removed: {nProcessed - nKept}')
        print(f'Output saved to: {sFilename_vector_out}')
    return

def main():
    """
    Main function for command line usage
    """
    import argparse

    parser = argparse.ArgumentParser(description='Remove small polygons from vector files')
    parser.add_argument('input', help='Input vector file')
    parser.add_argument('output', help='Output vector file')
    parser.add_argument('threshold', type=float, help='Area threshold in square kilometers')
    parser.add_argument('--quiet', action='store_true', help='Suppress progress output')
    parser.add_argument('--progress', type=int, default=1000, help='Progress reporting interval')
    parser.add_argument('--formats', action='store_true', help='Show supported formats and exit')

    args = parser.parse_args()

    if args.formats:
        print_supported_formats()
        return

    remove_small_polygon(
        args.input,
        args.output,
        args.threshold,
        verbose=not args.quiet,
        progress_interval=args.progress
    )

if __name__ == '__main__':
    main()