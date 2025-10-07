
import os, sys
import numpy as np
from osgeo import ogr, gdal, osr
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area  import calculate_polygon_area

from gcsbuffer.classes.vertex import pyvertex
from gcsbuffer.classes.edge import pyedge
from gcsbuffer.classes.polyline import pypolyline
from gcsbuffer.classes.polygon import pypolygon



def create_point_buffer_zone(sWkt, dBuffer_distance_in):
    #create a point geometry from WKT
    pGeometry = ogr.CreateGeometryFromWkt(sWkt)
    if pGeometry is None:
        print('Error: Invalid WKT input for point buffer creation.')
        return None
    sGeometry_type = pGeometry.GetGeometryName()
    if sGeometry_type != 'POINT':
        print('Error: Input geometry must be a POINT for buffer creation.')
        return None

    aCoords_gcs = get_geometry_coordinates(pGeometry)
    if len(aCoords_gcs) != 1:
        print('Error: Input geometry must be a single point for buffer creation.')
        return None

    point = dict()
    point['dLongitude_degree'] = aCoords_gcs[0][0]
    point['dLatitude_degree'] = aCoords_gcs[0][1]
    pVertex = pyvertex(point)

    sWkt_buffer_polygon = pVertex.calculate_buffer_zone(dBuffer_distance_in)

    return sWkt_buffer_polygon

def create_polyline_buffer_zone(sWkt, dBuffer_distance_in):

    #create a polyline geometry from WKT
    pGeometry = ogr.CreateGeometryFromWkt(sWkt)
    if pGeometry is None:
        print('Error: Invalid WKT input for polyline buffer creation.')
        return None
    sGeometry_type = pGeometry.GetGeometryName()
    if sGeometry_type == 'LINESTRING':
        aEdge = list()
        aCoords_gcs = get_geometry_coordinates(pGeometry)
        nPoint = len(aCoords_gcs) #remove the last point which is the same as the first point
        aVertex = list()
        point= dict()
        for i in range(0, nPoint):
            point['dLongitude_degree'] = aCoords_gcs[i][0]
            point['dLatitude_degree'] =  aCoords_gcs[i][1]
            pVertex = pyvertex(point)
            aVertex.append(pVertex)

        nVertex = len(aVertex)
        for i in range(0, nVertex-1):
            if aVertex[i] != aVertex[i+1]:
                pEdge = pyedge(aVertex[i], aVertex[i+1])
                aEdge.append(pEdge)
            else:
                pass

        ppolyline = pypolyline(aEdge)
        sWkt_buffer_polygon = ppolyline.calculate_buffer_zone(dBuffer_distance_in)
    else:
        print('Error: Input geometry must be a LINESTRING for buffer creation.')
        return None


    return sWkt_buffer_polygon

def create_buffer_zone_polygon_file(sFilename_polygon_in, sFilename_polygon_out,
                                    dThreshold_in=1.0E9,  # m2 to filter out small polygons
                                    dBuffer_distance_in=5000,  # m
                                    verbose=True):
    """
    Create buffer zones for polygons in a vector file

    :param sFilename_polygon_in: input vector file
    :param sFilename_polygon_out: output vector file
    :param dThreshold_in: area threshold in m² to filter out small polygons
    :param dBuffer_distance_in: buffer distance in meters
    :param verbose: whether to print progress information
    :return: None
    """

    if not os.path.exists(sFilename_polygon_in):
        print(f'Error: Input file {sFilename_polygon_in} does not exist!')
        return

    # Auto-detect formats from file extensions
    input_ext = os.path.splitext(sFilename_polygon_in)[1].lower()
    output_ext = os.path.splitext(sFilename_polygon_out)[1].lower()

    driver_map = get_supported_formats()
    output_driver_name = driver_map.get(output_ext, 'GeoJSON')

    if verbose:
        print(f'Input format: {input_ext}')
        print(f'Output format: {output_ext} -> {output_driver_name}')
        print(f'Buffer distance: {dBuffer_distance_in} meters')
        if dThreshold_in:
            print(f'Area threshold: {dThreshold_in} m²')

    # Open input dataset (auto-detect format)
    pDataSource = ogr.Open(sFilename_polygon_in, 0)
    if pDataSource is None:
        print(f'Error: Could not open input file {sFilename_polygon_in}')
        return

    pLayer = pDataSource.GetLayer()
    if pLayer is None:
        print('Error: No layer found in input file')
        pDataSource = None
        return

    # Get total feature count for progress tracking
    nTotal_features = pLayer.GetFeatureCount()
    if verbose:
        print(f'Processing {nTotal_features} features...')

    # Get output driver
    pDriver_out = ogr.GetDriverByName(output_driver_name)
    if pDriver_out is None:
        print(f'Error: Driver {output_driver_name} not available!')
        if verbose:
            print_supported_formats()
        pDataSource = None
        return

    # Prepare output (overwrite if exists)
    if os.path.exists(sFilename_polygon_out):
        if verbose:
            print(f'Removing existing output file: {sFilename_polygon_out}')
        pDriver_out.DeleteDataSource(sFilename_polygon_out)

    pOutDataSource = pDriver_out.CreateDataSource(sFilename_polygon_out)
    if pOutDataSource is None:
        print(f'Error: Could not create output file {sFilename_polygon_out}')
        pDataSource = None
        return

    # Handle layer creation based on output format
    if output_driver_name == 'ESRI Shapefile':
        # Shapefile requires simpler layer name
        layer_name = os.path.splitext(os.path.basename(sFilename_polygon_out))[0]
        pOutLayer = pOutDataSource.CreateLayer(layer_name, geom_type=ogr.wkbPolygon)
    else:
        pOutLayer = pOutDataSource.CreateLayer("buffer", geom_type=ogr.wkbPolygon)

    if pOutLayer is None:
        print('Error: Could not create output layer')
        pDataSource = None
        pOutDataSource = None
        return

    # Add fields for statistics
    pFieldDefn = ogr.FieldDefn('id', ogr.OFTInteger)
    pOutLayer.CreateField(pFieldDefn)
    pFieldDefn = ogr.FieldDefn('orig_area', ogr.OFTReal)
    pOutLayer.CreateField(pFieldDefn)
    pFieldDefn = ogr.FieldDefn('buffer_dist', ogr.OFTReal)
    pOutLayer.CreateField(pFieldDefn)

    pFeature = pLayer.GetNextFeature()
    nProcessed = 0
    nBuffered = 0
    feature_id = 1

    while pFeature:
        nProcessed += 1
        if verbose and nProcessed % 100 == 0:
            print(f'Processed {nProcessed}/{nTotal_features} features, created {nBuffered} buffers')

        pGeometry = pFeature.GetGeometryRef()
        if pGeometry is None:
            pFeature = pLayer.GetNextFeature()
            continue

        sGeometry_type = pGeometry.GetGeometryName()
        if sGeometry_type == 'MULTIPOLYGON':
            aaCoords_gcs = get_geometry_coordinates(pGeometry)
            nPart = len(aaCoords_gcs)
        elif sGeometry_type == 'POLYGON':
            aCoords_gcs = get_geometry_coordinates(pGeometry)
            aaCoords_gcs = [aCoords_gcs]
            nPart = 1
        else:
            if verbose:
                print(f'Warning: Skipping non-polygon geometry type {sGeometry_type} in feature {nProcessed}')
            pFeature = pLayer.GetNextFeature()
            continue

        for i in range(nPart):
            aCoords_gcs = aaCoords_gcs[i]
            if len(aCoords_gcs) < 3:  # Skip degenerate polygons
                continue

            aCoords_gcs = np.array(aCoords_gcs)
            # Calculate the area of the polygon
            dArea = calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1], iFlag_algorithm=2)

            if dThreshold_in is not None and dArea < dThreshold_in:
                continue

            nPoint = len(aCoords_gcs) - 1  # remove the last point which is the same as the first point
            aVertex = list()

            for j in range(nPoint):
                point = dict()
                point['dLongitude_degree'] = aCoords_gcs[j, 0]
                point['dLatitude_degree'] = aCoords_gcs[j, 1]
                pVertex = pyvertex(point)
                aVertex.append(pVertex)

            nVertex = len(aVertex)
            aEdge = list()
            for j in range(nVertex - 1):
                if aVertex[j] != aVertex[j + 1]:
                    pEdge = pyedge(aVertex[j], aVertex[j + 1])
                    aEdge.append(pEdge)

            # Close the polygon
            if nVertex > 0:
                pEdge = pyedge(aVertex[nVertex - 1], aVertex[0])
                aEdge.append(pEdge)

                try:
                    pPolygon = pypolygon(aEdge)
                    sWkt_buffer = pPolygon.calculate_buffer_zone(dBuffer_distance_in)

                    if sWkt_buffer:
                        buffer_geom = ogr.CreateGeometryFromWkt(sWkt_buffer)
                        if buffer_geom:
                            outFeature = ogr.Feature(pOutLayer.GetLayerDefn())
                            outFeature.SetGeometry(buffer_geom)
                            outFeature.SetField('id', feature_id)
                            outFeature.SetField('orig_area', dArea)
                            outFeature.SetField('buffer_dist', dBuffer_distance_in)
                            pOutLayer.CreateFeature(outFeature)
                            outFeature = None  # Free feature
                            nBuffered += 1
                            feature_id += 1
                except Exception as e:
                    if verbose:
                        print(f'Warning: Failed to create buffer for polygon part {i} in feature {nProcessed}: {str(e)}')

        pFeature = pLayer.GetNextFeature()

    # Cleanup
    pDataSource = None
    pOutDataSource = None

    if verbose:
        print(f'Processing complete!')
        print(f'Total features processed: {nProcessed}')
        print(f'Buffer zones created: {nBuffered}')
        print(f'Output saved to: {sFilename_polygon_out}')

    return

def main():
    """
    Main function for command line usage
    """
    import argparse

    parser = argparse.ArgumentParser(description='Create buffer zones for polygons in vector files')
    parser.add_argument('input', help='Input vector file')
    parser.add_argument('output', help='Output vector file')
    parser.add_argument('--buffer', type=float, default=5000, help='Buffer distance in meters (default: 5000)')
    parser.add_argument('--threshold', type=float, help='Area threshold in m² to filter small polygons')
    parser.add_argument('--quiet', action='store_true', help='Suppress progress output')
    parser.add_argument('--formats', action='store_true', help='Show supported formats and exit')

    args = parser.parse_args()

    if args.formats:
        print_supported_formats()
        return

    create_buffer_zone_polygon_file(
        args.input,
        args.output,
        dThreshold_in=args.threshold,
        dBuffer_distance_in=args.buffer,
        verbose=not args.quiet
    )

if __name__ == '__main__':
    main()



