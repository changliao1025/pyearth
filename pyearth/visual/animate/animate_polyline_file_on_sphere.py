import os
import logging
import traceback
import numpy as np
import math

import signal
import sys
import time
from osgeo import gdal, ogr, osr
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor, as_completed
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import extract_unique_vertices_and_connectivity
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
# Set up logger
logger = logging.getLogger(__name__)


def animate_polyline_file_on_sphere(
    sFilename_polyline_in,
    sFilename_animation_out=None,
    dLongitude_focus_in=0.0,
    dLatitude_focus_in=0.0,
    dZoom_factor=0.7,
    iFlag_show_coastlines=True,
    iFlag_show_graticule=True,
    sTitle_in=None,
    sColor_polyline_in='red',
    dLinewidth_in=2.0,
    iFigwidth_in=None,
    iFigheight_in=None
):
    """
    Visualize polyline data from a file on a 3D sphere using GeoVista.

    Creates an interactive or saved 3D visualization of polyline geometries
    mapped on a spherical globe with proper geographic context including
    coastlines and coordinate grid.

    Args:
        sFilename_polyline_in (str): Input polyline file path.
            Supports formats: .shp, .geojson, .kml, .gpx
        sFilename_animation_out (str, optional): Output screenshot file path.
            If None, displays interactive viewer. Supports: .png, .jpg, .svg
        dLongitude_focus_in (float, optional): Camera focal point longitude in degrees.
            Valid range: -180 to 180. Default is 0.0 (prime meridian).
        dLatitude_focus_in (float, optional): Camera focal point latitude in degrees.
            Valid range: -90 to 90. Default is 0.0 (equator).
        dZoom_factor (float, optional): Camera zoom level.
            Higher values zoom in. Default is 0.7.
        iFlag_show_coastlines (bool, optional): Show coastline overlay.
            Default is True.
        iFlag_show_graticule (bool, optional): Show coordinate grid with labels.
            Default is True.
        sTitle_in (str, optional): Title for the visualization.
            If None, uses filename.
        sColor_polyline_in (str, optional): Color for polylines.
            Default is 'red'. Supports named colors or hex codes.
        dLinewidth_in (float, optional): Line width for polylines.
            Default is 2.0.
        iFigwidth_in (int, optional): Figure width (for compatibility).
        iFigheight_in (int, optional): Figure height (for compatibility).

    Returns:
        bool: True if visualization successful, False otherwise

    Note:
        - Requires 'geovista' and 'pyvista' packages: pip install geovista pyvista
        - Requires 'gdal' for file reading: pip install gdal
        - Interactive mode requires display environment
    """

    # Validate input file
    if not isinstance(sFilename_polyline_in, str) or not sFilename_polyline_in.strip():
        logger.error('Input polyline filename must be a non-empty string')
        return False

    if not os.path.exists(sFilename_polyline_in):
        logger.error(f'Input polyline file not found: {sFilename_polyline_in}')
        return False



    # Validate focus coordinates
    dLongitude_focus = dLongitude_focus_in if dLongitude_focus_in is not None else 0.0
    dLatitude_focus = dLatitude_focus_in if dLatitude_focus_in is not None else 0.0

    if not (-180 <= dLongitude_focus <= 180):
        logger.warning(f'Longitude focus {dLongitude_focus} out of range [-180, 180], clamping')
        dLongitude_focus = np.clip(dLongitude_focus, -180, 180)

    if not (-90 <= dLatitude_focus <= 90):
        logger.warning(f'Latitude focus {dLatitude_focus} out of range [-90, 90], clamping')
        dLatitude_focus = np.clip(dLatitude_focus, -90, 90)

    # Validate zoom factor
    if dZoom_factor <= 0:
        logger.warning(f'Invalid zoom factor {dZoom_factor}, using default 0.7')
        dZoom_factor = 0.7

    # Validate line width
    if dLinewidth_in <= 0:
        logger.warning(f'Invalid line width {dLinewidth_in}, using default 2.0')
        dLinewidth_in = 2.0

    # Validate output file path if provided
    if sFilename_animation_out is not None:
        if not isinstance(sFilename_animation_out, str) or not sFilename_animation_out.strip():
            logger.error('Output filename must be a non-empty string')
            return False

        # Check output directory exists
        output_dir = os.path.dirname(sFilename_animation_out)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
                logger.info(f'Created output directory: {output_dir}')
            except Exception as e:
                logger.error(f'Cannot create output directory {output_dir}: {e}')
                return False

        # Check supported file extensions
        valid_output_extensions = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']
        output_ext = os.path.splitext(sFilename_animation_out)[1].lower()
        if output_ext not in valid_output_extensions:
            logger.warning(f'Output extension {output_ext} may not be supported. Recommended: .png, .jpg, .svg')


    try:
        import geovista as gv
        import pyvista as pv
        logger.info('GeoVista and PyVista libraries imported successfully')
    except ImportError as e:
        logger.error('GeoVista/PyVista libraries not available. Install with: pip install geovista pyvista')
        logger.error(f'Import error: {e}')
        return False

    try:
        logger.info('Loading polyline data...')
        logger.info(f'  - Input file: {sFilename_polyline_in}')
        logger.info(f'  - Focus: ({dLongitude_focus:.2f}°, {dLatitude_focus:.2f}°)')
        logger.info(f'  - Zoom factor: {dZoom_factor}')
        logger.info(f'  - Line color: {sColor_polyline_in}')
        logger.info(f'  - Line width: {dLinewidth_in}')

        # Load polyline data using GDAL/OGR
        try:
            # Open the vector dataset
            pDataset = ogr.Open(sFilename_polyline_in, gdal.GA_ReadOnly)
            if pDataset is None:
                logger.error(f'Could not open vector file: {sFilename_polyline_in}')
                return False

            # Get the first layer (assuming single layer)
            pLayer = pDataset.GetLayer(0)
            if pLayer is None:
                logger.error('Could not get layer from dataset')
                return False

            lFeatureCount = pLayer.GetFeatureCount()
            logger.info(f'Loaded {lFeatureCount} features from vector file')

            if lFeatureCount == 0:
                logger.error('No features found in the input file')
                return False

            # Get spatial reference
            pSpatialRef = pLayer.GetSpatialRef()
            logger.debug(f'Input spatial reference: {pSpatialRef}')

            # Check if we need coordinate transformation to WGS84
            wgs84_srs = osr.SpatialReference()
            wgs84_srs.ImportFromEPSG(4326)

            pTransform = None
            if pSpatialRef is not None and not pSpatialRef.IsSame(wgs84_srs):
                logger.info('Setting up coordinate transformation to WGS84')
                pTransform = osr.CoordinateTransformation(pSpatialRef, wgs84_srs)

        except Exception as e:
            logger.error(f'Failed to read polyline file: {e}')
            return False

        # Collect valid line geometries
        aLineGeometries = []

        # Reset layer reading to start from the beginning
        pLayer.ResetReading()

        for pFeature in pLayer:
            pGeometry = pFeature.GetGeometryRef()
            if pGeometry is None:
                continue

            geom_type = pGeometry.GetGeometryType()
            geom_name = pGeometry.GetGeometryName()

            # Check for LineString or MultiLineString
            if geom_type == ogr.wkbLineString or geom_type == ogr.wkbMultiLineString:
                # Clone geometry to avoid issues with feature lifecycle
                pGeomClone = pGeometry.Clone()

                # Transform to WGS84 if needed
                if pTransform is not None:
                    pGeomClone.Transform(pTransform)

                aLineGeometries.append(pGeomClone)
            else:
                logger.debug(f'Skipping geometry type: {geom_name}')

        if len(aLineGeometries) == 0:
            logger.error('No valid LineString or MultiLineString geometries found')
            return False

        logger.info(f'Found {len(aLineGeometries)} valid line geometries')


        # Create 3D plotter
        plotter = gv.GeoPlotter()

        # Set title
        if sTitle_in is not None:
            title = sTitle_in
        else:
            title = f'Polylines: {os.path.basename(sFilename_polyline_in)}'

        plotter.add_title(title, font_size=14)

        # Add polylines to the plotter
        logger.info('Adding polylines to 3D sphere...')

        # Collect all polyline coordinates and connectivity for batch processing
        all_points = []
        all_lines = []
        current_point_offset = 0
        polyline_count = 0

        logger.info('Processing polylines for batch visualization...')

        for geom_idx, pGeometry in enumerate(aLineGeometries):
            try:
                # Helper function to process individual LineString
                def process_linestring(line_geometry, geom_name):
                    nonlocal current_point_offset, polyline_count

                    coords = []
                    for i in range(line_geometry.GetPointCount()):
                        x, y, z = line_geometry.GetPoint(i)
                        coords.append([x, y, z if z != 0 else 0])

                    if len(coords) >= 2:  # Valid line needs at least 2 points
                        # Add points to global array
                        line_points = np.array(coords)
                        all_points.extend(line_points)

                        # Create line connectivity with global point indices
                        n_points = len(line_points)
                        line_segments = np.column_stack((
                            np.full(n_points - 1, 2),  # Each line segment has 2 points
                            np.arange(current_point_offset, current_point_offset + n_points - 1),  # Start indices
                            np.arange(current_point_offset + 1, current_point_offset + n_points)   # End indices
                        )).flatten()

                        all_lines.extend(line_segments)
                        current_point_offset += n_points
                        polyline_count += 1

                        logger.debug(f'Processed {geom_name} {geom_idx} with {n_points} points')

                # Process different geometry types
                if pGeometry.GetGeometryType() == ogr.wkbLineString:
                    process_linestring(pGeometry, 'LineString')

                elif pGeometry.GetGeometryType() == ogr.wkbMultiLineString:
                    # MultiLineString - process each component
                    for geom_part_idx in range(pGeometry.GetGeometryCount()):
                        pLineString = pGeometry.GetGeometryRef(geom_part_idx)
                        process_linestring(pLineString, f'MultiLineString[{geom_part_idx}]')

            except Exception as e:
                logger.warning(f'Error processing geometry at index {geom_idx}: {e}')
                continue

        if polyline_count == 0:
            logger.error('No valid polylines could be processed for visualization')
            return False

        logger.info(f'Collected {polyline_count} polylines with {len(all_points)} total points')

        # Create combined PyVista polydata for all polylines (batch processing)
        try:
            # Convert to numpy arrays
            combined_points = np.array(all_points)
            combined_lines = np.array(all_lines)

            # Create single PyVista polydata with all lines
            logger.info('Creating combined polydata for batch processing...')
            combined_polydata = pv.PolyData(combined_points, lines=combined_lines)

            # Transform entire dataset to sphere at once (much faster)
            logger.info('Transforming all polylines to sphere coordinates...')
            sphere_polydata = gv.Transform.from_1d(combined_polydata, crs='EPSG:4326')

            # Add single combined mesh to plotter (batch operation)
            logger.info('Adding combined polyline mesh to visualization...')
            plotter.add_mesh(
                sphere_polydata,
                color=sColor_polyline_in,
                line_width=dLinewidth_in,
                render_lines_as_tubes=True
            )

            logger.info(f'✓ Successfully added {polyline_count} polylines as single combined mesh')

        except Exception as e:
            logger.error(f'Failed to create combined polyline mesh: {e}')
            logger.error(f'Traceback: {traceback.format_exc()}')
            return False



        # Configure camera position and focus
        try:
            # Convert longitude/latitude to radians
            lon_rad = math.radians(dLongitude_focus)
            lat_rad = math.radians(dLatitude_focus)

            # Earth radius (approximately 6371 km, but use normalized units)
            earth_radius = 1.0
            camera_distance = earth_radius * 3.0  # Position camera 3x earth radius away

            # Convert spherical coordinates to Cartesian (x, y, z)
            x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
            y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
            z_focal = earth_radius * math.sin(lat_rad)

            x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
            y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
            z_camera = camera_distance * math.sin(lat_rad)

            focal_point = [x_focal, y_focal, z_focal]
            camera_position = [x_camera, y_camera, z_camera]

            plotter.camera.focal_point = focal_point
            plotter.camera.position = camera_position
            plotter.camera.zoom(dZoom_factor)

            logger.debug(f'Camera configured: focal={focal_point}, position={camera_position}')
        except Exception as e:
            logger.warning(f'Error setting camera position: {e}. Using default view.')

        # Add geographic context
        if iFlag_show_coastlines:
            try:
                plotter.add_coastlines()
                logger.debug('Added coastlines overlay')
            except Exception as e:
                logger.warning(f'Could not add coastlines: {e}')

        # Add coordinate axes
        try:
            plotter.add_axes()
            logger.debug('Added coordinate axes')
        except Exception as e:
            logger.warning(f'Could not add axes: {e}')

        # Add graticule (coordinate grid)
        if iFlag_show_graticule:
            try:
                plotter.add_graticule(show_labels=True)
                logger.debug('Added coordinate graticule with labels')
            except Exception as e:
                logger.warning(f'Could not add graticule: {e}')

        # Output or display
        if sFilename_animation_out is not None:
            # Save screenshot
            try:
                plotter.screenshot(sFilename_animation_out)
                logger.info(f'✓ Visualization saved to: {sFilename_animation_out}')

                # Verify file was created
                if os.path.exists(sFilename_animation_out):
                    file_size = os.path.getsize(sFilename_animation_out)
                    logger.info(f'  File size: {file_size / 1024:.1f} KB')
                else:
                    logger.warning(f'Screenshot command executed but file not found: {sFilename_animation_out}')

                plotter.close()
                return True

            except Exception as e:
                logger.error(f'Failed to save screenshot: {e}')
                logger.error(f'Traceback: {traceback.format_exc()}')
                plotter.close()
                return False
        else:
            # Interactive display
            try:
                logger.info('Opening interactive visualization window...')
                plotter.show()
                return True
            except Exception as e:
                logger.error(f'Failed to display interactive visualization: {e}')
                logger.error(f'Ensure display environment is available (X11, Wayland, etc.)')
                logger.error(f'Traceback: {traceback.format_exc()}')
                plotter.close()
                return False

    except ImportError as e:
        logger.error(f'Missing required dependencies: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during polyline visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False
