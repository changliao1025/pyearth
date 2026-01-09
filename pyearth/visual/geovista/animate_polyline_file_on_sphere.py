import os
import logging
import traceback
import geovista
import numpy as np
import math

import signal
import sys
import time
from osgeo import gdal, ogr, osr
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor, as_completed
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.extract_unique_vertices_and_connectivity import (
    extract_unique_vertices_and_connectivity,
)
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from geovista.geodesic import line as gv_line

# Set up logger
logger = logging.getLogger(__name__)


def transform_vectors_sph_to_cart(theta, phi, r, u, v, w):
    """Transform vectors from spherical (r, phi, theta) to cartesian coordinates (z, y, x).

    Note the "reverse" order of arrays's axes, commonly used in geosciences.

    Parameters
    ----------
    theta : array_like[float]
        Azimuthal angle in degrees ``[0, 360]`` of shape (M,)
    phi : array_like[float]
        Polar (zenith) angle in degrees ``[0, 180]`` of shape (N,)
    r : array_like[float]
        Distance (radius) from the point of origin of shape (P,)
    u : array_like[float]
        X-component of the vector of shape (P, N, M)
    v : array_like[float]
        Y-component of the vector of shape (P, N, M)
    w : array_like[float]
        Z-component of the vector of shape (P, N, M)

    Returns
    -------
    u_t, v_t, w_t : :class:`numpy.ndarray`
        Arrays of transformed x-, y-, z-components, respectively.

    """
    xx, yy, _ = np.meshgrid(np.radians(theta), np.radians(phi), r, indexing="ij")
    th, ph = xx.squeeze(), yy.squeeze()

    # Transform wind components from spherical to cartesian coordinates
    # https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
    u_t = np.sin(ph) * np.cos(th) * w + np.cos(ph) * np.cos(th) * v - np.sin(th) * u
    v_t = np.sin(ph) * np.sin(th) * w + np.cos(ph) * np.sin(th) * v + np.cos(th) * u
    w_t = np.cos(ph) * w - np.sin(ph) * v

    return u_t, v_t, w_t


def animate_polyline_file_on_sphere(
    sFilename_polyline_in,
    sFilename_animation_out=None,
    dLongitude_focus_in=0.0,
    dLatitude_focus_in=0.0,
    dZoom_factor=0.7,
    iFlag_show_coastlines=True,
    iFlag_show_graticule=True,
    sTitle_in=None,
    sColor_polyline_in="royalblue",
    dLinewidth_in=2.0,
    sLinewidth_attribute_in=None,
    dLinewidth_min_in=0.5,
    dLinewidth_max_in=3.0,
    iFigwidth_in=None,
    iFigheight_in=None,
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
            Default is 'royalblue'. Supports named colors or hex codes.
        dLinewidth_in (float, optional): Default line width for polylines.
            Default is 2.0. Used when no variable width data provided.
        sLinewidth_attribute_in (str, optional): Attribute name from input file to use for line width scaling.
            Values read from file attributes and scaled to [dLinewidth_min_in, dLinewidth_max_in].
        dLinewidth_min_in (float, optional): Minimum line width for scaling.
            Default is 0.5.
        dLinewidth_max_in (float, optional): Maximum line width for scaling.
            Default is 3.0.
        iFigwidth_in (int, optional): Figure width (for compatibility).
        iFigheight_in (int, optional): Figure height (for compatibility).

    Returns:
        bool: True if visualization successful, False otherwise

    Note:
        - Requires 'geovista' and 'pyvista' packages: pip install geovista pyvista
        - Requires 'gdal' for file reading: pip install gdal
        - Interactive mode requires display environment
        - Variable line width requires individual mesh processing (slower than batch)
    """

    # Validate input file
    if not isinstance(sFilename_polyline_in, str) or not sFilename_polyline_in.strip():
        logger.error("Input polyline filename must be a non-empty string")
        return False

    if not os.path.exists(sFilename_polyline_in):
        logger.error(f"Input polyline file not found: {sFilename_polyline_in}")
        return False

    # Validate focus coordinates
    dLongitude_focus = dLongitude_focus_in if dLongitude_focus_in is not None else 0.0
    dLatitude_focus = dLatitude_focus_in if dLatitude_focus_in is not None else 0.0

    if not (-180 <= dLongitude_focus <= 180):
        logger.warning(
            f"Longitude focus {dLongitude_focus} out of range [-180, 180], clamping"
        )
        dLongitude_focus = np.clip(dLongitude_focus, -180, 180)

    if not (-90 <= dLatitude_focus <= 90):
        logger.warning(
            f"Latitude focus {dLatitude_focus} out of range [-90, 90], clamping"
        )
        dLatitude_focus = np.clip(dLatitude_focus, -90, 90)

    # Validate zoom factor
    if dZoom_factor <= 0:
        logger.warning(f"Invalid zoom factor {dZoom_factor}, using default 0.7")
        dZoom_factor = 0.7

    # Validate line width
    if dLinewidth_in <= 0:
        logger.warning(f"Invalid line width {dLinewidth_in}, using default 2.0")
        dLinewidth_in = 2.0

    # Validate output file path if provided
    if sFilename_animation_out is not None:
        if (
            not isinstance(sFilename_animation_out, str)
            or not sFilename_animation_out.strip()
        ):
            logger.error("Output filename must be a non-empty string")
            return False

        # Check output directory exists
        output_dir = os.path.dirname(sFilename_animation_out)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
                logger.info(f"Created output directory: {output_dir}")
            except Exception as e:
                logger.error(f"Cannot create output directory {output_dir}: {e}")
                return False

        # Check supported file extensions
        valid_output_extensions = [
            ".png",
            ".jpg",
            ".jpeg",
            ".svg",
            ".tif",
            ".tiff",
            "mp4",
        ]
        output_ext = os.path.splitext(sFilename_animation_out)[1].lower()
        if output_ext not in valid_output_extensions:
            logger.warning(
                f"Output extension {output_ext} may not be supported. Recommended: .png, .jpg, .svg"
            )

    try:
        import geovista as gv
        import pyvista as pv

        logger.info("GeoVista and PyVista libraries imported successfully")
    except ImportError as e:
        logger.error(
            "GeoVista/PyVista libraries not available. Install with: pip install geovista pyvista"
        )
        logger.error(f"Import error: {e}")
        return False

    try:
        logger.info("Loading polyline data...")
        logger.info(f"  - Input file: {sFilename_polyline_in}")
        logger.info(f"  - Focus: ({dLongitude_focus:.2f}°, {dLatitude_focus:.2f}°)")
        logger.info(f"  - Zoom factor: {dZoom_factor}")
        logger.info(f"  - Line color: {sColor_polyline_in}")
        logger.info(f"  - Line width: {dLinewidth_in}")

        # Load polyline data using GDAL/OGR
        try:
            # Open the vector dataset
            pDataset = ogr.Open(sFilename_polyline_in, gdal.GA_ReadOnly)
            if pDataset is None:
                logger.error(f"Could not open vector file: {sFilename_polyline_in}")
                return False

            # Get the first layer (assuming single layer)
            pLayer = pDataset.GetLayer(0)
            if pLayer is None:
                logger.error("Could not get layer from dataset")
                return False

            lFeatureCount = pLayer.GetFeatureCount()
            logger.info(f"Loaded {lFeatureCount} features from vector file")

            if lFeatureCount == 0:
                logger.error("No features found in the input file")
                return False

            # Get spatial reference
            pSpatialRef = pLayer.GetSpatialRef()
            logger.debug(f"Input spatial reference: {pSpatialRef}")
            pSpatialRef_wkt = (
                pSpatialRef.ExportToWkt() if pSpatialRef is not None else "None"
            )

            # Check if we need coordinate transformation to WGS84
            wgs84_srs = osr.SpatialReference()
            wgs84_srs.ImportFromEPSG(4326)
            wgs84_wkt = wgs84_srs.ExportToWkt()

            pTransform = None
            if pSpatialRef_wkt != wgs84_wkt:
                logger.info("Setting up coordinate transformation to WGS84")
                pTransform = osr.CoordinateTransformation(pSpatialRef, wgs84_srs)

        except Exception as e:
            logger.error(f"Failed to read polyline file: {e}")
            return False

        # Collect valid line geometries and attributes
        aLineGeometries = []
        aLinewidthValues = []

        # Reset layer reading to start from the beginning
        pLayer.ResetReading()

        for pFeature in pLayer:
            pGeometry = pFeature.GetGeometryRef()
            if pGeometry is None:
                continue

            geom_type = pGeometry.GetGeometryType()
            geom_name = pGeometry.GetGeometryName()

            # Check for LineString or MultiLineString
            if geom_name == "LINESTRING" or geom_name == "MULTILINESTRING":
                # Clone geometry to avoid issues with feature lifecycle
                pGeomClone = pGeometry.Clone()

                # Transform to WGS84 if needed
                if pTransform is not None:
                    pGeomClone.Transform(pTransform)

                aLineGeometries.append(pGeomClone)

                # Extract line width attribute if specified
                if sLinewidth_attribute_in is not None:
                    try:
                        width_value = pFeature.GetField(sLinewidth_attribute_in)
                        if width_value is not None:
                            aLinewidthValues.append(float(width_value))
                        else:
                            aLinewidthValues.append(dLinewidth_in)
                    except Exception as e:
                        logger.warning(
                            f"Could not read attribute {sLinewidth_attribute_in}: {e}"
                        )
                        aLinewidthValues.append(dLinewidth_in)
                else:
                    aLinewidthValues.append(dLinewidth_in)
            else:
                logger.debug(f"Skipping geometry type: {geom_name}")

        if len(aLineGeometries) == 0:
            logger.error("No valid LineString or MultiLineString geometries found")
            return False

        logger.info(f"Found {len(aLineGeometries)} valid line geometries")

        # Handle line width data
        use_variable_width = False
        scaled_widths = []

        if sLinewidth_attribute_in is not None and len(aLinewidthValues) > 0:
            # Use attribute data from file
            use_variable_width = True
            # Scale the attribute values to the specified range
            data_min, data_max = min(aLinewidthValues), max(aLinewidthValues)
            if data_max > data_min:
                for val in aLinewidthValues:
                    scaled_width = dLinewidth_min_in + (val - data_min) / (
                        data_max - data_min
                    ) * (dLinewidth_max_in - dLinewidth_min_in)
                    scaled_widths.append(scaled_width)
                logger.info(
                    f'Using variable line widths from attribute "{sLinewidth_attribute_in}" (range: {data_min:.2f} to {data_max:.2f})'
                )
                logger.info(
                    f"Scaled to width range: {dLinewidth_min_in} to {dLinewidth_max_in}"
                )
            else:
                scaled_widths = [dLinewidth_in] * len(aLineGeometries)
        else:
            # Use uniform line width
            scaled_widths = [dLinewidth_in] * len(aLineGeometries)

        # Create 3D plotter
        if sFilename_animation_out is not None:
            plotter = gv.GeoPlotter(off_screen=True)
        else:
            plotter = gv.GeoPlotter()

        # Set title
        if sTitle_in is not None:
            title = sTitle_in
        else:
            title = f"Polylines: {os.path.basename(sFilename_polyline_in)}"

        plotter.add_title(title, font_size=14)

        # Add polylines to the plotter using optimized batch processing
        logger.info("Adding polylines to 3D sphere using batch processing...")

        # Method 1: Use GeoVista's multi-line support (most efficient)
        try:
            # Collect all coordinates for batch processing with GeoVista
            all_lons_list = []
            all_lats_list = []
            polyline_count = 0

            logger.info(
                "Processing polylines for GeoVista multi-line batch visualization..."
            )

            for geom_idx, pGeometry in enumerate(aLineGeometries):
                try:
                    points = pGeometry.GetPoints()
                    if points is None or len(points) < 2:
                        continue

                    lons = [point[0] for point in points]
                    lats = [point[1] for point in points]

                    all_lons_list.append(lons)
                    all_lats_list.append(lats)
                    polyline_count += 1

                except Exception as e:
                    logger.warning(
                        f"Error processing geometry at index {geom_idx}: {e}"
                    )
                    continue

            if polyline_count == 0:
                logger.error("No valid polylines could be processed for visualization")
                return False

            # Create combined mesh using GeoVista's multi-line support
            combined_flow_field = gv_line(lons=all_lons_list, lats=all_lats_list)

            if use_variable_width:
                # For variable width, we need to add scalars to the mesh
                try:
                    # Add line width scalars to the mesh
                    import numpy as np

                    # Create point-based scalars for line widths
                    point_widths = []
                    for i, width in enumerate(scaled_widths):
                        # Each line gets repeated width values for all its points
                        num_points = len(all_lons_list[i])
                        point_widths.extend([width] * num_points)

                    combined_flow_field.point_data["line_width"] = np.array(
                        point_widths
                    )

                    plotter.add_mesh(
                        combined_flow_field,
                        color=sColor_polyline_in,
                        scalars="line_width",
                        line_width=None,  # Let scalars control width
                        name="geovista_multilines_variable",
                    )
                    logger.info(
                        f"✓ Successfully added {polyline_count} polylines using GeoVista multi-line with variable widths (1 operation)"
                    )

                except Exception as width_error:
                    logger.warning(
                        f"Variable width with combined mesh failed: {width_error}"
                    )
                    logger.info(
                        "Falling back to individual mesh processing for variable widths..."
                    )

                    # Fallback: Individual processing for variable line widths
                    polyline_count = 0
                    for geom_idx, (pGeometry, width) in enumerate(
                        zip(aLineGeometries, scaled_widths)
                    ):
                        try:
                            points = pGeometry.GetPoints()
                            if points is None or len(points) < 2:
                                continue

                            lons = [point[0] for point in points]
                            lats = [point[1] for point in points]

                            flow_field = gv_line(lons=lons, lats=lats)
                            _ = plotter.add_mesh(
                                flow_field,
                                color=sColor_polyline_in,
                                line_width=width,
                                name=f"polyline_{geom_idx}",
                            )
                            polyline_count += 1

                        except Exception as e:
                            logger.warning(
                                f"Error processing geometry at index {geom_idx}: {e}"
                            )
                            continue

                    if polyline_count == 0:
                        logger.error(
                            "No valid polylines could be processed for visualization"
                        )
                        return False

                    logger.info(
                        f"✓ Successfully added {polyline_count} polylines with variable widths using individual processing"
                    )
            else:
                # Uniform width
                plotter.add_mesh(
                    combined_flow_field,
                    color=sColor_polyline_in,
                    line_width=dLinewidth_in,
                    name="geovista_multilines",
                )
                logger.info(
                    f"✓ Successfully added {polyline_count} polylines using GeoVista multi-line (1 operation)"
                )

        except Exception as e:
            logger.error(
                f"GeoVista multi-line failed, trying PyVista batch method: {e}"
            )

            # Method 2: Combine all polylines into a single PyVista PolyData object
            try:
                import pyvista as pv

                # Collect all line coordinates and connectivity
                all_points = []
                all_lines = []
                point_offset = 0
                polyline_count = 0

                logger.info("Processing polylines for PyVista batch visualization...")

                for geom_idx, pGeometry in enumerate(aLineGeometries):
                    try:
                        # Get points from the geometry
                        points = pGeometry.GetPoints()
                        if points is None or len(points) < 2:
                            continue

                        # Extract longitude and latitude, convert to 3D coordinates on sphere
                        for point in points:
                            lon, lat = point[0], point[1]
                            # Convert to Cartesian coordinates on unit sphere
                            lon_rad = math.radians(lon)
                            lat_rad = math.radians(lat)
                            x = math.cos(lat_rad) * math.cos(lon_rad)
                            y = math.cos(lat_rad) * math.sin(lon_rad)
                            z = math.sin(lat_rad)
                            all_points.append([x, y, z])

                        # Create line connectivity
                        num_points_in_line = len(points)
                        if num_points_in_line >= 2:
                            line = [
                                num_points_in_line
                            ]  # First element is number of points
                            line.extend(
                                range(point_offset, point_offset + num_points_in_line)
                            )
                            all_lines.extend(line)
                            point_offset += num_points_in_line
                            polyline_count += 1

                    except Exception as e:
                        logger.warning(
                            f"Error processing geometry at index {geom_idx}: {e}"
                        )
                        continue

                if polyline_count == 0:
                    logger.error(
                        "No valid polylines could be processed for visualization"
                    )
                    return False

                # Create single PyVista PolyData object with all polylines
                polydata = pv.PolyData(all_points, lines=all_lines)

                # Add the combined mesh to plotter - single operation instead of loop
                plotter.add_mesh(
                    polydata,
                    color=sColor_polyline_in,
                    line_width=dLinewidth_in,
                    name="combined_polylines",
                )

                logger.info(
                    f"✓ Successfully added {polyline_count} polylines as single PyVista combined mesh"
                )
                logger.info(f"  Total points: {len(all_points)}")
                logger.info(
                    f"  Combined into 1 mesh operation (vs {polyline_count} separate operations)"
                )

            except Exception as e:
                logger.error(
                    f"PyVista batch method failed, falling back to individual processing: {e}"
                )
                logger.error(f"Traceback: {traceback.format_exc()}")

                # Method 3: Fallback - Use individual mesh additions (original approach)
                logger.info("Using fallback method: individual mesh additions...")
                polyline_count = 0

                for geom_idx, pGeometry in enumerate(aLineGeometries):
                    try:
                        points = pGeometry.GetPoints()
                        if points is None:
                            continue

                        lons = [point[0] for point in points]
                        lats = [point[1] for point in points]

                        flow_field = gv_line(lons=lons, lats=lats)
                        _ = plotter.add_mesh(
                            flow_field,
                            color=sColor_polyline_in,
                            line_width=dLinewidth_in,
                        )
                        polyline_count += 1

                    except Exception as e:
                        logger.warning(
                            f"Error processing geometry at index {geom_idx}: {e}"
                        )
                        continue

                if polyline_count == 0:
                    logger.error(
                        "No valid polylines could be processed for visualization"
                    )
                    return False

                logger.info(
                    f"✓ Completed using fallback method: {polyline_count} individual mesh operations"
                )

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

            logger.debug(
                f"Camera configured: focal={focal_point}, position={camera_position}"
            )
        except Exception as e:
            logger.warning(f"Error setting camera position: {e}. Using default view.")

        # Add geographic context
        if iFlag_show_coastlines:
            try:
                plotter.add_coastlines()
                logger.debug("Added coastlines overlay")
            except Exception as e:
                logger.warning(f"Could not add coastlines: {e}")

        # Add coordinate axes
        try:
            plotter.add_axes()
            logger.debug("Added coordinate axes")
        except Exception as e:
            logger.warning(f"Could not add axes: {e}")

        # Add graticule (coordinate grid)
        if iFlag_show_graticule:
            try:
                plotter.add_graticule(show_labels=True)
                logger.debug("Added coordinate graticule with labels")
            except Exception as e:
                logger.warning(f"Could not add graticule: {e}")

        # Output or display
        if sFilename_animation_out is not None:
            output_ext = os.path.splitext(sFilename_animation_out)[1].lower()
            if output_ext == ".mp4":
                try:
                    success = _create_rotation_animation(
                        plotter,
                        sFilename_animation_out,
                        dLongitude_focus,
                        dLatitude_focus,
                        360,
                        1,
                        "mp4",
                    )
                    plotter.close()
                    return success
                except Exception as e:
                    logger.error(f"Failed to create animation: {e}")
                    logger.error(f"Traceback: {traceback.format_exc()}")
                    plotter.close()
                    return False
            else:
                # Save screenshot
                try:
                    plotter.screenshot(sFilename_animation_out)
                    logger.info(f"✓ Visualization saved to: {sFilename_animation_out}")

                    # Verify file was created
                    if os.path.exists(sFilename_animation_out):
                        file_size = os.path.getsize(sFilename_animation_out)
                        logger.info(f"  File size: {file_size / 1024:.1f} KB")
                    else:
                        logger.warning(
                            f"Screenshot command executed but file not found: {sFilename_animation_out}"
                        )

                    plotter.close()
                    return True

                except Exception as e:
                    logger.error(f"Failed to save screenshot: {e}")
                    logger.error(f"Traceback: {traceback.format_exc()}")
                    plotter.close()
                    return False
        else:
            # Interactive display
            try:
                logger.info("Opening interactive visualization window...")
                plotter.show()
                return True
            except Exception as e:
                logger.error(f"Failed to display interactive visualization: {e}")
                logger.error(
                    f"Ensure display environment is available (X11, Wayland, etc.)"
                )
                logger.error(f"Traceback: {traceback.format_exc()}")
                plotter.close()
                return False

    except ImportError as e:
        logger.error(f"Missing required dependencies: {e}")
        return False

    except Exception as e:
        logger.error(f"Unexpected error during polyline visualization: {e}")
        logger.error(f"Error type: {type(e).__name__}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return False


def _create_rotation_animation(
    plotter,
    sFilename_out,
    dLongitude_start,
    dLatitude_focus,
    iAnimation_frames,
    dAnimation_speed,
    sAnimation_format,
):
    """
    Create a rotating animation of the 3D globe visualization with sine wave latitude pattern.
    Generates frames by rotating the camera around the globe while varying latitude in a
    sine wave pattern that visits Arctic and Antarctic regions before returning to start.
    Args:
        plotter: GeoVista plotter instance with mesh already added
        sFilename_out (str): Output animation file path
        dLongitude_start (float): Starting longitude for rotation
        dLatitude_focus (float): Starting latitude focus point (center of sine wave)
        iAnimation_frames (int): Number of frames for 360° rotation
        dAnimation_speed (float): Degrees per frame for longitude
        sAnimation_format (str): Output format ('mp4', 'gif', 'avi')
    Returns:
        bool: True if animation created successfully, False otherwise
    Note:
        Animation pattern: longitude rotates 360°, latitude follows sine wave with ±75° amplitude
        Sequence: start → Arctic → start → Antarctic → start (complete cycle)
    """
    try:
        # Import required libraries
        try:
            import imageio

            logger.info("ImageIO library imported successfully for animation creation")
            # Check available plugins/backends
            # Check specifically for video codecs
            try:
                import imageio_ffmpeg

                logger.info("FFmpeg backend available for MP4 creation")
            except ImportError:
                logger.warning("FFmpeg backend not available. MP4 creation may fail.")
                logger.warning("Install with: pip install imageio[ffmpeg]")
        except ImportError as e:
            logger.error(
                "ImageIO library not available. Install with: pip install imageio[ffmpeg]"
            )
            logger.error(f"Import error: {e}")
            return False
        import math

        # Validate animation parameters
        if iAnimation_frames <= 0:
            logger.error(f"Invalid number of animation frames: {iAnimation_frames}")
            return False
        # Ensure proper animation speed calculation
        if dAnimation_speed <= 0 or dAnimation_speed is None:
            dAnimation_speed = 360.0 / iAnimation_frames
            logger.info(
                f"Auto-calculated animation speed: {dAnimation_speed:.2f}° per frame"
            )
        else:
            logger.info(
                f"Using provided animation speed: {dAnimation_speed:.2f}° per frame"
            )
        # Validate output format and check backend availability
        valid_formats = ["mp4", "gif", "avi"]
        original_format = sAnimation_format.lower()
        if original_format not in valid_formats:
            logger.warning(f"Unsupported format {original_format}, defaulting to gif")
            sAnimation_format = "gif"
        else:
            sAnimation_format = original_format
        # delete any existing animation file
        if os.path.exists(sFilename_out):
            try:
                os.remove(sFilename_out)
                logger.info(f"Deleted existing animation file: {sFilename_out}")
            except Exception as e:
                logger.error(
                    f"Cannot delete existing animation file {sFilename_out}: {e}"
                )
                return False
        # Prepare output filename
        base_name = os.path.splitext(sFilename_out)[0]
        animation_filename = f"{base_name}.{sAnimation_format.lower()}"
        if sAnimation_format != original_format:
            logger.info(
                f"Output format changed from {original_format} to {sAnimation_format}"
            )
        # Use PyVista's built-in movie functionality - no temporary files needed
        logger.info(
            f"Creating {iAnimation_frames} frames for 360° rotation animation..."
        )
        logger.info(f"Animation will be saved as: {animation_filename}")
        # Use PyVista's built-in movie functionality
        logger.info("Creating animation using PyVista movie writer...")
        try:
            # Open movie file for writing
            plotter.open_movie(animation_filename, framerate=30)
            logger.info(f"Opened movie file: {animation_filename}")
            # Generate animation frames directly to movie
            for i in range(iAnimation_frames):
                # Calculate current longitude (rotate around globe)
                longitude_increment = (360.0 * i) / iAnimation_frames
                current_longitude = (dLongitude_start + longitude_increment) % 360.0
                # Normalize to [-180, 180] range
                if current_longitude > 180.0:
                    current_longitude -= 360.0
                # Calculate sine wave latitude pattern
                progress = i / iAnimation_frames
                latitude_amplitude = 75.0
                current_latitude = dLatitude_focus + latitude_amplitude * math.sin(
                    2 * math.pi * progress
                )
                current_latitude = max(-90.0, min(90.0, current_latitude))
                # Convert to radians and calculate camera position
                lon_rad = math.radians(current_longitude)
                lat_rad = math.radians(current_latitude)
                earth_radius = 1.0
                camera_distance = earth_radius * 3.0
                # Convert spherical coordinates to Cartesian
                x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
                y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
                z_focal = earth_radius * math.sin(lat_rad)
                x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
                y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
                z_camera = camera_distance * math.sin(lat_rad)
                focal_point = [x_focal, y_focal, z_focal]
                camera_position = [x_camera, y_camera, z_camera]
                # Update camera position
                plotter.camera.focal_point = focal_point
                plotter.camera.position = camera_position
                # plotter.camera.zoom(1.0) #apply zoom if needed, but be careful of cumulative zooming
                plotter.add_axes()  # Re-add axes to ensure visibility
                plotter.render()  # Render the scene
                plotter.write_frame()  # Write current frame to movie
                # Progress reporting
                if (i + 1) % max(1, iAnimation_frames // 10) == 0:
                    progress_pct = (i + 1) / iAnimation_frames * 100
                    logger.info(
                        f"  Frame {i+1}/{iAnimation_frames} ({progress_pct:.1f}%) - Lon: {current_longitude:.1f}°, Lat: {current_latitude:.1f}°"
                    )
                # Force garbage collection every 10 frames
                if (i + 1) % 10 == 0:
                    import gc

                    gc.collect()
            # Close movie file
            logger.info("Movie file closed")
            # Verify animation file was created
            if os.path.exists(animation_filename):
                file_size = os.path.getsize(animation_filename)
                logger.info(f"✓ Animation created successfully: {animation_filename}")
                logger.info(f"  File size: {file_size / (1024*1024):.2f} MB")
                logger.info(f"  Frames: {iAnimation_frames}")
                logger.info(f"  Format: {sAnimation_format.upper()}")
                return True
            else:
                logger.error("Animation file was not created")
                return False
        except Exception as e:
            logger.error(f"Failed to create animation using PyVista movie writer: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            # Try to close movie file if it was opened
            try:
                plotter.close_movie()
            except:
                pass
            return False
        finally:
            # Final cleanup to prevent memory leaks
            import gc

            gc.collect()
            logger.debug("Performed final garbage collection after animation creation")
    except Exception as e:
        logger.error(f"Unexpected error during animation creation: {e}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return False
    finally:
        # Ensure cleanup even if exceptions occur
        try:
            import gc

            gc.collect()
        except:
            pass
