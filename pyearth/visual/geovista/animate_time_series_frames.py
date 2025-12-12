import os, sys
import logging
import numpy as np
from typing import Optional, List, Tuple, Union, Dict, Any
logger = logging.getLogger(__name__)
import traceback
from pyearth.visual.geovista.utility import (
    VisualizationConfig,
    AnimationConfig,
    setup_geovista_plotter,
    configure_camera,
    add_geographic_context,
)
def animate_multiple_frames(pPlotter, sFilename_out, dLongitude_start, dLatitude_focus,
                               dZoom_factor, iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose_in):
    """
    Create a rotating animation of the 3D globe visualization.

    Generates multiple frames by rotating the camera around the globe with enhanced
    camera movement patterns, then combines them into a video file.

    Args:
        pPlotter: GeoVista plotter instance with mesh already added
        sFilename_out (str): Output animation file path (e.g., 'animation.mp4')
        dLongitude_start (float): Starting longitude for rotation in degrees
        dLatitude_focus (float): Base latitude for camera focus in degrees
        dZoom_factor (float): Camera zoom level
        iAnimation_frames (int): Number of frames for 360° rotation
        dAnimation_speed (float): Degrees per frame
        sAnimation_format (str): Output format ('mp4', 'gif', 'avi')
        iFlag_verbose_in (bool): Enable verbose logging

    Returns:
        bool: True if animation created successfully, False otherwise
    """
    try:
        if iFlag_verbose_in:
            logger.info(f'Creating {iAnimation_frames}-frame rotation animation')
            logger.info(f'  - Starting longitude: {dLongitude_start:.1f}°')
            logger.info(f'  - Base latitude: {dLatitude_focus:.1f}°')
            logger.info(f'  - Rotation speed: {dAnimation_speed:.1f}°/frame')
            logger.info(f'  - Output format: {sAnimation_format}')

        # Animation parameters
        dEarth_radius = DEFAULT_EARTH_RADIUS
        dCamera_distance = dEarth_radius * DEFAULT_CAMERA_DISTANCE_MULTIPLIER
        dAmplitude_deg = 20.0  # Latitude oscillation amplitude
        dCycles = 1.0  # Number of sine cycles over full rotation
        dPhase = 0.0  # Phase shift for sine wave

        # Initialize movie recording
        pPlotter.open_movie(sFilename_out, framerate=30)

        if iFlag_verbose_in:
            logger.info('Generating animation frames...')

        for iFrame in range(iAnimation_frames):
            # Calculate current longitude with smooth rotation
            dLongitude_current = dLongitude_start + (iFrame * dAnimation_speed)
            dLongitude_current = dLongitude_current % 360.0  # Keep within [0, 360)
            if dLongitude_current > 180.0:
                dLongitude_current -= 360.0  # Convert to [-180, 180]

            # Enhanced latitude movement: sine-wave oscillation for dynamic viewing
            # This creates a more interesting camera path than fixed latitude
            dFrames_div = float(iAnimation_frames) if iAnimation_frames > 0 else 1.0
            dTheta = 2.0 * math.pi * (float(iFrame) / dFrames_div) * dCycles + dPhase
            dLatitude_current = float(dLatitude_focus) + dAmplitude_deg * math.sin(dTheta)

            # Clamp latitude to avoid pole singularities
            dLatitude_current = max(-89.9, min(89.9, dLatitude_current))

            # Convert to radians for calculations
            dLon_rad = math.radians(dLongitude_current)
            dLat_rad = math.radians(dLatitude_current)

            # Calculate focal point on Earth surface
            dX_focal = dEarth_radius * math.cos(dLat_rad) * math.cos(dLon_rad)
            dY_focal = dEarth_radius * math.cos(dLat_rad) * math.sin(dLon_rad)
            dZ_focal = dEarth_radius * math.sin(dLat_rad)

            # Calculate camera position away from Earth
            dX_camera = dCamera_distance * math.cos(dLat_rad) * math.cos(dLon_rad)
            dY_camera = dCamera_distance * math.cos(dLat_rad) * math.sin(dLon_rad)
            dZ_camera = dCamera_distance * math.sin(dLat_rad)

            aFocal_point = [dX_focal, dY_focal, dZ_focal]
            aCamera_position = [dX_camera, dY_camera, dZ_camera]

            # Update camera with smooth transitions
            pPlotter.camera.focal_point = aFocal_point
            pPlotter.camera.position = aCamera_position
            pPlotter.camera.up = [0, 0, 1]  # Maintain Z-up orientation

            # Apply zoom factor for consistent view
            pPlotter.camera.zoom(dZoom_factor)

            # Ensure axes remain visible throughout animation
            try:
                pPlotter.add_axes()
            except Exception:
                pass  # Axes may already exist

            # Render the current frame
            pPlotter.render()

            try:
                pPlotter.write_frame()

                if iFlag_verbose_in and (iFrame + 1) % max(1, iAnimation_frames // 10) == 0:
                    dProgress = ((iFrame + 1) / iAnimation_frames) * 100
                    logger.info(f'  Progress: {dProgress:.0f}% ({iFrame + 1}/{iAnimation_frames} frames)')

            except Exception as e:
                logger.error(f'Failed to render frame {iFrame + 1}: {e}')
                try:
                    pPlotter.close()
                except Exception:
                    pass
                return False

        # Close movie recording
        try:
            pPlotter.close()
        except Exception as e:
            logger.warning(f'Error closing plotter: {e}')

        # Validate output file creation
        if not os.path.exists(sFilename_out):
            logger.error('Animation file was not created')
            return False

        # Log success information
        iFile_size = os.path.getsize(sFilename_out)
        if iFlag_verbose_in:
            logger.info(f'✓ Animation created successfully: {sFilename_out}')
            logger.info(f'  File size: {iFile_size / (1024*1024):.2f} MB')
            logger.info(f'  Frames: {iAnimation_frames}')
            logger.info(f'  Format: {sAnimation_format.upper()}')
            logger.info(f'  Duration: ~{iAnimation_frames / 30:.1f} seconds at 30 FPS')

        return True

    except Exception as e:
        logger.error(f'Unexpected error during animation creation: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        try:
            pPlotter.close()
        except Exception:
            pass
        return False
