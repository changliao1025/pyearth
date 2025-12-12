import os, sys
import logging
import math
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


# Constants for visualization
DEFAULT_EARTH_RADIUS = 1.0
DEFAULT_CAMERA_DISTANCE_MULTIPLIER = 3.0
DEFAULT_ZOOM_FACTOR = 0.7
VALID_ANIMATION_FORMATS = ['mp4', 'gif', 'avi']
VALID_IMAGE_FORMATS = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']
COORDINATE_BOUNDS = {'longitude': (-180, 180), 'latitude': (-90, 90)}

def animate_rotating_frames(pMesh, aValid_cell_indices: np.ndarray, pConfig: VisualizationConfig, pConfig_anima: AnimationConfig, sFilename_out: Optional[str]) -> bool:
    """
    Create a rotating animation of the 3D globe visualization.

    Args:
        pMesh: GeoVista mesh object
        aValid_cell_indices: Indices of valid cells to display
        pConfig (VisualizationConfig): Visualization configuration object.
        aConfig (AnimationConfig): Animation configuration object.
        sFilename_out (Optional[str]): Output animation file path (e.g., 'animation.mp4').

    Returns:
        bool: True if animation created successfully, False otherwise.
    """
    try:
        if pConfig.verbose:
            logger.info(f'Creating {pConfig_anima.frames}-frame rotation animation')
            logger.info(f'  - Starting longitude: {pConfig_anima.longitude_start:.1f}°')
            logger.info(f'  - Base latitude: {pConfig_anima.latitude_focus:.1f}°')
            logger.info(f'  - Rotation speed: {pConfig_anima.speed:.1f}°/frame')
            logger.info(f'  - Output format: {pConfig_anima.format}')

        # Setup plotter
        pPlotter = setup_geovista_plotter(iFlag_off_screen=(sFilename_out is not None), iFlag_verbose_in=pConfig.verbose)
        if pPlotter is None:
            return False

        # Add mesh to plotter
        pMesh_valid = pMesh.extract_cells(aValid_cell_indices)
        pPlotter.add_mesh(pMesh_valid, scalars=pConfig_anima.scalars, cmap=pConfig.colormap)

        earth_radius = DEFAULT_EARTH_RADIUS
        camera_distance = earth_radius * DEFAULT_CAMERA_DISTANCE_MULTIPLIER
        dAmplitude_deg = pConfig_anima.latitude_amplitude
        dCycles = pConfig_anima.latitude_cycles
        dPhase = pConfig_anima.latitude_phase

        # Initialize movie recording
        pPlotter.open_movie(sFilename_out, framerate=pConfig_anima.framerate)

        if pConfig.verbose:
            logger.info('Generating animation frames...')

        for iFrame in range(pConfig_anima.frames):
            # Calculate current longitude with smooth rotation
            dLongitude_current = pConfig_anima.longitude_start + (iFrame * pConfig_anima.speed)
            dLongitude_current = dLongitude_current % 360.0  # Keep within [0, 360)
            if dLongitude_current > 180.0:
                dLongitude_current -= 360.0  # Convert to [-180, 180]

            # Enhanced latitude movement: sine-wave oscillation for dynamic viewing
            dFrames_div = float(pConfig_anima.frames) if pConfig_anima.frames > 0 else 1.0
            dTheta = 2.0 * math.pi * (float(iFrame) / dFrames_div) * dCycles + dPhase
            dLatitude_current = float(pConfig_anima.latitude_focus) + dAmplitude_deg * math.sin(dTheta)

            # Clamp latitude to avoid pole singularities
            dLatitude_current = max(-89.9, min(89.9, dLatitude_current))

            # Convert to radians for calculations
            dLon_rad = math.radians(dLongitude_current)
            dLat_rad = math.radians(dLatitude_current)

            # Calculate focal point on Earth surface
            dX_focal = earth_radius * math.cos(dLat_rad) * math.cos(dLon_rad)
            dY_focal = earth_radius * math.cos(dLat_rad) * math.sin(dLon_rad)
            dZ_focal = earth_radius * math.sin(dLat_rad)

            # Calculate camera position away from Earth
            dX_camera = camera_distance * math.cos(dLat_rad) * math.cos(dLon_rad)
            dY_camera = camera_distance * math.cos(dLat_rad) * math.sin(dLon_rad)
            dZ_camera = camera_distance * math.sin(dLat_rad)

            aFocal_point = [dX_focal, dY_focal, dZ_focal]
            aCamera_position = [dX_camera, dY_camera, dZ_camera]

            # Update camera with smooth transitions
            pPlotter.camera.focal_point = aFocal_point
            pPlotter.camera.position = aCamera_position
            pPlotter.camera.up = [0, 0, 1]  # Maintain Z-up orientation

            # Apply zoom factor for consistent view
            pPlotter.camera.zoom(pConfig.zoom_factor)

            # Ensure axes remain visible throughout animation
            try:
                pPlotter.add_axes()
            except Exception:
                pass  # Axes may already exist

            # Render the current frame
            pPlotter.render()

            try:
                pPlotter.write_frame()

                if pConfig.verbose and (iFrame + 1) % max(1, pConfig_anima.frames // 10) == 0:
                    dProgress = ((iFrame + 1) / pConfig_anima.frames) * 100
                    logger.info(f'  Progress: {dProgress:.0f}% ({iFrame + 1}/{pConfig_anima.frames} frames)')

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
        if pConfig.verbose:
            logger.info(f'✓ Animation created successfully: {sFilename_out}')
            logger.info(f'  File size: {iFile_size / (1024*1024):.2f} MB')
            logger.info(f'  Frames: {pConfig_anima.frames}')
            logger.info(f'  Format: {pConfig_anima.format.upper()}')
            logger.info(f'  Duration: ~{pConfig_anima.frames / pConfig_anima.framerate:.1f} seconds at {pConfig_anima.framerate} FPS')

        return True

    except Exception as e:
        logger.error(f'Unexpected error during animation creation: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        try:
            pPlotter.close()
        except Exception:
            pass
        return False
