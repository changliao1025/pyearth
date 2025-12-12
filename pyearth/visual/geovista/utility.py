"""
Visualization module for uraster class.

This module contains all visualization-related methods that were moved from the main uraster class
to reduce the size of the main uraster.py file and improve code organization.

Features:
- 3D mesh visualization using GeoVista
- Interactive and static rendering modes
- Animation support with rotation and camera movement
- Comprehensive error handling and validation
- Support for multiple output formats
"""

import os
import logging
import traceback
import math
from typing import Optional, List, Tuple, Union, Dict, Any
import numpy as np
from osgeo import gdal, ogr
gdal.UseExceptions()

# Set up logging
logger = logging.getLogger(__name__)
CRS = "EPSG:4326"

# Constants for visualization
DEFAULT_EARTH_RADIUS = 1.0
DEFAULT_CAMERA_DISTANCE_MULTIPLIER = 3.0
DEFAULT_ZOOM_FACTOR = 0.7
VALID_ANIMATION_FORMATS = ['mp4', 'gif', 'avi']
VALID_IMAGE_FORMATS = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']
COORDINATE_BOUNDS = {'longitude': (-180, 180), 'latitude': (-90, 90)}

class VisualizationConfig:
    """Configuration class for visualization parameters."""

    def __init__(self,
                 longitude_focus: float = 0.0,
                 latitude_focus: float = 0.0,
                 zoom_factor: float = DEFAULT_ZOOM_FACTOR,
                 show_coastlines: bool = True,
                 show_graticule: bool = True,
                 colormap: str = 'viridis',
                 coastline_color: str = 'black',
                 coastline_width: float = 1.0,
                 verbose: bool = False,
                 window_size: Tuple[int, int] = (800, 600),
                 image_scale: float = 1.0):
        self.longitude_focus = self._validate_longitude(longitude_focus)
        self.latitude_focus = self._validate_latitude(latitude_focus)
        self.zoom_factor = self._validate_zoom_factor(zoom_factor)
        self.show_coastlines = show_coastlines
        self.show_graticule = show_graticule
        self.colormap = colormap
        self.coastline_color = coastline_color
        self.coastline_width = coastline_width
        self.verbose = verbose
        self.window_size = window_size
        self.image_scale = image_scale

    def _validate_longitude(self, lon: float) -> float:
        """Validate and clamp longitude to valid range."""
        if not (-180 <= lon <= 180):
            logger.warning(f'Longitude {lon} out of range [-180, 180], clamping')
            return np.clip(lon, -180, 180)
        return lon

    def _validate_latitude(self, lat: float) -> float:
        """Validate and clamp latitude to valid range."""
        if not (-90 <= lat <= 90):
            logger.warning(f'Latitude {lat} out of range [-90, 90], clamping')
            return np.clip(lat, -90, 90)
        return lat

    def _validate_zoom_factor(self, zoom: float) -> float:
        """Validate zoom factor."""
        if zoom <= 0:
            logger.warning(f'Invalid zoom factor {zoom}, using default {DEFAULT_ZOOM_FACTOR}')
            return DEFAULT_ZOOM_FACTOR
        return zoom

class CameraController:
    """Handles camera positioning and movement calculations."""

    @staticmethod
    def calculate_camera_position(longitude: float, latitude: float,
                                zoom_factor: float = DEFAULT_ZOOM_FACTOR) -> Tuple[List[float], List[float]]:
        """
        Calculate camera position and focal point from geographic coordinates.

        Args:
            longitude: Longitude in degrees
            latitude: Latitude in degrees
            zoom_factor: Camera zoom level

        Returns:
            Tuple of (focal_point, camera_position) as [x, y, z] lists
        """
        # Convert to radians
        lon_rad = math.radians(longitude)
        lat_rad = math.radians(latitude)

        # Calculate positions
        earth_radius = DEFAULT_EARTH_RADIUS
        camera_distance = earth_radius * DEFAULT_CAMERA_DISTANCE_MULTIPLIER

        # Focal point on Earth surface
        x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
        y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
        z_focal = earth_radius * math.sin(lat_rad)

        # Camera position away from Earth
        x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
        y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
        z_camera = camera_distance * math.sin(lat_rad)

        focal_point = [x_focal, y_focal, z_focal]
        camera_position = [x_camera, y_camera, z_camera]

        return focal_point, camera_position

    @staticmethod
    def validate_camera_setup(focal_point: List[float], camera_position: List[float]) -> bool:
        """Validate camera setup to ensure proper positioning."""
        distance = math.sqrt(
            sum((c - f)**2 for c, f in zip(camera_position, focal_point))
        )
        return distance >= 0.1  # Minimum distance threshold

class AnimationConfig:
    """Configuration class for animation parameters."""

    def __init__(self,
                 frames: int = 36,
                 speed: float = 1.0,
                 format: str = 'mp4',
                 amplitude_deg: float = 20.0,
                 cycles: float = 1.0,
                 phase: float = 0.0,
                 dLongitude_start: float = 0.0,
                 dLatitude_start: float = 0.0):
        # Convert to int if string and validate
        try:
            frames_int = int(frames) if isinstance(frames, str) else frames
            self.frames = max(1, frames_int)  # Ensure at least 1 frame
        except (ValueError, TypeError):
            logger.warning(f'Invalid frames value {frames}, using default 36')
            self.frames = 36

        # Convert to float if string and validate
        try:
            speed_float = float(speed) if isinstance(speed, str) else speed
            self.speed = max(0.1, speed_float)  # Ensure positive speed
        except (ValueError, TypeError):
            logger.warning(f'Invalid speed value {speed}, using default 1.0')
            self.speed = 1.0

        self.format = format.lower()
        self.amplitude_deg = amplitude_deg
        self.cycles = cycles
        self.phase = phase
        self.dLongitude_start = dLongitude_start
        self.dLatitude_start = dLatitude_start

        # Validate format
        if self.format not in VALID_ANIMATION_FORMATS:
            logger.warning(f'Invalid animation format {format}, using mp4')
            self.format = 'mp4'

def setup_geovista_plotter(iFlag_off_screen: bool = False, iFlag_verbose_in: bool = False):
    """
    Set up GeoVista plotter with error handling.

    Args:
        iFlag_off_screen: Whether to create off-screen plotter
        iFlag_verbose_in: Enable verbose logging

    Returns:
        GeoVista plotter instance or None if failed
    """
    try:
        import geovista as gv
        if iFlag_verbose_in:
            logger.info('GeoVista library imported successfully')
    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return None

    try:
        if iFlag_off_screen:
            pPlotter = gv.GeoPlotter(off_screen=True)
            if iFlag_verbose_in:
                logger.debug('Created off-screen plotter')
        else:
            pPlotter = gv.GeoPlotter()
            if iFlag_verbose_in:
                logger.debug('Created interactive plotter')
        return pPlotter
    except Exception as e:
        logger.error(f'Failed to create GeoVista plotter: {e}')
        logger.error('This may be due to missing graphics context or display')
        if iFlag_off_screen:
            logger.error('For headless systems, ensure proper OpenGL/Mesa setup')
        else:
            logger.error('For interactive mode, ensure display environment (X11/Wayland) is available')
        return None

def add_geographic_context(pPlotter, pConfig: VisualizationConfig):
    """
    Add geographic context (coastlines, graticule, axes) to plotter.

    Args:
        pPlotter: GeoVista plotter instance
        pConfig: Visualization configuration
    """
    # Add coastlines
    if pConfig.show_coastlines:
        try:
            # You can set coastline color using the 'color' parameter
            # Common options: 'black', 'white', 'red', 'blue', 'gray', etc.
            # You can also use RGB tuples like (1.0, 0.0, 0.0) for red
            pPlotter.add_coastlines(color=pConfig.coastline_color, line_width=pConfig.coastline_width)
            if pConfig.verbose:
                logger.debug(f'Added coastlines overlay (color: {pConfig.coastline_color}, width: {pConfig.coastline_width})')
        except Exception as e:
            logger.warning(f'Could not add coastlines: {e}')

    # Add coordinate axes
    try:
        pPlotter.add_axes()
        if pConfig.verbose:
            logger.debug('Added coordinate axes')
    except Exception as e:
        logger.warning(f'Could not add axes: {e}')

    # Add graticule (coordinate grid)
    if pConfig.show_graticule:
        try:
            pPlotter.add_graticule(show_labels=True)
            if pConfig.verbose:
                logger.debug('Added coordinate graticule with labels')
        except Exception as e:
            logger.warning(f'Could not add graticule: {e}')

def configure_camera(pPlotter, pConfig: VisualizationConfig) -> bool:
    """
    Configure camera position and orientation.

    Args:
        pPlotter: GeoVista plotter instance
        pConfig: Visualization configuration

    Returns:
        bool: True if camera configured successfully, False otherwise
    """
    try:
        aFocal_point, aCamera_position = CameraController.calculate_camera_position(
            pConfig.longitude_focus, pConfig.latitude_focus, pConfig.zoom_factor
        )

        # Validate camera setup
        if not CameraController.validate_camera_setup(aFocal_point, aCamera_position):
            logger.warning('Camera and focal point are too close, using default view')
            raise ValueError('Invalid camera positioning')

        pPlotter.camera.focal_point = aFocal_point
        pPlotter.camera.position = aCamera_position
        pPlotter.camera.zoom(pConfig.zoom_factor)
        pPlotter.camera.up = [0, 0, 1]  # Z-up orientation

        if pConfig.verbose:
            logger.debug(f'Camera configured successfully:')
            logger.debug(f'  Focal point: {aFocal_point}')
            logger.debug(f'  Camera position: {aCamera_position}')

        return True

    except Exception as e:
        logger.warning(f'Error setting camera position: {e}. Using default view.')
        try:
            pPlotter.reset_camera()
        except Exception:
            pass
        return False

