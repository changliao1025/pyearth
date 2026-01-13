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
import sys
import logging
import traceback
import math
import subprocess
from typing import Optional, List, Tuple, Union, Dict, Any, NamedTuple
import numpy as np
from osgeo import gdal, ogr
from functools import lru_cache

gdal.UseExceptions()

# Set up logging
logger = logging.getLogger(__name__)
CRS = "EPSG:4326"

# Constants for visualization
DEFAULT_EARTH_RADIUS = 1.0
DEFAULT_CAMERA_DISTANCE_MULTIPLIER = 3.0
DEFAULT_ZOOM_FACTOR = 0.7
DEFAULT_FRAMERATE = 30
DEFAULT_ANIMATION_FRAMES = 36
DEFAULT_ANIMATION_SPEED = 10.0  # degrees per frame
VALID_ANIMATION_FORMATS = ["mp4", "gif", "avi"]
VALID_IMAGE_FORMATS = ["png", "jpg", "jpeg", "svg", "tif", "tiff", "pdf", "ps"]
COORDINATE_BOUNDS = {"longitude": (-180, 180), "latitude": (-90, 90)}


# Camera position cache for performance optimization
class CameraPosition(NamedTuple):
    """Immutable camera position data structure."""

    focal_point: Tuple[float, float, float]
    camera_position: Tuple[float, float, float]


class VisualizationConfig:
    """Configuration class for visualization parameters."""

    def __init__(
        self,
        longitude_focus: float = 0.0,
        latitude_focus: float = 0.0,
        zoom_factor: float = DEFAULT_ZOOM_FACTOR,
        show_coastlines: bool = True,
        show_graticule: bool = True,
        colormap: str = "viridis",
        coastline_color: str = "black",
        coastline_width: float = 1.0,
        verbose: bool = False,
        window_size: Tuple[int, int] = (800, 600),
        image_scale: float = 1.0,
        use_xvfb: Optional[bool] = None,
        force_xvfb: bool = False,
    ):

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
        self.use_xvfb = use_xvfb  # None=auto-detect, True=force use, False=don't use
        self.force_xvfb = force_xvfb

    def _validate_longitude(self, lon: float) -> float:
        """Validate and clamp longitude to valid range."""
        if not (-180 <= lon <= 180):
            logger.warning(f"Longitude {lon} out of range [-180, 180], clamping")
            return np.clip(lon, -180, 180)
        return lon

    def _validate_latitude(self, lat: float) -> float:
        """Validate and clamp latitude to valid range."""
        if not (-90 <= lat <= 90):
            logger.warning(f"Latitude {lat} out of range [-90, 90], clamping")
            return np.clip(lat, -90, 90)
        return lat

    def _validate_zoom_factor(self, zoom: float) -> float:
        """Validate zoom factor."""
        if zoom <= 0:
            logger.warning(
                f"Invalid zoom factor {zoom}, using default {DEFAULT_ZOOM_FACTOR}"
            )
            return DEFAULT_ZOOM_FACTOR
        return zoom


class MeshHandler:
    """Utility class for handling mesh operations."""

    @staticmethod
    def validate_mesh_data(mesh, valid_indices: np.ndarray) -> bool:
        """Validate mesh and indices for visualization."""
        try:
            if mesh is None:
                logger.error("Mesh is None")
                return False

            if valid_indices is None or len(valid_indices) == 0:
                logger.error("No valid cell indices provided")
                return False

            # Check if indices are within bounds
            if hasattr(mesh, "n_cells") and np.max(valid_indices) >= mesh.n_cells:
                logger.error("Cell indices exceed mesh size")
                return False

            return True
        except Exception as e:
            logger.error(f"Error validating mesh data: {e}")
            return False

    @staticmethod
    def extract_valid_mesh(mesh, valid_indices: np.ndarray):
        """Extract valid cells from mesh with error handling."""
        try:
            if not MeshHandler.validate_mesh_data(mesh, valid_indices):
                return None
            return mesh.extract_cells(valid_indices)
        except Exception as e:
            logger.error(f"Error extracting valid mesh cells: {e}")
            return None

    @staticmethod
    def get_scalar_info(mesh, scalar_name: str) -> Dict[str, Any]:
        """Get information about scalar data in mesh."""
        try:
            if not hasattr(mesh, "array_names") or scalar_name not in mesh.array_names:
                return {"exists": False, "range": None, "dtype": None}

            scalar_data = mesh[scalar_name]
            return {
                "exists": True,
                "range": (np.nanmin(scalar_data), np.nanmax(scalar_data)),
                "dtype": scalar_data.dtype,
                "has_nan": np.any(np.isnan(scalar_data)),
                "shape": scalar_data.shape,
            }
        except Exception as e:
            logger.error(f"Error getting scalar info: {e}")
            return {"exists": False, "range": None, "dtype": None}


class CameraController:
    """Enhanced camera positioning and movement calculations."""

    _position_cache = {}  # Cache for frequently used positions

    @staticmethod
    @lru_cache(maxsize=128)
    def calculate_camera_position(
        longitude: float, latitude: float, zoom_factor: float = DEFAULT_ZOOM_FACTOR
    ) -> CameraPosition:
        """
        Calculate camera position and focal point from geographic coordinates with caching.

        Args:
            longitude: Longitude in degrees [-180, 180]
            latitude: Latitude in degrees [-90, 90]
            zoom_factor: Camera zoom level (positive float)

        Returns:
            CameraPosition: Named tuple containing focal_point and camera_position as tuples

        Raises:
            ValueError: If coordinates are out of valid ranges
        """
        # Validate inputs
        if not (-180 <= longitude <= 180):
            raise ValueError(f"Longitude {longitude} must be in range [-180, 180]")
        if not (-90 <= latitude <= 90):
            raise ValueError(f"Latitude {latitude} must be in range [-90, 90]")
        if zoom_factor <= 0:
            raise ValueError(f"Zoom factor {zoom_factor} must be positive")

        # Convert to radians
        lon_rad = math.radians(longitude)
        lat_rad = math.radians(latitude)

        # Calculate positions
        earth_radius = DEFAULT_EARTH_RADIUS
        camera_distance = (
            earth_radius * DEFAULT_CAMERA_DISTANCE_MULTIPLIER / zoom_factor
        )

        # Focal point on Earth surface
        x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
        y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
        z_focal = earth_radius * math.sin(lat_rad)

        # Camera position away from Earth
        x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
        y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
        z_camera = camera_distance * math.sin(lat_rad)

        focal_point = (x_focal, y_focal, z_focal)
        camera_position = (x_camera, y_camera, z_camera)

        return CameraPosition(focal_point, camera_position)

    @staticmethod
    def calculate_animation_camera_position(
        longitude: float, latitude: float, frame_idx: int, config: "AnimationConfig"
    ) -> CameraPosition:
        """
        Calculate camera position for animation frame with enhanced latitude movement.

        Args:
            longitude: Base longitude in degrees
            latitude: Base latitude in degrees
            frame_idx: Current frame index
            config: Animation configuration

        Returns:
            CameraPosition: Camera position for this frame
        """
        # Calculate current longitude with rotation
        longitude_current = config.longitude_start + (frame_idx * config.speed)
        longitude_current = longitude_current % 360.0  # Keep within [0, 360)
        if longitude_current > 180.0:
            longitude_current -= 360.0  # Convert to [-180, 180]

        # Enhanced latitude movement: sine-wave oscillation for dynamic viewing
        frames_div = float(config.frames) if config.frames > 0 else 1.0
        theta = (
            2.0 * math.pi * (float(frame_idx) / frames_div) * config.cycles
            + config.phase
        )
        latitude_current = float(
            config.latitude_start
        ) + config.amplitude_deg * math.sin(theta)

        # Clamp latitude to avoid pole singularities
        latitude_current = max(-89.9, min(89.9, latitude_current))

        # Convert to radians for calculations
        lon_rad = math.radians(longitude_current)
        lat_rad = math.radians(latitude_current)

        # Calculate positions
        earth_radius = DEFAULT_EARTH_RADIUS
        camera_distance = earth_radius * DEFAULT_CAMERA_DISTANCE_MULTIPLIER

        # Calculate focal point on Earth surface
        x_focal = earth_radius * math.cos(lat_rad) * math.cos(lon_rad)
        y_focal = earth_radius * math.cos(lat_rad) * math.sin(lon_rad)
        z_focal = earth_radius * math.sin(lat_rad)

        # Calculate camera position away from Earth
        x_camera = camera_distance * math.cos(lat_rad) * math.cos(lon_rad)
        y_camera = camera_distance * math.cos(lat_rad) * math.sin(lon_rad)
        z_camera = camera_distance * math.sin(lat_rad)

        focal_point = (x_focal, y_focal, z_focal)
        camera_position = (x_camera, y_camera, z_camera)

        return CameraPosition(focal_point, camera_position)

    @staticmethod
    def validate_camera_setup(
        focal_point: Union[List[float], Tuple[float, ...]],
        camera_position: Union[List[float], Tuple[float, ...]],
    ) -> bool:
        """
        Validate camera setup to ensure proper positioning.

        Args:
            focal_point: Focal point coordinates [x, y, z]
            camera_position: Camera position coordinates [x, y, z]

        Returns:
            bool: True if camera setup is valid
        """
        try:
            if len(focal_point) != 3 or len(camera_position) != 3:
                return False

            distance = math.sqrt(
                sum((c - f) ** 2 for c, f in zip(camera_position, focal_point))
            )
            return distance >= 0.1  # Minimum distance threshold
        except Exception:
            return False

    @staticmethod
    def apply_camera_to_plotter(
        plotter, camera_pos: CameraPosition, zoom_factor: float = 1.0
    ) -> bool:
        """
        Apply camera position to plotter with error handling.

        Args:
            plotter: GeoVista plotter instance
            camera_pos: Camera position data
            zoom_factor: Additional zoom factor

        Returns:
            bool: True if successful
        """
        try:
            plotter.camera.focal_point = camera_pos.focal_point
            plotter.camera.position = camera_pos.camera_position
            plotter.camera.up = [0, 0, 1]  # Z-up orientation
            if zoom_factor != 1.0:
                plotter.camera.zoom(zoom_factor)
            return True
        except Exception as e:
            logger.error(f"Failed to apply camera settings: {e}")
            return False


class ScalarBarConfig:
    """Configuration class for scalar bar parameters."""

    def __init__(
        self,
        title: str = "",
        title_font_size: int = 12,
        label_font_size: int = 10,
        fmt: str = "%.2f",
        n_labels: int = 5,
        shadow: bool = True,
        orientation: str = "horizontal",
        position_x: float = None,
        position_y: float = None,
        width: float = None,
        height: float = None,
    ):
        self.title = str(title)
        self.title_font_size = max(8, int(title_font_size))
        self.label_font_size = max(6, int(label_font_size))
        self.fmt = str(fmt)
        self.n_labels = max(2, int(n_labels))
        self.shadow = bool(shadow)

        # Validate and set orientation
        self.orientation = self._validate_orientation(orientation)

        # Set default positions based on orientation
        self.position_x, self.position_y, self.width, self.height = (
            self._set_default_dimensions(
                position_x, position_y, width, height, self.orientation
            )
        )

    def _validate_orientation(self, orientation: str) -> str:
        """Validate and normalize orientation parameter."""
        valid_orientations = ["horizontal", "vertical"]
        orientation_lower = str(orientation).lower()

        if orientation_lower not in valid_orientations:
            logger.warning(f'Invalid orientation "{orientation}", using "horizontal"')
            return "horizontal"

        return orientation_lower

    def _set_default_dimensions(
        self,
        position_x: Optional[float],
        position_y: Optional[float],
        width: Optional[float],
        height: Optional[float],
        orientation: str,
    ) -> Tuple[float, float, float, float]:
        """Set default dimensions based on orientation."""
        if orientation == "vertical":
            # Default values for vertical scalar bar
            default_position_x = 0.9
            default_position_y = 0.5
            default_width = 0.05
            default_height = 0.4
        else:  # horizontal
            # Default values for horizontal scalar bar
            default_position_x = 0.5
            default_position_y = 0.1
            default_width = 0.4
            default_height = 0.1

        # Use provided values or defaults
        final_position_x = np.clip(
            position_x if position_x is not None else default_position_x, 0.0, 1.0
        )
        final_position_y = np.clip(
            position_y if position_y is not None else default_position_y, 0.0, 1.0
        )
        final_width = np.clip(width if width is not None else default_width, 0.01, 0.5)
        final_height = np.clip(
            height if height is not None else default_height, 0.01, 1.0
        )

        return final_position_x, final_position_y, final_width, final_height

    def to_dict(self, scalar_name: str = "", unit: str = "") -> Dict[str, Any]:
        """Convert to dictionary format for GeoVista scalar_bar_args."""
        title = (
            f"{scalar_name} / {unit}"
            if unit and scalar_name
            else (scalar_name or self.title)
        )
        return {
            "title": title,
            "shadow": self.shadow,
            "title_font_size": self.title_font_size,
            "label_font_size": self.label_font_size,
            "fmt": self.fmt,
            "n_labels": self.n_labels,
            "position_x": self.position_x,
            "position_y": self.position_y,
            "width": self.width,
            "height": self.height,
            "vertical": (self.orientation == "vertical"),
        }


class AnimationConfig:
    """Enhanced configuration class for animation parameters."""

    def __init__(
        self,
        frames: Union[int, str] = DEFAULT_ANIMATION_FRAMES,
        speed: Union[float, str] = DEFAULT_ANIMATION_SPEED,
        format: str = "mp4",
        framerate: Union[int, str] = DEFAULT_FRAMERATE,
        amplitude_deg: float = 20.0,
        cycles: float = 1.0,
        phase: float = 0.0,
        longitude_start: float = 0.0,
        latitude_start: float = 0.0,
        latitude_focus: float = 0.0,
        # Legacy parameter support (deprecated)
        dLongitude_start: Optional[float] = None,
        dLatitude_start: Optional[float] = None,
    ):

        # Handle legacy parameters
        if dLongitude_start is not None:
            logger.warning(
                "Parameter 'dLongitude_start' is deprecated, use 'longitude_start'"
            )
            longitude_start = dLongitude_start
        if dLatitude_start is not None:
            logger.warning(
                "Parameter 'dLatitude_start' is deprecated, use 'latitude_start'"
            )
            latitude_start = dLatitude_start

        # Validate and set frames
        self.frames = self._validate_int_parameter(
            frames, "frames", DEFAULT_ANIMATION_FRAMES, min_val=1
        )

        # Validate and set speed (degrees per frame)
        self.speed = self._validate_float_parameter(
            speed, "speed", DEFAULT_ANIMATION_SPEED, min_val=0.1
        )

        # Validate and set framerate
        self.framerate = self._validate_int_parameter(
            framerate, "framerate", DEFAULT_FRAMERATE, min_val=1, max_val=120
        )

        # Validate format
        self.format = format.lower()
        if self.format not in VALID_ANIMATION_FORMATS:
            logger.warning(f"Invalid animation format {format}, using mp4")
            self.format = "mp4"

        # Animation motion parameters
        self.amplitude_deg = self._validate_float_parameter(
            amplitude_deg, "amplitude_deg", 20.0, min_val=0.0, max_val=90.0
        )
        self.cycles = self._validate_float_parameter(cycles, "cycles", 1.0, min_val=0.1)
        self.phase = math.radians(phase) if isinstance(phase, (int, float)) else 0.0

        # Starting position parameters
        self.longitude_start = self._validate_coordinate(longitude_start, "longitude")
        self.latitude_start = self._validate_coordinate(latitude_start, "latitude")
        self.latitude_focus = self._validate_coordinate(latitude_focus, "latitude")

        # Computed properties
        self.total_rotation = self.frames * self.speed
        self.estimated_duration = self.frames / self.framerate

    def _validate_int_parameter(
        self,
        value: Union[int, str],
        param_name: str,
        default: int,
        min_val: Optional[int] = None,
        max_val: Optional[int] = None,
    ) -> int:
        """Validate and convert integer parameters."""
        try:
            int_val = int(value) if isinstance(value, str) else int(value)
            if min_val is not None and int_val < min_val:
                logger.warning(
                    f"{param_name} {int_val} below minimum {min_val}, using {max(min_val, default)}"
                )
                return max(min_val, default)
            if max_val is not None and int_val > max_val:
                logger.warning(
                    f"{param_name} {int_val} above maximum {max_val}, using {min(max_val, default)}"
                )
                return min(max_val, default)
            return int_val
        except (ValueError, TypeError):
            logger.warning(
                f"Invalid {param_name} value {value}, using default {default}"
            )
            return default

    def _validate_float_parameter(
        self,
        value: Union[float, str],
        param_name: str,
        default: float,
        min_val: Optional[float] = None,
        max_val: Optional[float] = None,
    ) -> float:
        """Validate and convert float parameters."""
        try:
            float_val = float(value) if isinstance(value, str) else float(value)
            if min_val is not None and float_val < min_val:
                logger.warning(
                    f"{param_name} {float_val} below minimum {min_val}, using {max(min_val, default)}"
                )
                return max(min_val, default)
            if max_val is not None and float_val > max_val:
                logger.warning(
                    f"{param_name} {float_val} above maximum {max_val}, using {min(max_val, default)}"
                )
                return min(max_val, default)
            return float_val
        except (ValueError, TypeError):
            logger.warning(
                f"Invalid {param_name} value {value}, using default {default}"
            )
            return default

    def _validate_coordinate(self, coord: float, coord_type: str) -> float:
        """Validate coordinate values."""
        bounds = COORDINATE_BOUNDS.get(coord_type, (-180, 180))
        if not (bounds[0] <= coord <= bounds[1]):
            logger.warning(f"{coord_type} {coord} out of range {bounds}, clamping")
            return np.clip(coord, bounds[0], bounds[1])
        return coord

    @property
    def dLongitude_start(self) -> float:
        """Legacy property for backward compatibility."""
        return self.longitude_start

    @property
    def dLatitude_start(self) -> float:
        """Legacy property for backward compatibility."""
        return self.latitude_start

    def get_animation_info(self) -> str:
        """Get formatted animation information string."""
        return (
            f"Animation: {self.frames} frames, {self.framerate} FPS, "
            f"{self.estimated_duration:.1f}s duration, {self.total_rotation:.1f}° total rotation"
        )


def setup_geovista_plotter(
    iFlag_off_screen: bool = False, iFlag_verbose_in: bool = False
):
    """
    Legacy wrapper for backward compatibility.

    Args:
        iFlag_off_screen: Whether to create off-screen plotter
        iFlag_verbose_in: Enable verbose logging

    Returns:
        GeoVista plotter instance or None if failed
    """
    logger.warning(
        "setup_geovista_plotter is deprecated, use PlotterManager.setup_geovista_plotter"
    )
    return PlotterManager.setup_geovista_plotter(
        off_screen=iFlag_off_screen,
        verbose=iFlag_verbose_in,
        use_xvfb=None,  # Auto-detect
        force_xvfb=False,
    )


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
            pPlotter.add_coastlines(
                color=pConfig.coastline_color, line_width=pConfig.coastline_width
            )
            if pConfig.verbose:
                logger.debug(
                    f"Added coastlines overlay (color: {pConfig.coastline_color}, width: {pConfig.coastline_width})"
                )
        except Exception as e:
            logger.warning(f"Could not add coastlines: {e}")

    # Add coordinate axes
    try:
        pPlotter.add_axes()
        if pConfig.verbose:
            logger.debug("Added coordinate axes")
    except Exception as e:
        logger.warning(f"Could not add axes: {e}")

    # Add graticule (coordinate grid)
    if pConfig.show_graticule:
        try:
            pPlotter.add_graticule(show_labels=True)
            if pConfig.verbose:
                logger.debug("Added coordinate graticule with labels")
        except Exception as e:
            logger.warning(f"Could not add graticule: {e}")


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
            logger.warning("Camera and focal point are too close, using default view")
            raise ValueError("Invalid camera positioning")

        pPlotter.camera.focal_point = aFocal_point
        pPlotter.camera.position = aCamera_position
        pPlotter.camera.zoom(pConfig.zoom_factor)
        pPlotter.camera.up = [0, 0, 1]  # Z-up orientation

        if pConfig.verbose:
            logger.debug(f"Camera configured successfully:")
            logger.debug(f"  Focal point: {aFocal_point}")
            logger.debug(f"  Camera position: {aCamera_position}")

        return True

    except Exception as e:
        logger.warning(f"Error setting camera position: {e}. Using default view.")
        try:
            pPlotter.reset_camera()
        except Exception:
            pass
        return False


class PlotterManager:
    """Enhanced plotter management with error recovery."""

    @staticmethod
    def setup_geovista_plotter(
        off_screen: bool = False,
        verbose: bool = False,
        window_size: Tuple[int, int] = (800, 600),
        theme: str = "default",
        use_xvfb: Optional[bool] = None,
        force_xvfb: bool = False,
    ) -> Optional[Any]:
        """
        Set up GeoVista plotter with comprehensive error handling and recovery.

        Args:
            off_screen: Whether to create off-screen plotter
            verbose: Enable verbose logging
            window_size: Window size for interactive mode (width, height)
            theme: Plotter theme ('default', 'dark', 'paraview')
            use_xvfb: None=auto-detect, True=force use xvfb, False=don't use xvfb
            force_xvfb: Force xvfb even if not detected as needed

        Returns:
            GeoVista plotter instance or None if failed
        """
        # Handle xvfb setup for headless environments
        if use_xvfb is None:
            # Auto-detect if xvfb is needed
            should_use_xvfb = detect_headless_environment() or force_xvfb
        else:
            should_use_xvfb = use_xvfb

        if should_use_xvfb:
            if verbose:
                logger.info("Setting up Xvfb for headless environment...")
            xvfb_success = setup_xvfb_if_needed(force_xvfb=force_xvfb, verbose=verbose)
            if not xvfb_success:
                logger.warning(
                    "Failed to setup xvfb, attempting to continue without it"
                )

        try:
            import geovista as gv

            if verbose:
                logger.info("GeoVista library imported successfully")
        except ImportError as e:
            logger.error(
                "GeoVista library not available. Install with: pip install geovista"
            )
            logger.error(f"Import error: {e}")
            return None

        # Try multiple plotter creation strategies
        creation_strategies = [
            (
                "standard",
                lambda: PlotterManager._create_standard_plotter(
                    gv, off_screen, window_size, theme
                ),
            ),
            (
                "fallback",
                lambda: PlotterManager._create_fallback_plotter(gv, off_screen),
            ),
            ("minimal", lambda: PlotterManager._create_minimal_plotter(gv, off_screen)),
        ]

        for strategy_name, create_func in creation_strategies:
            try:
                plotter = create_func()
                if verbose:
                    logger.info(
                        f"Successfully created plotter using {strategy_name} strategy"
                    )
                return plotter
            except Exception as e:
                if verbose:
                    logger.warning(
                        f"Failed to create plotter with {strategy_name} strategy: {e}"
                    )
                continue

        logger.error("All plotter creation strategies failed")
        return None

    @staticmethod
    def _create_standard_plotter(
        gv, off_screen: bool, window_size: Tuple[int, int], theme: str
    ):
        """Create plotter with standard configuration."""
        kwargs = {"off_screen": off_screen}
        kwargs["window_size"] = window_size
        if theme != "default":
            kwargs["theme"] = theme
        return gv.GeoPlotter(**kwargs)

    @staticmethod
    def _create_fallback_plotter(gv, off_screen: bool):
        """Create plotter with minimal configuration."""
        return gv.GeoPlotter(off_screen=off_screen)

    @staticmethod
    def _create_minimal_plotter(gv, off_screen: bool):
        """Create the most basic plotter possible."""
        if off_screen:
            return gv.GeoPlotter(off_screen=True)
        else:
            return gv.GeoPlotter()


def add_geographic_context_enhanced(
    plotter, config: VisualizationConfig, retry_on_failure: bool = True
) -> Dict[str, bool]:
    """
    Add geographic context with comprehensive error handling and recovery.

    Args:
        plotter: GeoVista plotter instance
        config: Visualization configuration
        retry_on_failure: Whether to retry failed operations with alternative methods

    Returns:
        Dict indicating success/failure of each context element
    """
    results = {"coastlines": False, "axes": False, "graticule": False}

    # Add coastlines with fallback options
    if config.show_coastlines:
        coastline_strategies = [
            (
                "standard",
                lambda: plotter.add_coastlines(
                    color=config.coastline_color, line_width=config.coastline_width
                ),
            ),
            ("simple", lambda: plotter.add_coastlines()),
            ("basic", lambda: plotter.add_coastlines(color="black")),
        ]

        for strategy_name, add_func in coastline_strategies:
            try:
                add_func()
                results["coastlines"] = True
                if config.verbose:
                    logger.debug(f"Added coastlines using {strategy_name} method")
                break
            except Exception as e:
                if config.verbose:
                    logger.warning(f"Coastlines {strategy_name} method failed: {e}")
                if not retry_on_failure:
                    break
                continue

    # Add coordinate axes with fallback
    axes_strategies = [
        ("standard", lambda: plotter.add_axes()),
        (
            "simple",
            lambda: (
                plotter.add_axes(interactive=False)
                if hasattr(plotter, "add_axes")
                else None
            ),
        ),
    ]

    for strategy_name, add_func in axes_strategies:
        try:
            result = add_func()
            if result is not None or not hasattr(plotter, "add_axes"):
                results["axes"] = True
                if config.verbose:
                    logger.debug(f"Added axes using {strategy_name} method")
                break
        except Exception as e:
            if config.verbose:
                logger.warning(f"Axes {strategy_name} method failed: {e}")
            if not retry_on_failure:
                break
            continue

    # Add graticule with fallback options
    if config.show_graticule:
        graticule_strategies = [
            ("labeled", lambda: plotter.add_graticule(show_labels=True)),
            ("unlabeled", lambda: plotter.add_graticule(show_labels=False)),
            ("basic", lambda: plotter.add_graticule()),
        ]

        for strategy_name, add_func in graticule_strategies:
            try:
                add_func()
                results["graticule"] = True
                if config.verbose:
                    logger.debug(f"Added graticule using {strategy_name} method")
                break
            except Exception as e:
                if config.verbose:
                    logger.warning(f"Graticule {strategy_name} method failed: {e}")
                if not retry_on_failure:
                    break
                continue

    return results


def configure_camera_enhanced(
    plotter, config: VisualizationConfig, use_enhanced_controller: bool = True
) -> bool:
    """
    Configure camera position with enhanced error handling and recovery.

    Args:
        plotter: GeoVista plotter instance
        config: Visualization configuration
        use_enhanced_controller: Whether to use enhanced camera controller

    Returns:
        bool: True if camera configured successfully, False otherwise
    """
    if use_enhanced_controller:
        try:
            camera_pos = CameraController.calculate_camera_position(
                config.longitude_focus, config.latitude_focus, config.zoom_factor
            )

            if CameraController.validate_camera_setup(
                camera_pos.focal_point, camera_pos.camera_position
            ):
                success = CameraController.apply_camera_to_plotter(
                    plotter, camera_pos, config.zoom_factor
                )
                if success and config.verbose:
                    logger.debug(f"Enhanced camera configured successfully")
                    logger.debug(f"  Focal point: {camera_pos.focal_point}")
                    logger.debug(f"  Camera position: {camera_pos.camera_position}")
                return success
            else:
                logger.warning(
                    "Enhanced camera validation failed, trying legacy method"
                )
        except Exception as e:
            logger.warning(
                f"Enhanced camera configuration failed: {e}, trying legacy method"
            )

    # Legacy/fallback camera configuration
    try:
        # Calculate positions using legacy method for backward compatibility
        lon_rad = math.radians(config.longitude_focus)
        lat_rad = math.radians(config.latitude_focus)

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

        # Apply camera settings
        plotter.camera.focal_point = focal_point
        plotter.camera.position = camera_position
        plotter.camera.up = [0, 0, 1]  # Z-up orientation
        plotter.camera.zoom(config.zoom_factor)

        if config.verbose:
            logger.debug("Legacy camera configured successfully")

        return True

    except Exception as e:
        logger.error(
            f"All camera configuration methods failed: {e}. Using default view."
        )
        try:
            plotter.reset_camera()
            return True  # Default view is still usable
        except Exception:
            logger.error("Could not even reset camera to default view")
            return False


def add_mesh_to_plotter(
    plotter,
    mesh,
    valid_indices: np.ndarray,
    scalar_name: Optional[str] = None,
    scalar_config: Optional[ScalarBarConfig] = None,
    style: str = "surface",
    color: str = "white",
    colormap: str = "viridis",
    unit: str = "",
    opacity: float = 1.0,
    show_edges: bool = False,
    edge_color: str = "black",
    validate_data: bool = True,
) -> bool:
    """
    Add mesh to plotter with comprehensive validation and error handling.

    Args:
        plotter: GeoVista plotter instance
        mesh: Mesh object to visualize
        valid_indices: Indices of valid cells
        scalar_name: Name of scalar field to visualize
        scalar_config: Scalar bar configuration
        style: Mesh rendering style ('surface', 'wireframe', 'points')
        color: Mesh color (used when no scalars)
        colormap: Colormap name for scalar visualization
        unit: Unit string for scalar bar
        opacity: Mesh opacity (0.0 to 1.0)
        show_edges: Whether to show mesh edges
        edge_color: Color of mesh edges
        validate_data: Whether to validate mesh data before processing

    Returns:
        bool: True if successful
    """
    try:
        # Validate inputs
        if validate_data and not MeshHandler.validate_mesh_data(mesh, valid_indices):
            logger.error("Mesh validation failed")
            return False

        # Extract valid mesh
        mesh_valid = MeshHandler.extract_valid_mesh(mesh, valid_indices)
        if mesh_valid is None:
            logger.error("Failed to extract valid mesh")
            return False

        # Configure scalar visualization
        if scalar_name:
            # Get scalar information
            scalar_info = MeshHandler.get_scalar_info(mesh_valid, scalar_name)
            if not scalar_info["exists"]:
                logger.warning(
                    f"Scalar '{scalar_name}' not found in mesh, adding without scalars"
                )
                plotter.add_mesh(mesh_valid)
                return True

            # Configure scalar bar
            if scalar_config is None:
                scalar_config = ScalarBarConfig()

            scalar_args = scalar_config.to_dict(scalar_name, unit)

            # Add mesh with scalars
            plotter.add_mesh(
                mesh_valid,
                scalars=scalar_name,
                style=style,
                scalar_bar_args=scalar_args,
                cmap=colormap,
                opacity=opacity,
                show_edges=show_edges,
                edge_color=edge_color,
            )

            if scalar_info["has_nan"]:
                logger.warning(f"Scalar field '{scalar_name}' contains NaN values")
        else:
            # Add mesh without scalars
            plotter.add_mesh(
                mesh_valid,
                style=style,
                color=color,
                opacity=opacity,
                show_edges=show_edges,
                edge_color=edge_color,
            )

        return True

    except Exception as e:
        logger.error(f"Failed to add mesh to plotter: {e}")
        return False


# Utility functions for file validation and format checking
def validate_output_filename(
    filename: str, expected_formats: List[str] = None
) -> Tuple[bool, str]:
    """
    Validate output filename and format.

    Args:
        filename: Output filename to validate
        expected_formats: List of expected file extensions

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not filename or not isinstance(filename, str):
        return False, "Filename must be a non-empty string"

    if expected_formats is None:
        expected_formats = VALID_IMAGE_FORMATS + VALID_ANIMATION_FORMATS

    # Check file extension
    file_ext = os.path.splitext(filename.lower())[1].lstrip(".")
    if expected_formats and file_ext not in expected_formats:
        return False, f"Invalid file format '{file_ext}'. Expected: {expected_formats}"

    # Check directory exists or can be created
    directory = os.path.dirname(filename)
    if directory:
        try:
            os.makedirs(directory, exist_ok=True)
        except Exception as e:
            return False, f"Cannot create output directory: {e}"

    return True, ""


def detect_headless_environment() -> bool:
    """
    Detect if we're running in a headless Linux environment without graphics.

    Returns:
        bool: True if headless environment detected (needs xvfb), False otherwise
    """
    # Check if we're on Linux
    if sys.platform != "linux":
        return False

    # Check if DISPLAY environment variable is set
    display = os.environ.get("DISPLAY")
    if display:
        # DISPLAY is set, but check if it's accessible
        try:
            # Try to connect to the display
            import subprocess

            result = subprocess.run(["xdpyinfo"], capture_output=True, timeout=5)
            if result.returncode == 0:
                return False  # Display is accessible
        except (subprocess.TimeoutExpired, FileNotFoundError, ImportError):
            pass

    # Check for common headless environment indicators
    headless_indicators = [
        not display,  # No DISPLAY variable
        os.environ.get("SSH_CLIENT") is not None,  # SSH connection
        os.environ.get("SSH_TTY") is not None,  # SSH TTY
        os.environ.get("CI") == "true",  # CI environment
        os.environ.get("DEBIAN_FRONTEND") == "noninteractive",  # Non-interactive mode
    ]

    return any(headless_indicators)


def setup_xvfb_if_needed(force_xvfb: bool = False, verbose: bool = False) -> bool:
    """
    Set up Xvfb (X Virtual Framebuffer) if running in headless Linux environment.

    Args:
        force_xvfb: Force xvfb setup even if not detected as headless
        verbose: Enable verbose logging

    Returns:
        bool: True if xvfb was set up or already available, False if setup failed
    """
    try:
        import pyvista as pv
    except ImportError:
        if verbose:
            logger.warning("PyVista not available, cannot setup xvfb")
        return False

    # Check if we need xvfb
    needs_xvfb = force_xvfb or detect_headless_environment()

    if not needs_xvfb:
        if verbose:
            logger.debug("Display available, xvfb not needed")
        return True

    # Check if xvfb is already running (DISPLAY set by xvfb)
    if os.environ.get("DISPLAY") and not force_xvfb:
        if verbose:
            logger.debug("DISPLAY already set, assuming xvfb is running")
        return True

    try:
        if verbose:
            logger.info("🖥️ Setting up Xvfb for headless environment...")

        # Start xvfb using PyVista
        pv.start_xvfb()

        if verbose:
            logger.info("✅ Xvfb started successfully")

        return True

    except Exception as e:
        if verbose:
            logger.error(f"❌ Failed to start xvfb: {e}")
        logger.warning(f"Could not start xvfb: {e}")
        return False


def get_system_info() -> Dict[str, Any]:
    """Get system information for troubleshooting visualization issues."""
    info = {
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "platform": sys.platform,
        "os_name": os.name,
        "cwd": os.getcwd(),
        "display_available": os.environ.get("DISPLAY") is not None,
        "headless_detected": detect_headless_environment(),
        "geovista_available": False,
        "pyvista_available": False,
        "vtk_available": False,
    }

    try:
        import geovista

        info["geovista_available"] = True
        info["geovista_version"] = getattr(geovista, "__version__", "unknown")
    except ImportError:
        pass

    try:
        import pyvista

        info["pyvista_available"] = True
        info["pyvista_version"] = getattr(pyvista, "__version__", "unknown")
    except ImportError:
        pass

    try:
        import vtk

        info["vtk_available"] = True
        info["vtk_version"] = getattr(
            vtk, "vtkVersion", lambda: vtk.VTK_VERSION
        )().GetVTKVersion()
    except ImportError:
        pass

    return info
