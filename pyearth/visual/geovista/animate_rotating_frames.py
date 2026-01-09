"""
Enhanced rotating animation module for GeoVista.

This module provides improved animation capabilities with comprehensive error handling,
validation, progress tracking, and enhanced user experience.
"""

import os
import sys
import logging
import traceback
from typing import Optional, Dict, Any, List, Tuple, Callable, Union
import numpy as np

from pyearth.visual.geovista.utility import (
    VisualizationConfig,
    AnimationConfig,
    ScalarBarConfig,
    PlotterManager,
    MeshHandler,
    CameraController,
    configure_camera_enhanced,
    add_geographic_context_enhanced,
    add_mesh_to_plotter,
    validate_output_filename,
    get_system_info,
    VALID_ANIMATION_FORMATS
)

# Set up logging
logger = logging.getLogger(__name__)

class AnimationResult:
    """Result object for animation operations with comprehensive information."""

    def __init__(self, success: bool, message: str = "",
                 file_info: Optional[Dict[str, Any]] = None,
                 animation_info: Optional[Dict[str, Any]] = None,
                 frame_info: Optional[Dict[str, Any]] = None,
                 system_info: Optional[Dict[str, Any]] = None):
        self.success = success
        self.message = message
        self.file_info = file_info or {}
        self.animation_info = animation_info or {}
        self.frame_info = frame_info or {}
        self.system_info = system_info or {}

    def __str__(self) -> str:
        return f"AnimationResult(success={self.success}, message='{self.message}')"

    def __repr__(self) -> str:
        return self.__str__()

    def get_summary(self) -> str:
        """Get a formatted summary of the animation result."""
        if not self.success:
            return f"❌ Animation failed: {self.message}"

        lines = ["✅ Animation completed successfully"]

        if self.file_info:
            lines.append(f"📁 File: {self.file_info.get('filename', 'N/A')}")
            lines.append(f"📏 Size: {self.file_info.get('size_mb', 0):.2f} MB")

        if self.animation_info:
            lines.append(f"🎬 Frames: {self.animation_info.get('total_frames', 0)}")
            lines.append(f"⏱️ Duration: {self.animation_info.get('duration', 0):.1f}s")
            lines.append(f"🎯 FPS: {self.animation_info.get('framerate', 0)}")

        return "\n".join(lines)

def animate_rotating_frames(
    pMesh,
    aValid_cell_indices: np.ndarray,
    pConfig: VisualizationConfig,
    pConfig_anima: AnimationConfig,
    sScalar: Optional[str] = None,
    sUnit: Optional[str] = None,
    sFilename_out: Optional[str] = None,
    # Enhanced parameters (optional for backward compatibility)
    scalar_config: Optional[ScalarBarConfig] = None,
    validate_inputs: bool = True,
    retry_on_failure: bool = True,
    progress_callback: Optional[Callable[[int, int, float], None]] = None,
    return_detailed_result: bool = False
) -> Union[bool, AnimationResult]:
    """
    Create enhanced rotating animation of the 3D globe visualization.

    Args:
        pMesh: GeoVista mesh object to animate
        aValid_cell_indices: Indices of valid cells to display
        pConfig: Visualization configuration object
        pConfig_anima: Animation configuration object
        sScalar: Name of scalar field to visualize (None for geometry-only animation)
        sUnit: Unit string for scalar bar display (e.g., 'm/s', 'kg/m³')
        sFilename_out: Path for output animation file (required for animation)
        scalar_config: Custom scalar bar configuration (uses default if None)
        validate_inputs: Whether to validate inputs before processing
        retry_on_failure: Whether to attempt recovery strategies on failure
        progress_callback: Optional callback function for progress updates: callback(frame, total_frames, percentage)
        return_detailed_result: Whether to return detailed result object

    Returns:
        bool or AnimationResult: Success status (bool for backward compatibility) or detailed result object

    Raises:
        ValueError: If input validation fails and validate_inputs is True
        TypeError: If mesh, config objects are of wrong type, or required parameters are missing

    Example:
        >>> viz_config = VisualizationConfig(longitude_focus=-100, latitude_focus=40)
        >>> anim_config = AnimationConfig(frames=60, speed=6.0, framerate=30)
        >>> # Simple usage (backward compatible)
        >>> success = animate_rotating_frames(
        ...     mesh, indices, viz_config, anim_config,
        ...     sScalar="temperature",
        ...     sUnit="°C",
        ...     sFilename_out="animation.mp4"
        ... )
        >>> # Advanced usage with detailed results
        >>> result = animate_rotating_frames(
        ...     mesh, indices, viz_config, anim_config,
        ...     sScalar="temperature",
        ...     sUnit="°C",
        ...     sFilename_out="animation.mp4",
        ...     return_detailed_result=True
        ... )
        >>> print(result.get_summary())
    """

    # Input validation
    if validate_inputs:
        validation_errors = _validate_animation_inputs(
            pMesh, aValid_cell_indices, pConfig, pConfig_anima, sFilename_out
        )
        if validation_errors:
            error_msg = f"Animation input validation failed: {'; '.join(validation_errors)}"
            logger.error(error_msg)
            result = AnimationResult(False, error_msg, system_info=get_system_info())
            return result if return_detailed_result else False

    # Initialize progress tracking
    progress_info = {
        'plotter_created': False,
        'mesh_added': False,
        'camera_configured': False,
        'context_added': False,
        'movie_initialized': False,
        'frames_completed': 0,
        'animation_saved': False
    }
    import geovista as gv

    plotter = None
    start_time = None

    try:
        import time
        start_time = time.time()

        # Validate output filename
        if not sFilename_out:
            error_msg = "Output filename is required for animation"
            logger.error(error_msg)
            result = AnimationResult(False, error_msg)
            return result if return_detailed_result else False

        is_valid, validation_msg = validate_output_filename(sFilename_out, VALID_ANIMATION_FORMATS)
        if not is_valid:
            error_msg = f"Animation filename validation failed: {validation_msg}"
            logger.error(error_msg)
            result = AnimationResult(False, error_msg)
            return result if return_detailed_result else False

        # Display animation information
        if pConfig.verbose:
            logger.info("🎬 Starting enhanced rotating animation...")
            logger.info(f"📊 {pConfig_anima.get_animation_info()}")
            logger.info(f"🎯 Target file: {sFilename_out}")

        # Create plotter with enhanced manager
        plotter = PlotterManager.setup_geovista_plotter(
            off_screen=True,  # Always off-screen for animation
            verbose=pConfig.verbose,
            window_size=pConfig.window_size,
            use_xvfb=pConfig.use_xvfb,
            force_xvfb=pConfig.force_xvfb
        )

        if plotter is None:
            error_msg = "Failed to create GeoVista plotter for animation"
            logger.error(error_msg)
            result = AnimationResult(False, error_msg, system_info=get_system_info())
            return result if return_detailed_result else False

        plotter.add_base_layer(texture= gv.natural_earth_hypsometric())

        progress_info['plotter_created'] = True

        # Add mesh to plotter using enhanced handler
        if pConfig.verbose:
            logger.info("📐 Adding mesh data to animation plotter...")

        mesh_success = add_mesh_to_plotter(
            plotter=plotter,
            mesh=pMesh,
            valid_indices=aValid_cell_indices,
            scalar_name=sScalar,
            scalar_config=scalar_config,
            colormap=pConfig.colormap,
            unit=sUnit or "",
            validate_data=validate_inputs
        )

        if not mesh_success:
            error_msg = "Failed to add mesh to animation plotter"
            logger.error(error_msg)
            result = AnimationResult(False, error_msg, progress_info)
            return result if return_detailed_result else False

        progress_info['mesh_added'] = True

        # Configure initial camera position
        if pConfig.verbose:
            logger.info("📷 Configuring initial camera position...")

        camera_success = configure_camera_enhanced(
            plotter=plotter,
            config=pConfig,
            use_enhanced_controller=True
        )

        progress_info['camera_configured'] = camera_success

        # Add geographic context
        if pConfig.verbose:
            logger.info("🌍 Adding geographic context...")

        context_results = add_geographic_context_enhanced(
            plotter=plotter,
            config=pConfig,
            retry_on_failure=retry_on_failure
        )

        progress_info['context_added'] = any(context_results.values())

        # Initialize movie recording
        if pConfig.verbose:
            logger.info(f"🎥 Initializing movie recording at {pConfig_anima.framerate} FPS...")

        try:
            plotter.open_movie(sFilename_out, framerate=pConfig_anima.framerate)
            progress_info['movie_initialized'] = True
        except Exception as e:
            error_msg = f"Failed to initialize movie recording: {e}"
            logger.error(error_msg)
            result = AnimationResult(False, error_msg, progress_info)
            return result if return_detailed_result else False

        # Generate animation frames
        if pConfig.verbose:
            logger.info(f"🔄 Generating {pConfig_anima.frames} animation frames...")

        failed_frames = []

        for frame_idx in range(pConfig_anima.frames):
            try:
                # Calculate camera position for this frame using enhanced controller
                camera_pos = CameraController.calculate_animation_camera_position(
                    longitude=pConfig.longitude_focus,
                    latitude=pConfig.latitude_focus,
                    frame_idx=frame_idx,
                    config=pConfig_anima
                )

                # Apply camera position
                success = CameraController.apply_camera_to_plotter(plotter, camera_pos)
                if not success:
                    logger.warning(f"Camera positioning failed for frame {frame_idx + 1}, using fallback")

                # Render the frame
                plotter.render()
                plotter.write_frame()

                progress_info['frames_completed'] = frame_idx + 1

                # Progress reporting
                if progress_callback:
                    percentage = ((frame_idx + 1) / pConfig_anima.frames) * 100
                    progress_callback(frame_idx + 1, pConfig_anima.frames, percentage)

                # Periodic progress logging
                if pConfig.verbose and (frame_idx + 1) % max(1, pConfig_anima.frames // 10) == 0:
                    progress = ((frame_idx + 1) / pConfig_anima.frames) * 100
                    elapsed = time.time() - start_time
                    est_total = elapsed * pConfig_anima.frames / (frame_idx + 1)
                    remaining = est_total - elapsed
                    logger.info(f"  📈 Progress: {progress:.0f}% ({frame_idx + 1}/{pConfig_anima.frames} frames)")
                    logger.info(f"  ⏱️ Elapsed: {elapsed:.1f}s, Remaining: ~{remaining:.1f}s")

            except Exception as e:
                failed_frames.append(frame_idx + 1)
                logger.error(f"Failed to render frame {frame_idx + 1}: {e}")

                if not retry_on_failure:
                    error_msg = f"Frame rendering failed and retry disabled"
                    logger.error(error_msg)
                    result = AnimationResult(False, error_msg, progress_info)
                    return result if return_detailed_result else False

                # Continue with next frame for retry strategy
                continue

        # Close movie recording
        try:
            plotter.close()
            progress_info['animation_saved'] = True
        except Exception as e:
            logger.warning(f"Warning during plotter close: {e}")

        # Validate output file creation
        if not os.path.exists(sFilename_out):
            error_msg = f"Animation file was not created: {sFilename_out}"
            logger.error(error_msg)
            result = AnimationResult(False, error_msg, progress_info)
            return result if return_detailed_result else False

        # Collect file and animation information
        file_size = os.path.getsize(sFilename_out)
        end_time = time.time()
        total_time = end_time - start_time

        file_info = {
            'filename': sFilename_out,
            'size_bytes': file_size,
            'size_kb': file_size / 1024,
            'size_mb': file_size / (1024 * 1024),
            'exists': True,
            'format': pConfig_anima.format.upper()
        }

        animation_info = {
            'total_frames': pConfig_anima.frames,
            'completed_frames': progress_info['frames_completed'],
            'failed_frames': failed_frames,
            'framerate': pConfig_anima.framerate,
            'duration': pConfig_anima.frames / pConfig_anima.framerate,
            'generation_time': total_time,
            'frames_per_second_generated': pConfig_anima.frames / total_time if total_time > 0 else 0,
            'total_rotation_degrees': pConfig_anima.total_rotation,
            'speed_degrees_per_frame': pConfig_anima.speed
        }

        frame_info = {
            'successful': progress_info['frames_completed'] - len(failed_frames),
            'failed': len(failed_frames),
            'success_rate': ((progress_info['frames_completed'] - len(failed_frames)) / pConfig_anima.frames) * 100,
            'failed_frame_numbers': failed_frames
        }

        # Success logging
        if pConfig.verbose:
            logger.info("✅ Animation generation completed!")
            logger.info(f"📁 File: {sFilename_out}")
            logger.info(f"📏 Size: {file_info['size_mb']:.2f} MB")
            logger.info(f"🎬 Frames: {animation_info['completed_frames']}/{animation_info['total_frames']}")
            logger.info(f"⏱️ Duration: {animation_info['duration']:.1f} seconds at {animation_info['framerate']} FPS")
            logger.info(f"⚡ Generation rate: {animation_info['frames_per_second_generated']:.1f} frames/sec")

            if failed_frames:
                logger.warning(f"⚠️ Failed frames: {len(failed_frames)} out of {pConfig_anima.frames}")

        success_msg = f"Animation completed successfully with {frame_info['successful']}/{animation_info['total_frames']} frames"

        result = AnimationResult(
            success=True,
            message=success_msg,
            file_info=file_info,
            animation_info=animation_info,
            frame_info=frame_info,
            system_info=get_system_info() if return_detailed_result else {}
        )

        return result if return_detailed_result else True

    except Exception as e:
        error_msg = f"Unexpected error in animation generation: {e}"
        logger.error(error_msg)
        logger.error(f"Traceback: {traceback.format_exc()}")

        result = AnimationResult(False, error_msg, progress_info, get_system_info())
        return result if return_detailed_result else False

    finally:
        # Cleanup resources
        if plotter is not None:
            try:
                plotter.close()
                if pConfig.verbose:
                    logger.debug("Animation plotter resources cleaned up")
            except Exception as e:
                logger.warning(f"Error during animation plotter cleanup: {e}")

def _validate_animation_inputs(mesh, valid_cell_indices: np.ndarray,
                              viz_config: VisualizationConfig,
                              anim_config: AnimationConfig,
                              output_filename: Optional[str]) -> List[str]:
    """
    Validate inputs for animation generation.

    Returns:
        List of validation error messages (empty if all valid)
    """
    errors = []

    # Validate mesh
    if mesh is None:
        errors.append("Mesh cannot be None")

    # Validate indices
    if valid_cell_indices is None:
        errors.append("Valid cell indices cannot be None")
    elif not isinstance(valid_cell_indices, np.ndarray):
        errors.append("Valid cell indices must be a numpy array")
    elif len(valid_cell_indices) == 0:
        errors.append("Valid cell indices array cannot be empty")

    # Validate visualization config
    if not isinstance(viz_config, VisualizationConfig):
        errors.append("viz_config must be an instance of VisualizationConfig")

    # Validate animation config
    if not isinstance(anim_config, AnimationConfig):
        errors.append("anim_config must be an instance of AnimationConfig")
    else:
        # Validate animation parameters
        if anim_config.frames <= 0:
            errors.append("Animation frames must be positive")
        if anim_config.speed <= 0:
            errors.append("Animation speed must be positive")
        if anim_config.framerate <= 0:
            errors.append("Animation framerate must be positive")

    # Validate output filename
    if output_filename is not None:
        if not isinstance(output_filename, str) or not output_filename.strip():
            errors.append("Output filename must be a non-empty string")
        else:
            # Check extension
            ext = os.path.splitext(output_filename.lower())[1].lstrip('.')
            if ext not in VALID_ANIMATION_FORMATS:
                errors.append(f"Invalid animation format '{ext}'. Valid formats: {VALID_ANIMATION_FORMATS}")

    return errors

def create_animation_preview(
    mesh,
    valid_cell_indices: np.ndarray,
    viz_config: VisualizationConfig,
    anim_config: AnimationConfig,
    preview_frames: int = 5,
    scalar_name: Optional[str] = None,
    unit: Optional[str] = None,
    output_dir: str = "preview_frames"
) -> Dict[str, Any]:
    """
    Create preview frames for animation validation.

    Args:
        mesh: GeoVista mesh object
        valid_cell_indices: Valid cell indices
        viz_config: Visualization configuration
        anim_config: Animation configuration
        preview_frames: Number of preview frames to generate
        scalar_name: Optional scalar field name
        unit: Optional unit string
        output_dir: Directory to save preview frames

    Returns:
        Dict with preview information and generated file paths
    """
    preview_info = {
        'success': False,
        'generated_files': [],
        'errors': [],
        'frame_positions': []
    }

    try:
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Calculate frame positions to preview
        if preview_frames >= anim_config.frames:
            frame_indices = list(range(anim_config.frames))
        else:
            step = anim_config.frames // preview_frames
            frame_indices = [i * step for i in range(preview_frames)]
            if frame_indices[-1] != anim_config.frames - 1:
                frame_indices.append(anim_config.frames - 1)

        # Import here to avoid circular imports
        from .map_single_frame import map_single_frame_enhanced

        # Generate preview frames
        for i, frame_idx in enumerate(frame_indices):
            try:
                # Calculate camera position for this frame
                camera_pos = CameraController.calculate_animation_camera_position(
                    longitude=viz_config.longitude_focus,
                    latitude=viz_config.latitude_focus,
                    frame_idx=frame_idx,
                    config=anim_config
                )

                # Create modified config for this frame position
                frame_config = VisualizationConfig(
                    longitude_focus=camera_pos.focal_point[0],
                    latitude_focus=camera_pos.focal_point[1],
                    zoom_factor=viz_config.zoom_factor,
                    show_coastlines=viz_config.show_coastlines,
                    show_graticule=viz_config.show_graticule,
                    colormap=viz_config.colormap,
                    coastline_color=viz_config.coastline_color,
                    coastline_width=viz_config.coastline_width,
                    verbose=viz_config.verbose
                )

                # Generate preview frame
                output_file = os.path.join(output_dir, f"preview_frame_{i:03d}_of_{frame_idx:03d}.png")

                result = map_single_frame_enhanced(
                    mesh=mesh,
                    valid_cell_indices=valid_cell_indices,
                    config=frame_config,
                    scalar_name=scalar_name,
                    unit=unit,
                    output_filename=output_file,
                    validate_inputs=True,
                    return_detailed_result=True
                )

                if result.success:
                    preview_info['generated_files'].append(output_file)
                    preview_info['frame_positions'].append({
                        'frame_index': frame_idx,
                        'preview_number': i,
                        'filename': output_file,
                        'camera_position': camera_pos
                    })
                else:
                    preview_info['errors'].append(f"Failed to generate frame {frame_idx}: {result.message}")

            except Exception as e:
                preview_info['errors'].append(f"Error generating preview frame {frame_idx}: {e}")

        preview_info['success'] = len(preview_info['generated_files']) > 0

        if viz_config.verbose:
            logger.info(f"Generated {len(preview_info['generated_files'])} preview frames in {output_dir}")

    except Exception as e:
        preview_info['errors'].append(f"Preview generation failed: {e}")
        logger.error(f"Animation preview generation failed: {e}")

    return preview_info
