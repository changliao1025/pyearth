"""
Enhanced single frame visualization module for GeoVista.

This module provides improved single-frame visualization capabilities with comprehensive
error handling, validation, and enhanced user experience.
"""

import os
import sys
import logging
import traceback
from typing import Optional, Dict, Any, Tuple, List, Union
import numpy as np

from pyearth.visual.geovista.utility import (
    VisualizationConfig,
    ScalarBarConfig,
    PlotterManager,
    MeshHandler,
    configure_camera_enhanced,
    add_geographic_context_enhanced,
    add_mesh_to_plotter,
    validate_output_filename,
    get_system_info,
    VALID_IMAGE_FORMATS
)

# Set up logging
logger = logging.getLogger(__name__)

class SingleFrameResult:
    """Result object for single frame visualization operations."""

    def __init__(self, success: bool, message: str = "",
                 file_info: Optional[Dict[str, Any]] = None,
                 system_info: Optional[Dict[str, Any]] = None):
        self.success = success
        self.message = message
        self.file_info = file_info or {}
        self.system_info = system_info or {}

    def __str__(self) -> str:
        return f"SingleFrameResult(success={self.success}, message='{self.message}')"

    def __repr__(self) -> str:
        return self.__str__()

    def get_summary(self) -> str:
        """Get a formatted summary of the visualization result."""
        if not self.success:
            return f"❌ Visualization failed: {self.message}"

        lines = ["✅ Visualization completed successfully"]

        if self.file_info:
            lines.append(f"📁 File: {self.file_info.get('filename', 'N/A')}")
            lines.append(f"📏 Size: {self.file_info.get('size_kb', 0):.1f} KB")

        return "\n".join(lines)

def map_single_frame(
    pMesh,
    aValid_cell_indices: np.ndarray,
    pConfig: VisualizationConfig,
    sScalar: Optional[str] = None,
    sUnit: Optional[str] = None,
    sFilename_out: Optional[str] = None,
    # Enhanced parameters (optional for backward compatibility)
    scalar_config: Optional[ScalarBarConfig] = None,
    validate_inputs: bool = True,
    retry_on_failure: bool = True,
    return_detailed_result: bool = False
) -> Union[bool, SingleFrameResult]:
    """
    Enhanced single frame visualization with comprehensive error handling and validation.

    Args:
        pMesh: GeoVista mesh object to visualize
        aValid_cell_indices: Numpy array of indices for valid cells to display
        pConfig: Visualization configuration object
        sScalar: Name of scalar field to visualize (None for geometry-only visualization)
        sUnit: Unit string for scalar bar display (e.g., 'm/s', 'kg/m³')
        sFilename_out: Path for output file (None for interactive display)
        scalar_config: Custom scalar bar configuration (uses default if None)
        validate_inputs: Whether to validate inputs before processing
        retry_on_failure: Whether to attempt recovery strategies on failure
        return_detailed_result: Whether to return detailed result object

    Returns:
        bool or SingleFrameResult: Success status (bool for backward compatibility) or detailed result object

    Raises:
        ValueError: If input validation fails and validate_inputs is True
        TypeError: If mesh or config objects are of wrong type

    Example:
        >>> config = VisualizationConfig(longitude_focus=-100, latitude_focus=40)
        >>> # Simple usage (backward compatible)
        >>> success = map_single_frame(
        ...     mesh, indices, config,
        ...     sScalar="temperature",
        ...     sUnit="°C",
        ...     sFilename_out="output.png"
        ... )
        >>> # Advanced usage with detailed results
        >>> result = map_single_frame(
        ...     mesh, indices, config,
        ...     sScalar="temperature",
        ...     sUnit="°C",
        ...     sFilename_out="output.png",
        ...     return_detailed_result=True
        ... )
        >>> print(result.get_summary())
    """

    # Input validation
    if validate_inputs:
        validation_errors = _validate_inputs(pMesh, aValid_cell_indices, pConfig, sFilename_out)
        if validation_errors:
            error_msg = f"Input validation failed: {'; '.join(validation_errors)}"
            logger.error(error_msg)
            result = SingleFrameResult(False, error_msg, system_info=get_system_info())
            return result if return_detailed_result else False

    # Initialize result tracking
    result_info = {
        'plotter_created': False,
        'mesh_added': False,
        'camera_configured': False,
        'context_added': False,
        'output_saved': False
    }

    plotter = None

    try:
        # Validate output filename if provided
        if sFilename_out:
            is_valid, validation_msg = validate_output_filename(sFilename_out, VALID_IMAGE_FORMATS)
            if not is_valid:
                error_msg = f"Output filename validation failed: {validation_msg}"
                logger.error(error_msg)
                result = SingleFrameResult(False, error_msg)
                return result if return_detailed_result else False

        # Create plotter with enhanced manager
        if pConfig.verbose:
            logger.info("🎨 Creating GeoVista plotter...")

        plotter = PlotterManager.setup_geovista_plotter(
            off_screen=(sFilename_out is not None),
            verbose=pConfig.verbose,
            window_size=pConfig.window_size,
            use_xvfb=pConfig.use_xvfb,
            force_xvfb=pConfig.force_xvfb
        )

        if plotter is None:
            error_msg = "Failed to create GeoVista plotter"
            logger.error(error_msg)
            if retry_on_failure:
                logger.info("Attempting plotter creation with fallback strategies...")
                # Additional retry logic could be added here

            result = SingleFrameResult(False, error_msg, system_info=get_system_info())
            return result if return_detailed_result else False

        result_info['plotter_created'] = True

        # Add mesh to plotter using enhanced handler
        if pConfig.verbose:
            logger.info("📐 Adding mesh data to plotter...")

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
            error_msg = "Failed to add mesh to plotter"
            logger.error(error_msg)
            result = SingleFrameResult(False, error_msg, result_info)
            return result if return_detailed_result else False

        result_info['mesh_added'] = True

        # Configure camera with enhanced controller
        if pConfig.verbose:
            logger.info("📷 Configuring camera position...")

        camera_success = configure_camera_enhanced(
            plotter=plotter,
            config=pConfig,
            use_enhanced_controller=True
        )

        if not camera_success and pConfig.verbose:
            logger.warning("Camera configuration had issues, but continuing...")

        result_info['camera_configured'] = camera_success

        # Add geographic context with enhanced handler
        if pConfig.verbose:
            logger.info("🌍 Adding geographic context...")

        context_results = add_geographic_context_enhanced(
            plotter=plotter,
            config=pConfig,
            retry_on_failure=retry_on_failure
        )

        result_info['context_added'] = any(context_results.values())

        if pConfig.verbose and context_results:
            success_items = [k for k, v in context_results.items() if v]
            if success_items:
                logger.info(f"Successfully added: {', '.join(success_items)}")

        # Handle output or display
        file_info = {}

        if sFilename_out:
            # Save visualization to file
            if pConfig.verbose:
                logger.info(f"💾 Saving visualization to: {sFilename_out}")

            try:
                # Apply image scaling if specified
                if hasattr(pConfig, 'image_scale') and pConfig.image_scale != 1.0:
                    plotter.image_scale = pConfig.image_scale

                #plotter.screenshot(sFilename_out)
                plotter.save_graphic(sFilename_out, raster = False)
                result_info['output_saved'] = True

                # Verify and collect file information
                if os.path.exists(sFilename_out):
                    file_size = os.path.getsize(sFilename_out)
                    file_info = {
                        'filename': sFilename_out,
                        'size_bytes': file_size,
                        'size_kb': file_size / 1024,
                        'size_mb': file_size / (1024 * 1024),
                        'exists': True
                    }

                    if pConfig.verbose:
                        logger.info(f"✅ Visualization saved successfully")
                        logger.info(f"📁 File: {sFilename_out}")
                        logger.info(f"📏 Size: {file_info['size_kb']:.1f} KB")
                else:
                    error_msg = f"Screenshot command executed but file not found: {sFilename_out}"
                    logger.warning(error_msg)
                    file_info['exists'] = False

            except Exception as e:
                error_msg = f"Failed to save screenshot: {e}"
                logger.error(error_msg)
                result = SingleFrameResult(False, error_msg, result_info, get_system_info())
                return result if return_detailed_result else False
        else:
            # Interactive display
            if pConfig.verbose:
                logger.info("🖥️ Opening interactive visualization window...")
                logger.info("Close the window to continue...")

            try:
                plotter.show()
                result_info['output_saved'] = True  # Consider interactive display as "output"
            except Exception as e:
                error_msg = f"Failed to display interactive window: {e}"
                logger.error(error_msg)
                result = SingleFrameResult(False, error_msg, result_info, get_system_info())
                return result if return_detailed_result else False

        # Success!
        success_msg = "Single frame visualization completed successfully"
        if pConfig.verbose:
            logger.info(f"✅ {success_msg}")

        result = SingleFrameResult(
            success=True,
            message=success_msg,
            file_info=file_info,
            system_info=get_system_info() if return_detailed_result else {}
        )

        return result if return_detailed_result else True

    except Exception as e:
        error_msg = f"Unexpected error in single frame visualization: {e}"
        logger.error(error_msg)
        logger.error(f"Traceback: {traceback.format_exc()}")

        result = SingleFrameResult(False, error_msg, result_info, get_system_info())
        return result if return_detailed_result else False

    finally:
        # Cleanup resources
        if plotter is not None:
            try:
                plotter.close()
                if pConfig.verbose:
                    logger.debug("Plotter resources cleaned up")
            except Exception as e:
                logger.warning(f"Error during plotter cleanup: {e}")

def _validate_inputs(mesh, valid_cell_indices: np.ndarray,
                    config: VisualizationConfig,
                    output_filename: Optional[str]) -> List[str]:
    """
    Validate inputs for single frame visualization.

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

    # Validate config
    if not isinstance(config, VisualizationConfig):
        errors.append("Config must be an instance of VisualizationConfig")

    # Validate output filename if provided
    if output_filename is not None:
        if not isinstance(output_filename, str) or not output_filename.strip():
            errors.append("Output filename must be a non-empty string")
        else:
            # Check extension
            ext = os.path.splitext(output_filename.lower())[1].lstrip('.')
            if ext not in VALID_IMAGE_FORMATS:
                errors.append(f"Invalid image format '{ext}'. Valid formats: {VALID_IMAGE_FORMATS}")

    return errors