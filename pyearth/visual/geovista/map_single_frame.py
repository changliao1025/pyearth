import os, sys
import logging
import numpy as np
from typing import Optional, List, Tuple, Union, Dict, Any
logger = logging.getLogger(__name__)
import traceback
from pyearth.visual.geovista.utility import (
    VisualizationConfig,
    setup_geovista_plotter,
    configure_camera,
    add_geographic_context,
)

def map_single_frame(pMesh, aValid_cell_indices: np.ndarray,
    pConfig: VisualizationConfig,
                     sScalar: Optional[str],
                                     sUnit: Optional[str],
                                     sFilename_out: Optional[str]) -> bool:
    """
    Handle single frame visualization (static image or interactive).

    Args:
        pMesh: GeoVista mesh object
        sScalars: Name of scalar field to visualize (optional, None for no color bar)
        aValid_cell_indices: Indices of valid cells to display
        sUnit: Unit string for colorbar (optional, None for no unit)
        pConfig: Visualization configuration
        sFilename: Output filename or None for interactive

    Returns:
        bool: True if successful, False otherwise
    """

    # Setup plotter
    pPlotter = setup_geovista_plotter(iFlag_off_screen=(sFilename_out is not None), iFlag_verbose_in=pConfig.verbose)
    if pPlotter is None:
        return False

    if sScalar:
        # Configure scalar bar
        dSargs = {
            "title": f"{sScalar} / {sUnit}" if sUnit else sScalar,
            "shadow": True,
            "title_font_size": 12,
            "label_font_size": 10,
            "fmt": "%.2f",
            "n_labels": 5,
        }
        # Add mesh with scalars to plotter
        pMesh_valid = pMesh.extract_cells(aValid_cell_indices)
        pPlotter.add_mesh(pMesh_valid, scalars=sScalar, scalar_bar_args=dSargs, cmap=pConfig.colormap)
    else:
        # Add mesh without scalars
        pMesh_valid = pMesh.extract_cells(aValid_cell_indices)
        pPlotter.add_mesh(pMesh_valid)

    # Configure camera and add geographic context
    configure_camera(pPlotter, pConfig)
    add_geographic_context(pPlotter, pConfig)

    try:
        if sFilename_out is not None:
            # Save screenshot
            pPlotter.screenshot(sFilename_out)
            if pConfig.verbose:
                logger.info(f'✓ Visualization saved to: {sFilename_out}')
                # Verify file was created
                if os.path.exists(sFilename_out):
                    iFile_size = os.path.getsize(sFilename_out)
                    logger.info(f'  File size: {iFile_size / 1024:.1f} KB')
                else:
                    logger.warning(f'Screenshot command executed but file not found: {sFilename_out}')
            pPlotter.close()
            return True
        else:
            # Interactive display
            if pConfig.verbose:
                logger.info('Opening interactive visualization window...')
            pPlotter.show()
            return True
    except Exception as e:
        logger.error(f'Error in single frame visualization: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False