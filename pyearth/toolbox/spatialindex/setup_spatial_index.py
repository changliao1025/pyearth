
"""
Setup spatial indexing library using rtree.

This module provides a utility for selecting and initializing the rtree spatial indexing library.
"""

import logging
import importlib.util
from typing import Tuple, Any

logger = logging.getLogger(__name__)


def setup_spatial_index() -> Any:
    """
    Setup spatial indexing library using rtree.

    Returns
    -------
    Any
        The spatial index class (rtree.index.Index)

    Raises
    ------
    ImportError
        If rtree is not available.

    Example
    -------
    >>> index_class = setup_spatial_index()
    >>> spatial_index = index_class()
    """
    try:
        import rtree.index
        logger.info("Using rtree for spatial indexing")
        return rtree.index.Index
    except ImportError:
        raise ImportError(
            "No spatial indexing library available. Please install 'rtree'.\n"
            "Install with: pip install rtree"
        )
