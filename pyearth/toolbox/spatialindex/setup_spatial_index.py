"""
Setup spatial indexing library with fallback options.

This module provides utilities for selecting and initializing spatial indexing
libraries (tinyr or rtree) with automatic fallback support.
"""

import logging
import importlib.util
from typing import Tuple, Any

logger = logging.getLogger(__name__)


def setup_spatial_index() -> Tuple[Any, bool]:
    """
    Setup spatial indexing library with fallback options.

    This function attempts to use tinyr first (which requires Cython and is faster),
    then falls back to rtree if tinyr is not available. If neither library is
    available, it raises an ImportError with installation instructions.

    Returns:
        Tuple[Any, bool]: A tuple containing:
            - The spatial index class (RTree from tinyr or Index from rtree)
            - A boolean indicating if tinyr is being used (True) or rtree (False)

    Raises:
        ImportError: If neither tinyr nor rtree is available

    Example:
        >>> index_class, is_tinyr = setup_spatial_index()
        >>> if is_tinyr:
        ...     spatial_index = index_class(interleaved=True, max_cap=5, min_cap=2)
        ... else:
        ...     spatial_index = index_class()
    """
    try:
        # Try tinyr first (fast C implementation, requires Cython)
        if importlib.util.find_spec("tinyr") is not None:
            from tinyr import RTree
            logger.info("Using tinyr for spatial indexing")
            return RTree, True
    except ImportError:
        pass

    try:
        # Fallback to rtree
        if importlib.util.find_spec("rtree") is not None:
            import rtree.index
            logger.info("Using rtree for spatial indexing")
            return rtree.index.Index, False
    except ImportError:
        pass

    raise ImportError(
        "No spatial indexing library available. Please install either 'tinyr' (with Cython) or 'rtree'.\n"
        "Install with: pip install cython tinyr  OR  pip install rtree"
    )
