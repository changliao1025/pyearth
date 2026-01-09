"""
Calculate ground resolution for Google Maps/Web Mercator tile pyramids.

This module provides utilities for calculating the meters per pixel (ground resolution)
at different zoom levels in the Google Maps tiling scheme (Web Mercator projection).
"""

import numpy as np

from pyearth.system.define_global_variables import earth_radius


def Google_MetersPerPixel(zoomLevel: int) -> float:
    """Calculate ground resolution (meters per pixel) at a given zoom level.

    Computes the ground resolution for Google Maps / Web Mercator tile pyramids
    at the equator. The resolution decreases (more detail) as zoom level increases.

    Parameters
    ----------
    zoomLevel : int
        Zoom level for the tile pyramid. Valid range is [0, 25].
        - Level 0: Entire world in a single 256×256 pixel tile
        - Level 1: World in 2×2 tiles (512×512 pixels total)
        - Higher levels: More tiles, more detail
        Common range: 0-20, with 21-25 for very high detail.

    Returns
    -------
    float
        Ground resolution in meters per pixel at the equator.
        At zoom level 0: ~156,543 meters/pixel
        At zoom level 20: ~0.15 meters/pixel

    Raises
    ------
    ValueError
        If zoom level is negative or excessively large (> 30).
    TypeError
        If zoom level cannot be converted to integer.

    Notes
    -----
    - Formula: resolution = (2π × R) / (256 × 2^zoom)
      where R is Earth's radius (WGS84 equatorial: 6,378,137 meters)
    - Resolution is calculated at the equator; actual resolution varies with latitude
    - At latitude φ, multiply by cos(φ) to get local resolution
    - Uses Web Mercator projection (EPSG:3857)
    - Standard tile size is 256×256 pixels
    - Google Maps typically supports zoom levels 0-22

    Resolution at Common Zoom Levels (at equator)
    ----------------------------------------------
    - Zoom 0: ~156,543 m/pixel (world view)
    - Zoom 5: ~4,892 m/pixel (country view)
    - Zoom 10: ~153 m/pixel (city view)
    - Zoom 15: ~4.77 m/pixel (street view)
    - Zoom 20: ~0.15 m/pixel (building detail)

    Examples
    --------
    >>> # World view (zoom 0)
    >>> Google_MetersPerPixel(0)
    156543.03392804097

    >>> # City view (zoom 10)
    >>> Google_MetersPerPixel(10)
    152.8740565703525

    >>> # Street view (zoom 15)
    >>> Google_MetersPerPixel(15)
    4.777314267823516

    >>> # Building detail (zoom 20)
    >>> Google_MetersPerPixel(20)
    0.14929107087573174

    >>> # Adjust for latitude (e.g., 45° N)
    >>> resolution_equator = Google_MetersPerPixel(10)
    >>> resolution_45N = resolution_equator * np.cos(np.radians(45))
    >>> resolution_45N
    108.10293268405542

    See Also
    --------
    numpy.power : Power function used for zoom calculation

    References
    ----------
    - Google Maps/Bing Maps Tile System:
      https://docs.microsoft.com/en-us/bingmaps/articles/bing-maps-tile-system
    - Web Mercator projection: https://epsg.io/3857
    """
    # Validate and convert zoom level
    try:
        zoom = int(zoomLevel)
    except (TypeError, ValueError) as e:
        raise TypeError(f"Zoom level must be an integer. Got {zoomLevel}. Error: {e}")

    # Validate zoom level range
    if zoom < 0:
        raise ValueError(f"Zoom level must be non-negative. Got {zoom}")

    if zoom > 30:
        raise ValueError(
            f"Zoom level {zoom} is excessively large (> 30). "
            "Typical maximum is 22-25 for most mapping services."
        )

    # Constants
    TILE_SIZE_PIXELS = 256  # Standard tile size in pixels

    # Calculate meters per pixel at equator
    # Formula: (2π × Earth_radius) / (tile_size × 2^zoom)
    # This is the circumference of Earth divided by total pixels at this zoom level
    meters_per_pixel = (2 * np.pi * earth_radius) / (
        TILE_SIZE_PIXELS * np.power(2, zoom)
    )

    return float(meters_per_pixel)
