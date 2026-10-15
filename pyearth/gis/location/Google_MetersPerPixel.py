"""Calculate ground resolution for Google Maps/Web Mercator tile pyramids."""

from typing import Optional, Union
import numpy as np

from pyearth.system.define_global_variables import earth_radius

__all__ = ["google_meters_per_pixel"]


def google_meters_per_pixel(
    zoom_level: Union[int, float], latitude: Optional[float] = None
) -> float:
    """Calculate ground resolution (meters per pixel) at a given zoom level.

    Computes ground resolution for Google Maps / Web Mercator tile pyramids
    (EPSG:3857). When latitude is provided, calculates the local ground resolution
    at that latitude; otherwise returns the resolution at the equator.

    Parameters
    ----------
    zoom_level : int or float
        Zoom level for the tile pyramid. Valid range is [0, 30].
        - Level 0: Entire world in a single 256×256 pixel tile (~156,543 m/pixel at equator)
        - Level 20: High resolution building detail (~0.15 m/pixel at equator)
    latitude : float, optional
        Latitude in decimal degrees [-85.0511, 85.0511] for local resolution.
        If None, computes resolution at the equator.

    Returns
    -------
    float
        Ground resolution in meters per pixel.

    Raises
    ------
    ValueError
        If zoom_level is negative, > 30, or if latitude is outside [-90, 90].
    TypeError
        If inputs are not numeric.
    """
    try:
        zoom = float(zoom_level)
    except (TypeError, ValueError) as e:
        raise TypeError(f"Zoom level must be numeric. Got {zoom_level}. Error: {e}")

    if zoom < 0:
        raise ValueError(f"Zoom level must be non-negative. Got {zoom}")

    if zoom > 30:
        raise ValueError(
            f"Zoom level {zoom} is excessively large (> 30). "
            "Typical maximum is 22-25 for most mapping services."
        )

    tile_size_pixels = 256.0

    # Ground resolution at equator: (2π * R) / (256 * 2^zoom)
    meters_per_pixel = (2.0 * np.pi * earth_radius) / (
        tile_size_pixels * np.power(2.0, zoom)
    )

    if latitude is not None:
        try:
            lat = float(latitude)
        except (TypeError, ValueError) as e:
            raise TypeError(f"Latitude must be numeric. Got {latitude}. Error: {e}")
        if not -90.0 <= lat <= 90.0:
            raise ValueError(f"Latitude must be in range [-90, 90]. Got {lat}")
        meters_per_pixel *= np.cos(np.radians(lat))

    return float(meters_per_pixel)

