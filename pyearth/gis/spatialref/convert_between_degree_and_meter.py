"""Coordinate resolution conversions between angular degrees and linear meters."""

from typing import Union
import numpy as np

from pyearth.system.define_global_variables import earth_radius

__all__ = ["degree_to_meter", "meter_to_degree"]


def degree_to_meter(
    resolution_degree: Union[float, np.ndarray] = None,
    latitude: Union[float, np.ndarray] = 0.0,
    **kwargs,
) -> Union[float, np.ndarray]:
    """Convert angular resolution in degrees to linear resolution in meters.

    Converts degrees to meters at a given latitude, accounting for the Earth's
    curvature. The conversion varies with latitude because lines of longitude
    converge at the poles.

    Parameters
    ----------
    resolution_degree : float or np.ndarray
        Resolution in degrees to convert.
    latitude : float or np.ndarray, optional
        Latitude in decimal degrees where the conversion applies.
        Default is 0.0 (equator).

    Returns
    -------
    float or np.ndarray
        Resolution in meters at the specified latitude.
    """
    if resolution_degree is None:
        resolution_degree = kwargs.get("dResolution_degree_in")
        if resolution_degree is None:
            raise ValueError("resolution_degree must be provided")

    if "dLatitude_mean_in" in kwargs:
        latitude = kwargs["dLatitude_mean_in"]

    res_arr = np.asarray(resolution_degree, dtype=float)
    lat_arr = np.asarray(latitude, dtype=float)

    latitude_rad = np.radians(lat_arr)
    radius_at_latitude = earth_radius * np.cos(latitude_rad)
    resolution_meter = res_arr * (2.0 * np.pi * radius_at_latitude) / 360.0

    if np.ndim(resolution_meter) == 0:
        return float(resolution_meter)
    return resolution_meter


def meter_to_degree(
    resolution_meter: Union[float, np.ndarray] = None,
    latitude: Union[float, np.ndarray] = 0.0,
    **kwargs,
) -> Union[float, np.ndarray]:
    """Convert linear resolution in meters to angular resolution in degrees.

    Converts meters to degrees at a given latitude, accounting for the Earth's
    curvature. The conversion varies with latitude because lines of longitude
    converge at the poles.

    Parameters
    ----------
    resolution_meter : float or np.ndarray
        Resolution in meters to convert.
    latitude : float or np.ndarray, optional
        Latitude in decimal degrees where the conversion applies.
        Default is 0.0 (equator).

    Returns
    -------
    float or np.ndarray
        Resolution in degrees at the specified latitude.

    Raises
    ------
    ValueError
        If latitude is too close to the poles (within 1 meter radius of axis).
    """
    if resolution_meter is None:
        resolution_meter = kwargs.get("dResolution_meter_in")
        if resolution_meter is None:
            raise ValueError("resolution_meter must be provided")

    if "dLatitude_mean_in" in kwargs:
        latitude = kwargs["dLatitude_mean_in"]

    res_arr = np.asarray(resolution_meter, dtype=float)
    lat_arr = np.asarray(latitude, dtype=float)

    latitude_abs = np.abs(lat_arr)
    latitude_rad = np.radians(latitude_abs)
    radius_at_latitude = earth_radius * np.cos(latitude_rad)

    # Check for near-zero radius (very close to poles)
    if np.any(radius_at_latitude < 1.0):
        min_radius = np.min(radius_at_latitude)
        raise ValueError(
            f"Latitude is too close to poles for accurate conversion. "
            f"Effective radius at this latitude: {min_radius:.2f} meters."
        )

    resolution_degree = res_arr * 360.0 / (2.0 * np.pi * radius_at_latitude)

    if np.ndim(resolution_degree) == 0:
        return float(resolution_degree)
    return resolution_degree

