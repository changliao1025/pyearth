"""Great-circle distance calculations on spherical Earth."""

from typing import Optional, Union

import numpy as np

from pyearth.system.define_global_variables import earth_radius

__all__ = [
    "calculate_distance_based_on_longitude_latitude",
]


def calculate_distance_based_on_longitude_latitude(
    aLongitude_from: Union[float, list, np.ndarray],
    aLatitude_from: Union[float, list, np.ndarray],
    aLongitude_to: Union[float, list, np.ndarray],
    aLatitude_to: Union[float, list, np.ndarray],
    bUnits_are_radians: bool = False,
    dRadius_in: Optional[float] = None,
    iFlag_radian: Optional[bool] = None,
    radius: Optional[float] = None,
) -> Union[float, np.ndarray]:
    """Calculate the great circle distance using the haversine formula.

    Supports both scalar coordinates and NumPy arrays.

    Parameters
    ----------
    aLongitude_from : float or array_like
        Longitudes of the start points.
    aLatitude_from : float or array_like
        Latitudes of the start points.
    aLongitude_to : float or array_like
        Longitudes of the end points.
    aLatitude_to : float or array_like
        Latitudes of the end points.
    bUnits_are_radians : bool, optional
        If True, coordinates are in radians. If False (default), in degrees.
    dRadius_in : float, optional
        Sphere radius in meters. If None, default earth_radius is used.

    Returns
    -------
    float or np.ndarray
        Great circle distances in meters (or in radians if bUnits_are_radians=True).
    """
    lon_from = np.asarray(aLongitude_from)
    lat_from = np.asarray(aLatitude_from)
    lon_to = np.asarray(aLongitude_to)
    lat_to = np.asarray(aLatitude_to)

    units_radians = iFlag_radian if iFlag_radian is not None else bUnits_are_radians

    if not units_radians:
        lon1 = np.radians(lon_from)
        lat1 = np.radians(lat_from)
        lon2 = np.radians(lon_to)
        lat2 = np.radians(lat_to)
    else:
        lon1 = lon_from
        lat1 = lat_from
        lon2 = lon_to
        lat2 = lat_to

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    c = 2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))

    if units_radians:
        return float(c) if c.ndim == 0 else c

    rad_val = radius if radius is not None else dRadius_in
    eff_radius = float(rad_val) if rad_val is not None else float(earth_radius)
    d = c * eff_radius
    return float(d) if d.ndim == 0 else d
