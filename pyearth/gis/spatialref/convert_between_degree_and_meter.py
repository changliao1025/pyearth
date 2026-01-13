import numpy as np
from pyearth.system.define_global_variables import earth_radius


def degree_to_meter(
    dResolution_degree_in: float, dLatitude_mean_in: float = 0.0
) -> float:
    """Convert angular resolution in degrees to linear resolution in meters.

    Converts degrees to meters at a given latitude, accounting for the Earth's
    curvature. The conversion varies with latitude because lines of longitude
    converge at the poles.

    Parameters
    ----------
    dResolution_degree_in : float
        Resolution in degrees to convert.
    dLatitude_mean_in : float, optional
        Mean latitude in degrees where the conversion applies.
        Default is 0.0 (equator).

    Returns
    -------
    float
        Resolution in meters at the specified latitude.

    Notes
    -----
    - Uses WGS84 Earth radius (6378137.0 meters)
    - At the equator (0°), one degree ≈ 111,319 meters
    - At 60° latitude, one degree ≈ 55,660 meters
    - The conversion accounts for the cosine of latitude to adjust for
      meridian convergence

    Examples
    --------
    >>> degree_to_meter(1.0, 0.0)  # One degree at equator
    111319.49079327358
    >>> degree_to_meter(1.0, 60.0)  # One degree at 60° latitude
    55659.74539663679
    """
    latitude_rad = np.radians(dLatitude_mean_in)
    radius_at_latitude = earth_radius * np.cos(latitude_rad)
    resolution_meter = dResolution_degree_in * (2 * np.pi * radius_at_latitude) / 360.0

    return float(resolution_meter)


def meter_to_degree(
    dResolution_meter_in: float, dLatitude_mean_in: float = 0.0
) -> float:
    """Convert linear resolution in meters to angular resolution in degrees.

    Converts meters to degrees at a given latitude, accounting for the Earth's
    curvature. The conversion varies with latitude because lines of longitude
    converge at the poles.

    Parameters
    ----------
    dResolution_meter_in : float
        Resolution in meters to convert.
    dLatitude_mean_in : float, optional
        Mean latitude in degrees where the conversion applies.
        Default is 0.0 (equator). The absolute value is used.

    Returns
    -------
    float
        Resolution in degrees at the specified latitude.

    Raises
    ------
    ValueError
        If the calculated radius at latitude is too close to zero
        (near the poles), which would cause division by zero.

    Notes
    -----
    - Uses WGS84 Earth radius (6378137.0 meters)
    - At the equator (0°), 111,319 meters ≈ one degree
    - At 60° latitude, 55,660 meters ≈ one degree
    - Takes the absolute value of latitude since the conversion is
      symmetric about the equator

    Examples
    --------
    >>> meter_to_degree(111319.0, 0.0)  # ~1 degree at equator
    0.9999955926593086
    >>> meter_to_degree(55660.0, 60.0)  # ~1 degree at 60° latitude
    1.0000046073406915
    """
    latitude_abs = abs(dLatitude_mean_in)
    latitude_rad = np.radians(latitude_abs)
    radius_at_latitude = earth_radius * np.cos(latitude_rad)

    # Check for near-zero radius (very close to poles)
    if radius_at_latitude < 1.0:  # Less than 1 meter radius
        raise ValueError(
            f"Latitude {dLatitude_mean_in}° is too close to poles for accurate conversion. "
            f"Effective radius at this latitude: {radius_at_latitude:.2f} meters."
        )

    resolution_degree = dResolution_meter_in * 360.0 / (2 * np.pi * radius_at_latitude)

    return float(resolution_degree)
