"""UTM (Universal Transverse Mercator) zone and EPSG code utilities.

References
----------
- UTM zones: 1-60, each covering 6° of longitude
- EPSG codes: 326XX for Northern Hemisphere, 327XX for Southern Hemisphere
"""

from osgeo import osr

osr.UseExceptions()

__all__ = [
    "get_utm_zone",
    "get_utm_epsg_code",
    "get_utm_spatial_reference_wkt",
]


def get_utm_zone(longitude: float, latitude: float) -> int:
    """Calculate the UTM zone number from longitude and latitude coordinates.

    The Universal Transverse Mercator (UTM) system divides the Earth into 60
    zones, each covering 6 degrees of longitude. Zone 1 starts at 180°W.

    Parameters
    ----------
    longitude : float
        Longitude in decimal degrees, range [-180, 180].
    latitude : float
        Latitude in decimal degrees, range [-80, 84].
        Note: UTM is not defined for latitudes beyond [-80, 84] (use UPS instead).

    Returns
    -------
    int
        UTM zone number (1-60).

    Raises
    ------
    ValueError
        If longitude is outside [-180, 180] range.
        If latitude is outside [-80, 84] range (UTM validity range).
    """
    try:
        lon = float(longitude)
        lat = float(latitude)
    except (TypeError, ValueError) as e:
        raise ValueError(
            f"Invalid coordinate values: longitude={longitude}, latitude={latitude}. Error: {e}"
        ) from e

    # Validate longitude
    if not -180.0 <= lon <= 180.0:
        raise ValueError(
            f"Longitude {lon}° is outside valid range [-180, 180]. "
            "Please provide longitude in decimal degrees."
        )

    # Validate latitude (UTM is defined for latitudes between 80°S and 84°N)
    if not -80.0 <= lat <= 84.0:
        raise ValueError(
            f"Latitude {lat}° is outside UTM validity range [-80, 84]. "
            "For polar regions (|lat| > 80°), use UPS (Universal Polar Stereographic) instead."
        )

    # Calculate zone: each zone is 6° wide, starting from -180°
    # Zone 1 is from -180° to -174°, Zone 31 is from 0° to 6°, etc.
    zone = int((lon + 180.0) / 6.0) + 1
    zone = max(1, min(60, zone))

    return zone


def get_utm_epsg_code(longitude: float, latitude: float) -> int:
    """Calculate the EPSG code for the UTM zone at given coordinates.

    EPSG codes for UTM zones follow a standard pattern:
    - Northern Hemisphere: 326XX (e.g., 32631 for Zone 31N)
    - Southern Hemisphere: 327XX (e.g., 32731 for Zone 31S)
    where XX is the zero-padded zone number (01-60).

    Parameters
    ----------
    longitude : float
        Longitude in decimal degrees, range [-180, 180].
    latitude : float
        Latitude in decimal degrees, range [-80, 84].

    Returns
    -------
    int
        EPSG code for the UTM zone (32601-32660 for North, 32701-32760 for South).

    Raises
    ------
    ValueError
        If longitude or latitude is outside valid UTM range.
    """
    utm_zone = get_utm_zone(longitude, latitude)
    zone_str = f"{utm_zone:02d}"

    if float(latitude) >= 0.0:
        return int("326" + zone_str)
    return int("327" + zone_str)


def get_utm_spatial_reference_wkt(longitude: float, latitude: float) -> str:
    """Create a WKT spatial reference string for the appropriate UTM zone.

    Generates a Well-Known Text (WKT) representation of the UTM spatial reference
    system for the given coordinates. The function automatically determines the
    correct UTM zone and hemisphere based on the input coordinates.

    Parameters
    ----------
    longitude : float
        Longitude in decimal degrees, range [-180, 180].
    latitude : float
        Latitude in decimal degrees, range [-80, 84].

    Returns
    -------
    str
        WKT (Well-Known Text) representation of the UTM spatial reference system.
        Based on WGS84 datum.

    Raises
    ------
    ValueError
        If longitude or latitude is outside valid UTM range.
    RuntimeError
        If the EPSG code cannot be imported or WKT export fails.
    """
    epsg_code = get_utm_epsg_code(longitude, latitude)

    spatial_ref = osr.SpatialReference()
    spatial_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    try:
        result = spatial_ref.ImportFromEPSG(epsg_code)
        if result != 0:
            raise RuntimeError(
                f"Failed to import EPSG code {epsg_code} for coordinates "
                f"(lon={longitude}, lat={latitude})."
            )
    except Exception as e:
        raise RuntimeError(
            f"Error importing EPSG:{epsg_code} for UTM zone at "
            f"(lon={longitude}, lat={latitude}): {e}"
        ) from e

    try:
        wkt_projection = spatial_ref.ExportToWkt()
        if not wkt_projection:
            raise RuntimeError(
                f"WKT export returned empty string for EPSG:{epsg_code}."
            )
        return wkt_projection
    except Exception as e:
        raise RuntimeError(f"Failed to export WKT for EPSG:{epsg_code}: {e}") from e

