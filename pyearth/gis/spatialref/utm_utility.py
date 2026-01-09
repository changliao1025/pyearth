"""
UTM (Universal Transverse Mercator) zone and EPSG code utilities.

References
----------
- https://stackoverflow.com/a/40140326/4556479
- UTM zones: 1-60, each covering 6° of longitude
- EPSG codes: 326XX for Northern Hemisphere, 327XX for Southern Hemisphere
"""

from osgeo import osr


def get_utm_zone(lon: float, lat: float) -> int:
    """Calculate the UTM zone number from longitude and latitude coordinates.

    The Universal Transverse Mercator (UTM) system divides the Earth into 60
    zones, each covering 6 degrees of longitude. Zone 1 starts at 180°W.

    Parameters
    ----------
    lon : float
        Longitude in decimal degrees, range [-180, 180].
    lat : float
        Latitude in decimal degrees, range [-80, 84].
        Note: UTM is not defined for latitudes beyond ±80° (use UPS instead).

    Returns
    -------
    int
        UTM zone number (1-60).

    Raises
    ------
    ValueError
        If longitude is outside [-180, 180] range.
        If latitude is outside [-80, 84] range (UTM validity range).

    Notes
    -----
    - Special zones exist (e.g., Norway, Svalbard) but are not handled here
    - For latitudes above 84°N or below 80°S, use UPS (Universal Polar Stereographic)
    - Zone boundaries are at multiples of 6° starting from -180°

    Examples
    --------
    >>> get_utm_zone(0.0, 0.0)  # Prime meridian, equator
    31
    >>> get_utm_zone(-122.4194, 37.7749)  # San Francisco
    10
    >>> get_utm_zone(139.6917, 35.6895)  # Tokyo
    54
    """
    # Validate longitude
    if not -180 <= lon <= 180:
        raise ValueError(
            f"Longitude {lon}° is outside valid range [-180, 180]. "
            "Please provide longitude in decimal degrees."
        )

    # Validate latitude (UTM is defined for latitudes between 80°S and 84°N)
    if not -80 <= lat <= 84:
        raise ValueError(
            f"Latitude {lat}° is outside UTM validity range [-80, 84]. "
            "For polar regions (|lat| > 80°), use UPS (Universal Polar Stereographic) instead."
        )

    # Calculate zone: each zone is 6° wide, starting from -180°
    # Zone 1 is from -180° to -174°, Zone 31 is from 0° to 6°, etc.
    zone = int((lon + 180) / 6) + 1

    # Ensure zone is within valid range [1, 60]
    zone = max(1, min(60, zone))

    return zone


def get_utm_epsg_code(lon: float, lat: float) -> int:
    """Calculate the EPSG code for the UTM zone at given coordinates.

    EPSG codes for UTM zones follow a standard pattern:
    - Northern Hemisphere: 326XX (e.g., 32631 for Zone 31N)
    - Southern Hemisphere: 327XX (e.g., 32731 for Zone 31S)
    where XX is the zero-padded zone number (01-60).

    Parameters
    ----------
    lon : float
        Longitude in decimal degrees, range [-180, 180].
    lat : float
        Latitude in decimal degrees, range [-80, 84].

    Returns
    -------
    int
        EPSG code for the UTM zone (32601-32660 for North, 32701-32760 for South).

    Raises
    ------
    ValueError
        If longitude or latitude is outside valid UTM range.

    Notes
    -----
    - Uses WGS84 datum (EPSG:4326) as reference
    - Northern Hemisphere uses EPSG 326XX series
    - Southern Hemisphere uses EPSG 327XX series
    - Equator (lat=0) is assigned to Northern Hemisphere

    Examples
    --------
    >>> get_utm_epsg_code(0.0, 51.5074)  # London (Zone 31N)
    32631
    >>> get_utm_epsg_code(-122.4194, 37.7749)  # San Francisco (Zone 10N)
    32610
    >>> get_utm_epsg_code(-58.3816, -34.6037)  # Buenos Aires (Zone 21S)
    32721
    >>> get_utm_epsg_code(139.6917, 35.6895)  # Tokyo (Zone 54N)
    32654
    """
    # Get UTM zone (validates lon/lat internally)
    utm_zone = get_utm_zone(lon, lat)

    # Format zone number with zero-padding (01-60)
    zone_str = f"{utm_zone:02d}"

    # Determine hemisphere: Northern (326XX) or Southern (327XX)
    if lat >= 0:
        epsg_code = int("326" + zone_str)
    else:
        epsg_code = int("327" + zone_str)

    return epsg_code


def get_utm_spatial_reference_wkt(dLongitude_in: float, dLatitude_in: float) -> str:
    """Create a WKT spatial reference string for the appropriate UTM zone.

    Generates a Well-Known Text (WKT) representation of the UTM spatial reference
    system for the given coordinates. The function automatically determines the
    correct UTM zone and hemisphere based on the input coordinates.

    Parameters
    ----------
    dLongitude_in : float
        Longitude in decimal degrees, range [-180, 180].
    dLatitude_in : float
        Latitude in decimal degrees, range [-80, 84].
        Note: UTM is valid only within this latitude range.

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

    Notes
    -----
    - Uses WGS84 datum (EPSG:4326) as the geographic coordinate system
    - Automatically selects the correct UTM zone (1-60) based on longitude
    - Automatically selects hemisphere (North/South) based on latitude
    - UTM zones are only defined for latitudes between 80°S and 84°N

    Examples
    --------
    >>> wkt = get_utm_spatial_reference_wkt(-122.4194, 37.7749)  # San Francisco
    >>> 'UTM zone 10N' in wkt or '32610' in wkt
    True
    >>> wkt = get_utm_spatial_reference_wkt(139.6917, 35.6895)  # Tokyo
    >>> 'UTM zone 54N' in wkt or '32654' in wkt
    True
    """
    # Validate and convert inputs to float
    try:
        longitude = float(dLongitude_in)
        latitude = float(dLatitude_in)
    except (TypeError, ValueError) as e:
        raise ValueError(
            f"Invalid coordinate values: longitude={dLongitude_in}, latitude={dLatitude_in}. "
            f"Error: {e}"
        )

    # Get EPSG code (validates lon/lat internally)
    epsg_code = get_utm_epsg_code(longitude, latitude)

    # Create spatial reference from EPSG code
    spatial_ref = osr.SpatialReference()

    try:
        result = spatial_ref.ImportFromEPSG(epsg_code)
        if result != 0:
            raise RuntimeError(
                f"Failed to import EPSG code {epsg_code} for coordinates "
                f"(lon={longitude}, lat={latitude}). "
                "This may indicate an invalid EPSG code or OSR configuration issue."
            )
    except Exception as e:
        raise RuntimeError(
            f"Error importing EPSG:{epsg_code} for UTM zone at "
            f"(lon={longitude}, lat={latitude}): {e}"
        )

    # Export to WKT format
    try:
        wkt_projection = spatial_ref.ExportToWkt()
        if not wkt_projection:
            raise RuntimeError(
                f"WKT export returned empty string for EPSG:{epsg_code}."
            )
        return wkt_projection
    except Exception as e:
        raise RuntimeError(f"Failed to export WKT for EPSG:{epsg_code}: {e}")
