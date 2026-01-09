"""
Determine HydroSHEDS continent codes from geographic extents.

This module provides utilities for identifying which HydroSHEDS continent region(s)
overlap with a given geographic bounding box. HydroSHEDS uses 2-letter continent
codes (af, as, au, etc.) for organizing global hydrological datasets.
"""

from typing import Tuple, List, Optional, Union


def _box_area(box: Tuple[float, float, float, float]) -> float:
    """Calculate area of an axis-aligned bounding box.

    Parameters
    ----------
    box : Tuple[float, float, float, float]
        Bounding box as (minx, miny, maxx, maxy) in degrees.

    Returns
    -------
    float
        Area in square degrees (width × height).

    Notes
    -----
    This is a simple planar approximation suitable for continent-scale
    bounding boxes. For precise area calculations, use spherical geometry.
    """
    minx, miny, maxx, maxy = box
    width = max(0.0, maxx - minx)
    height = max(0.0, maxy - miny)
    return width * height


def _intersection_area(
    b1: Tuple[float, float, float, float], b2: Tuple[float, float, float, float]
) -> float:
    """Compute intersection area of two axis-aligned bounding boxes.

    Parameters
    ----------
    b1 : Tuple[float, float, float, float]
        First bounding box as (minx, miny, maxx, maxy).
    b2 : Tuple[float, float, float, float]
        Second bounding box as (minx, miny, maxx, maxy).

    Returns
    -------
    float
        Intersection area in square degrees, or 0.0 if boxes don't overlap.

    Notes
    -----
    Uses standard axis-aligned bounding box intersection algorithm.
    """
    minx = max(b1[0], b2[0])
    miny = max(b1[1], b2[1])
    maxx = min(b1[2], b2[2])
    maxy = min(b1[3], b2[3])

    # No intersection if boxes don't overlap
    if maxx <= minx or maxy <= miny:
        return 0.0

    return (maxx - minx) * (maxy - miny)


def get_hydrosheds_continent_from_extent(
    aExtent: Tuple[float, float, float, float],
    threshold_fraction: float = 0.10,
    return_all: bool = False,
) -> Optional[Union[str, List[str]]]:
    """Determine HydroSHEDS continent code(s) overlapping a geographic extent.

    Uses coarse continent bounding boxes to identify which HydroSHEDS continent
    region(s) overlap with the input extent. Returns 2-letter continent codes
    used by HydroSHEDS datasets.

    Parameters
    ----------
    aExtent : Tuple[float, float, float, float]
        Bounding box as (min_lon, max_lon, min_lat, max_lat) in WGS84 degrees.
        Longitude range: [-180, 180] or [0, 360] (auto-normalized).
        Latitude range: [-90, 90].
    threshold_fraction : float, optional
        When return_all=True, minimum overlap fraction required to include
        a continent. Overlap area must be >= threshold_fraction × input_area.
        Default is 0.10 (10% overlap required).
    return_all : bool, optional
        If True, return list of all continents meeting threshold.
        If False (default), return only the continent with largest overlap.

    Returns
    -------
    str, List[str], or None
        If return_all=False: 2-letter continent code (e.g., "na", "af") or None.
        If return_all=True: List of 2-letter codes sorted by overlap area, or None.
        Returns None if no continent overlaps or extent is invalid.

    Raises
    ------
    ValueError
        If threshold_fraction is not in range [0, 1].
        If latitude values are outside [-90, 90].

    Notes
    -----
    HydroSHEDS Continent Codes:
    - af: Africa
    - ar: Arctic
    - as: Asia
    - au: Australia
    - eu: Europe
    - gr: Greenland
    - na: North America
    - sa: South America
    - si: Siberia

    Limitations:
    - Uses coarse bounding boxes; not precise polygon boundaries
    - Arctic region box is approximate (western hemisphere only)
    - Some regions may overlap (e.g., Europe/Asia boundary)
    - Handles dateline crossing (longitude wrapping)
    - Uses planar approximation for area calculations

    Examples
    --------
    >>> # Texas, USA
    >>> extent = (-106.0, -93.5, 25.8, 36.5)  # (min_lon, max_lon, min_lat, max_lat)
    >>> get_hydrosheds_continent_from_extent(extent)
    'na'

    >>> # Western Europe
    >>> extent = (-10.0, 30.0, 35.0, 60.0)
    >>> get_hydrosheds_continent_from_extent(extent)
    'eu'

    >>> # Region spanning Europe and Asia
    >>> extent = (30.0, 80.0, 40.0, 60.0)
    >>> get_hydrosheds_continent_from_extent(extent, return_all=True)
    ['eu', 'as']  # or ['as', 'eu'] depending on overlap

    >>> # Dateline crossing (Pacific)
    >>> extent = (170.0, -170.0, -20.0, 10.0)  # wraps around dateline
    >>> get_hydrosheds_continent_from_extent(extent)
    'au'

    See Also
    --------
    HydroSHEDS : https://www.hydrosheds.org/
    """
    # Validate inputs
    if aExtent is None:
        return None

    if not isinstance(aExtent, (tuple, list)) or len(aExtent) != 4:
        raise ValueError(
            f"Extent must be a 4-element tuple/list (min_lon, max_lon, min_lat, max_lat). "
            f"Got {aExtent}"
        )

    if not 0.0 <= threshold_fraction <= 1.0:
        raise ValueError(
            f"threshold_fraction must be in range [0, 1]. Got {threshold_fraction}"
        )

    min_lon, max_lon, min_lat, max_lat = aExtent

    # Validate latitude range
    if not -90.0 <= min_lat <= 90.0:
        raise ValueError(f"min_lat must be in range [-90, 90]. Got {min_lat}")
    if not -90.0 <= max_lat <= 90.0:
        raise ValueError(f"max_lat must be in range [-90, 90]. Got {max_lat}")
    if min_lat > max_lat:
        raise ValueError(f"min_lat ({min_lat}) must be <= max_lat ({max_lat})")

    # Normalize longitudes to [-180, 180]
    def _norm_lon(lon: float) -> float:
        """Normalize longitude to [-180, 180] range."""
        while lon < -180.0:
            lon += 360.0
        while lon > 180.0:
            lon -= 360.0
        return lon

    min_lon = _norm_lon(min_lon)
    max_lon = _norm_lon(max_lon)

    # Handle dateline crossing by creating one or two boxes
    # Box format: (minx, miny, maxx, maxy)
    if min_lon <= max_lon:
        # Normal case: doesn't cross dateline
        input_boxes = [(min_lon, min_lat, max_lon, max_lat)]
    else:
        # Crosses dateline: split into western and eastern boxes
        input_boxes = [
            (min_lon, min_lat, 180.0, max_lat),  # Western box
            (-180.0, min_lat, max_lon, max_lat),  # Eastern box
        ]

    # Calculate total input area
    input_area = sum(_box_area(b) for b in input_boxes)
    if input_area == 0.0:
        return None

    # HydroSHEDS continent abbreviations
    continent_abbrev = {
        "Africa": "af",
        "Arctic": "ar",
        "Asia": "as",
        "Australia": "au",
        "Europe": "eu",
        "Greenland": "gr",
        "North America": "na",
        "South America": "sa",
        "Siberia": "si",
    }

    # Coarse continent bounding boxes
    # Format: (min_lon, max_lon, min_lat, max_lat)
    continent_boxes = {
        "Africa": (-20.0, 52.0, -40.0, 37.5),
        "Arctic": (-180.0, -60.0, 51.0, 90.0),
        "Asia": (55.0, 145.0, 0.0, 55.0),
        "Australia": (110.0, 180.0, -50.0, -10.0),
        "Europe": (-20.0, 75.0, 12.0, 82.0),
        "Greenland": (-73.0, -12.0, 59.0, 83.0),
        "North America": (-140.0, -50.0, 5.0, 63.0),
        "South America": (-90.0, -30.0, -60.0, 15.0),
        "Siberia": (60.0, 180.0, 45.0, 81.0),
    }

    # Convert continent boxes to standard format (minx, miny, maxx, maxy)
    cont_boxes_conv = {
        name: (box[0], box[2], box[1], box[3]) for name, box in continent_boxes.items()
    }

    # Calculate overlap scores for each continent
    scores = {}
    for cont_name, cont_box in cont_boxes_conv.items():
        inter_area = 0.0
        for input_box in input_boxes:
            inter_area += _intersection_area(input_box, cont_box)
        if inter_area > 0.0:
            scores[cont_name] = inter_area

    if not scores:
        return None

    # Sort continents by overlap area (descending)
    sorted_scores = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)

    if return_all:
        # Return all continents meeting threshold
        min_area = threshold_fraction * input_area
        significant_continents = [
            continent_abbrev.get(name, name)
            for name, area in sorted_scores
            if area >= min_area
        ]
        return significant_continents if significant_continents else None
    else:
        # Return only the continent with largest overlap
        best_continent = sorted_scores[0][0]
        return continent_abbrev.get(best_continent, best_continent)
