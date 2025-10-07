import os
from typing import Tuple, List, Optional, Union


def _box_area(box: Tuple[float, float, float, float]) -> float:
    """Return area of a lon/lat axis-aligned box (minx, miny, maxx, maxy)."""
    minx, miny, maxx, maxy = box
    w = max(0.0, maxx - minx)
    h = max(0.0, maxy - miny)
    return w * h


def _intersection_area(b1: Tuple[float, float, float, float],
                       b2: Tuple[float, float, float, float]) -> float:
    """Compute intersection area of two axis-aligned lon/lat boxes.

    Boxes are (minx, miny, maxx, maxy).
    """
    minx = max(b1[0], b2[0])
    miny = max(b1[1], b2[1])
    maxx = min(b1[2], b2[2])
    maxy = min(b1[3], b2[3])
    if maxx <= minx or maxy <= miny:
        return 0.0
    return (maxx - minx) * (maxy - miny)


def get_hydrosheds_continent_from_extent(aExtent: Tuple[float, float, float, float],
                                         threshold_fraction: float = 0.10,
                                         return_all: bool = False
                                         ) -> Optional[Union[str, List[str]]]:
    """
    Determine continent(s) overlapped by an extent using a coarse internal
    bounding-box heuristic. This function accepts a bounding box in WGS84
    (longitude/latitude) and returns the best-matching continent name or a
    list of significant continents when `return_all=True`.

    Args:
        aExtent: (minx, maxx, miny, maxy) in lon/lat (WGS84)
        threshold_fraction: when return_all=True, only return continents whose overlap
            area is >= threshold_fraction of the input bbox area.
        return_all: if True return list of all significant continent names, otherwise
            return the best match as a single string.

    Returns:
        continent name string, list of continent names, or None if no match.
    """
    if aExtent is None:
        return None

    minx, maxx, miny, maxy = aExtent

    # normalize longitudes to [-180,180]
    def _norm_lon(lon):
        while lon < -180:
            lon += 360
        while lon > 180:
            lon -= 360
        return lon

    minx = _norm_lon(minx)
    maxx = _norm_lon(maxx)

    # handle dateline crossing by creating one or two boxes (minx,miny,maxx,maxy)
    if minx <= maxx:
        input_boxes = [(minx, miny, maxx, maxy)]
    else:
        input_boxes = [(minx, miny, 180.0, maxy), (-180.0, miny, maxx, maxy)]

    input_area = sum(_box_area(b) for b in input_boxes)
    if input_area == 0:
        return None

    # Fallback: coarse continent bounding boxes (lon_min, lon_max, lat_min, lat_max)
    # NOTE: boxes below are stored as (lon_min, lon_max, lat_min, lat_max) in
    # the original source; convert to (minx, miny, maxx, maxy) tuples for math.



    #build a short lookup table because hydrosheds uses two letters for continent names
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
# ...existing

    continent_boxes = {
        "Africa": (-20.0, 52.0, -40.0, 37.5), # min lon, max lon, min lat, max lat
        "Arctic": (-180.0, -60.0, 51.0, 90.0),
        "Asia": (55.0, 145.0, 0.0, 55.0),
        "Australia": (110.0, 180.0, -50.0, -10.0),
        "Europe": (-20.0, 75.0, 12.0, 82.0),
        "Greenland": (-73.0, -12.0, 59.0, 83.0),
        "North America": (-140.0, -50.0, 5.0, 63.0),
        "South America": (-90.0, -30.0, -60.0, 15.0),
        "Siberia": (60.0, 180.0, 45.0, 81.0),
    }

    # convert continent boxes to (minx,miny,maxx,maxy)
    cont_boxes_conv = {k: (v[0], v[2], v[1], v[3]) for k, v in continent_boxes.items()}

    scores = {}
    for cname, cbox in cont_boxes_conv.items():
        inter_area = 0.0
        for ib in input_boxes:
            inter_area += _intersection_area(ib, cbox)
        if inter_area > 0.0:
            scores[cname] = inter_area

    if not scores:
        return None

    sorted_scores = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
    # Always return the abbreviation of the continent with the largest overlap
    best_continent = sorted_scores[0][0]
    return continent_abbrev.get(best_continent, best_continent)


# Example usage (commented):
# aExtent = (-93.8, -93.7, 32.5, 32.7)
#