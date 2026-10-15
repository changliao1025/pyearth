"""Check whether a spatial reference matches the WGS84 geographic system."""

from typing import Union, Optional
from osgeo import osr

osr.UseExceptions()

__all__ = ["is_wgs84_projection"]


def _get_wgs84_spatial_reference() -> osr.SpatialReference:
    """Create and return a cached WGS84 spatial reference with traditional GIS axis order."""
    sr = osr.SpatialReference()
    sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    sr.ImportFromEPSG(4326)
    return sr


_WGS84_SR: Optional[osr.SpatialReference] = None


def is_wgs84_projection(projection: Union[str, osr.SpatialReference]) -> bool:
    """Return True if the provided projection matches WGS84 (EPSG:4326).

    Parameters
    ----------
    projection : str or osr.SpatialReference
        Spatial reference to check. Can be a WKT string, PROJ string,
        EPSG identifier (e.g. 'EPSG:4326'), or an osr.SpatialReference instance.

    Returns
    -------
    bool
        True if the projection corresponds to EPSG:4326 WGS84, False otherwise.
    """
    if not projection:
        return False

    global _WGS84_SR
    if _WGS84_SR is None:
        _WGS84_SR = _get_wgs84_spatial_reference()

    if isinstance(projection, osr.SpatialReference):
        target_sr = projection.Clone()
        target_sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        return bool(target_sr.IsSame(_WGS84_SR))

    if not isinstance(projection, str):
        return False

    proj_str = projection.strip()
    if not proj_str:
        return False

    # Check common EPSG shorthand
    proj_upper = proj_str.upper()
    if proj_upper in ("EPSG:4326", "4326"):
        return True

    target_sr = osr.SpatialReference()
    target_sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    # Try importing as WKT first
    try:
        if target_sr.ImportFromWkt(proj_str) == 0:
            if target_sr.IsSame(_WGS84_SR) or (
                target_sr.IsGeographic()
                and target_sr.GetAuthorityCode(None) == "4326"
            ):
                return True
    except Exception:
        pass

    # Try importing as PROJ string
    try:
        if target_sr.ImportFromProj4(proj_str) == 0:
            if target_sr.IsSame(_WGS84_SR):
                return True
    except Exception:
        pass

    # Try importing from user input (e.g. EPSG:4326, urn, etc.)
    try:
        if target_sr.SetFromUserInput(proj_str) == 0:
            if target_sr.IsSame(_WGS84_SR) or (
                target_sr.IsGeographic()
                and target_sr.GetAuthorityCode(None) == "4326"
            ):
                return True
    except Exception:
        pass

    return False

