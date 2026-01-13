"""
Get output geometry type for OGR vector operations.

This module provides functionality to map input geometry types to appropriate
output geometry types, handling multi-geometry types and 2.5D/3D variants.
"""

from osgeo import ogr

# Geometry type mapping for common conversions
GEOMETRY_TYPE_MAPPING = {
    ogr.wkbPoint: ogr.wkbPoint,
    ogr.wkbLineString: ogr.wkbLineString,
    ogr.wkbLineString25D: ogr.wkbLineString,
    ogr.wkbPolygon: ogr.wkbPolygon,
    ogr.wkbMultiPoint: ogr.wkbPoint,
    ogr.wkbMultiLineString: ogr.wkbLineString,
    ogr.wkbMultiLineString25D: ogr.wkbLineString,
    ogr.wkbMultiPolygon: ogr.wkbPolygon,
}


def get_output_geometry_type(input_geom_type: int) -> int:
    """
    Map input geometry type to appropriate output geometry type.

    This function normalizes geometry types by:
    - Converting multi-geometry types to their single equivalents
    - Converting 2.5D/3D variants to their 2D equivalents
    - Preserving standard geometry types as-is

    Args:
        input_geom_type: OGR geometry type code (e.g., ogr.wkbLineString)

    Returns:
        int: Output geometry type code suitable for layer creation

    Raises:
        ValueError: If geometry type is not supported

    Example:
        >>> from osgeo import ogr
        >>> get_output_geometry_type(ogr.wkbMultiLineString)
        2  # ogr.wkbLineString
        >>> get_output_geometry_type(ogr.wkbLineString25D)
        2  # ogr.wkbLineString

    Note:
        Common mappings:
        - wkbMultiPoint → wkbPoint
        - wkbMultiLineString → wkbLineString
        - wkbMultiPolygon → wkbPolygon
        - wkbLineString25D → wkbLineString
        - wkbMultiLineString25D → wkbLineString
    """
    if input_geom_type in GEOMETRY_TYPE_MAPPING:
        return GEOMETRY_TYPE_MAPPING[input_geom_type]
    else:
        geom_name = ogr.GeometryTypeToName(input_geom_type)
        raise ValueError(f"Unsupported geometry type: {geom_name} ({input_geom_type})")
