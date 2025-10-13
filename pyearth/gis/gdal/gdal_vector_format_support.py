import os
def get_supported_formats():
    """
    Return a dictionary of supported vector formats
    """
    return {
        '.geojson': 'GeoJSON',
        '.json': 'GeoJSON',
        '.shp': 'ESRI Shapefile',
        '.gpkg': 'GPKG (GeoPackage)',
        '.kml': 'KML (Google Earth)',
        '.gml': 'GML (Geography Markup Language)',
        '.csv': 'CSV (with geometry column)',
        '.parquet': 'Parquet',
        '.geoparquet': 'Parquet'
    }

def print_supported_formats():
    """
    Print all supported vector formats
    """
    print("Supported vector formats:")
    formats = get_supported_formats()
    for ext, desc in formats.items():
        print(f"  {ext}: {desc}")

def get_format_from_extension(sFilename_in):
    """
    Determine the OGR format string from file extension.

    Args:
        filename: Input filename with extension

    Returns:
        Format string for OGR driver
    """
    _, ext = os.path.splitext(sFilename_in.lower())

    format_map = {
        '.geojson': 'GeoJSON',
        '.json': 'GeoJSON',
        '.shp': 'ESRI Shapefile',
        '.gpkg': 'GPKG',
        '.kml': 'KML',
        '.gml': 'GML',
        '.sqlite': 'SQLite',
        '.parquet': 'Parquet',
        '.geoparquet': 'Parquet'
    }

    return format_map.get(ext, 'GeoJSON')