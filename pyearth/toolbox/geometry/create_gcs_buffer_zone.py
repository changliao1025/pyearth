"""
Geographic Coordinate System (GCS) Buffer Zone Creation

This module provides functionality for creating buffer zones around geographic
features (points, polylines, polygons) in geographic coordinate systems (WGS84).
Buffer zones are calculated using geodesic distances on the WGS84 ellipsoid,
ensuring accurate areas regardless of latitude.

Main Functions
--------------
create_point_buffer_zone : Create circular buffer around a point in GCS
create_polyline_buffer_zone : Create buffer corridor along a polyline in GCS
create_buffer_zone_polygon_file : Batch process polygons to create buffer zones

Key Features
------------
- Geodesic buffer calculation (accurate at any latitude)
- WKT (Well-Known Text) geometry input/output
- Support for points, polylines, and polygons
- Batch processing with area-based filtering
- Multi-format vector file support (Shapefile, GeoJSON, etc.)
- Progress tracking for large datasets
- Professional error handling and validation

Use Cases
---------
1. **Proximity Analysis**: Create zones around points of interest
2. **Corridor Mapping**: Buffer transportation routes or rivers
3. **Environmental Studies**: Define impact zones around facilities
4. **Urban Planning**: Analyze accessibility within walking distance
5. **Conservation**: Create protection zones around habitats
6. **Risk Assessment**: Model hazard zones around critical infrastructure
7. **Spatial Analysis**: Generate service areas or catchment zones

Technical Details
-----------------
Buffer calculations are performed on the WGS84 ellipsoid using geodesic
distances, which accounts for Earth's curvature. This provides accurate
results across different latitudes, unlike simple Euclidean buffers.

For point buffers:
- Creates approximate circle with specified radius
- Vertices calculated at regular angular intervals
- More vertices at higher latitudes for accuracy

For polyline buffers:
- Creates corridor with parallel offset lines
- Handles line segment endpoints with rounded caps
- Connects segments with appropriate joins

Dependencies
------------
- osgeo (GDAL/OGR): Geometry operations and file I/O
- numpy: Numerical operations
- pyearth.toolbox.mesh: point, edge, polyline, polygon classes
- pyearth.gis: Coordinate extraction and area calculation

See Also
--------
- pypoint.calculate_buffer_zone: Point buffer implementation
- pypolyline.calculate_buffer_zone: Polyline buffer implementation
- pypolygon.calculate_buffer_zone: Polygon buffer implementation
"""

import os
import sys
import logging
import numpy as np
from typing import Optional
from osgeo import ogr, gdal, osr

from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_format_from_filename
from pyearth.toolbox.mesh.point import pypoint
from pyearth.toolbox.mesh.line import pyline
from pyearth.toolbox.mesh.polyline import pypolyline
from pyearth.toolbox.mesh.polygon import pypolygon

# Configure logging
logger = logging.getLogger(__name__)


def create_point_buffer_zone(sWkt: str, dBuffer_distance_in: float) -> Optional[str]:
    """
    Create a circular buffer zone around a point in geographic coordinates.

    This function generates a geodesic buffer (approximate circle) around a point
    specified in WKT format. The buffer is calculated using geodesic distances
    on the WGS84 ellipsoid, ensuring accuracy across different latitudes.

    Parameters
    ----------
    sWkt : str
        Well-Known Text (WKT) representation of a POINT geometry.
        Must be a valid WKT string in the format: "POINT (longitude latitude)"

        Examples:
        - "POINT (-122.4194 37.7749)" for San Francisco
        - "POINT (0.0 51.5074)" for London

    dBuffer_distance_in : float
        Buffer distance in meters. Must be positive.

        - Typical values: 100-10000 meters for local analysis
        - Maximum: Limited by point class implementation
        - Minimum: Should be > 0 for meaningful results

    Returns
    -------
    str or None
        WKT representation of the buffer polygon (approximate circle).
        Returns None if:
        - Input WKT is invalid
        - Geometry is not a POINT
        - Point extraction fails
        - Buffer calculation fails

        Output format: "POLYGON ((lon1 lat1, lon2 lat2, ...))"

    Raises
    ------
    TypeError
        If sWkt is None or not a string.
        If dBuffer_distance_in is None or not numeric.
    ValueError
        If dBuffer_distance_in is not positive.

    Notes
    -----
    1. **Geodesic Calculation**: Uses WGS84 ellipsoid for accurate distances
    2. **Coordinate Order**: WKT format is "longitude latitude" (x y)
    3. **Buffer Shape**: Creates approximate circle with multiple points
    4. **Point Count**: More points used at higher latitudes for accuracy
    5. **Return Format**: Closed polygon in WKT format (first = last point)
    6. **Coordinate System**: Assumes input/output in WGS84 (EPSG:4326)
    7. **Error Handling**: Returns None on errors with logged messages
    8. **Distance Units**: Input in meters, calculation on ellipsoid
    9. **Latitude Range**: Valid for -90° to +90° latitude
    10. **Longitude Range**: Valid for -180° to +180° longitude

    Examples
    --------
    Create 1km buffer around a point in New York City:

    >>> wkt_point = "POINT (-74.0060 40.7128)"
    >>> buffer_distance = 1000  # 1 kilometer
    >>> wkt_buffer = create_point_buffer_zone(wkt_point, buffer_distance)
    >>> print(wkt_buffer[:50])
    POLYGON ((-74.0060 40.7218, -74.0070 40.7216, ...

    Create 5km buffer around the North Pole (high latitude):

    >>> wkt_pole = "POINT (0.0 89.0)"
    >>> buffer_5km = create_point_buffer_zone(wkt_pole, 5000)
    >>> # Buffer will be nearly circular despite high latitude

    Create buffer around the equator:

    >>> wkt_equator = "POINT (0.0 0.0)"
    >>> buffer_10km = create_point_buffer_zone(wkt_equator, 10000)
    >>> # Buffer calculation is most accurate near equator

    Handle invalid input:

    >>> invalid_wkt = "LINESTRING (0 0, 1 1)"
    >>> result = create_point_buffer_zone(invalid_wkt, 1000)
    Error: Input geometry must be a POINT for buffer creation.
    >>> print(result)
    None

    Create buffer for spatial analysis:

    >>> # Define point of interest (Statue of Liberty)
    >>> poi = "POINT (-74.0445 40.6892)"
    >>> # Create 500m accessibility zone
    >>> accessibility_zone = create_point_buffer_zone(poi, 500)
    >>> # Use with OGR to find features within zone

    See Also
    --------
    create_polyline_buffer_zone : Create buffer around polyline
    create_buffer_zone_polygon_file : Batch buffer creation for polygons
    pypoint.calculate_buffer_zone : Underlying buffer calculation method

    References
    ----------
    .. [1] Karney, C. F. F. (2013). Algorithms for geodesics.
           Journal of Geodesy, 87(1), 43-55.
    .. [2] OGC Simple Features Specification - Well-Known Text
           https://www.ogc.org/standards/sfa
    """
    # Validate input parameters
    if sWkt is None:
        error_msg = "WKT string cannot be None"
        logger.error(error_msg)
        raise TypeError(error_msg)

    if not isinstance(sWkt, str):
        error_msg = f"WKT must be a string, got {type(sWkt).__name__}"
        logger.error(error_msg)
        raise TypeError(error_msg)

    if sWkt.strip() == "":
        error_msg = "WKT string cannot be empty"
        logger.error(error_msg)
        raise ValueError(error_msg)

    if dBuffer_distance_in is None:
        error_msg = "Buffer distance cannot be None"
        logger.error(error_msg)
        raise TypeError(error_msg)

    try:
        dBuffer_distance_in = float(dBuffer_distance_in)
    except (ValueError, TypeError):
        error_msg = (
            f"Buffer distance must be numeric, got {type(dBuffer_distance_in).__name__}"
        )
        logger.error(error_msg)
        raise TypeError(error_msg)

    if dBuffer_distance_in <= 0:
        error_msg = f"Buffer distance must be positive, got {dBuffer_distance_in}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Create point geometry from WKT
    try:
        pGeometry = ogr.CreateGeometryFromWkt(sWkt)
    except Exception as e:
        error_msg = f"Failed to parse WKT string: {str(e)}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    if pGeometry is None:
        error_msg = f"Invalid WKT input for point buffer creation: {sWkt}"
        logger.error(error_msg)
        return None

    # Validate geometry type
    sGeometry_type = pGeometry.GetGeometryName()
    if sGeometry_type != "POINT":
        error_msg = (
            f"Input geometry must be a POINT for buffer creation, got {sGeometry_type}"
        )
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return None

    # Extract coordinates
    try:
        aCoords_gcs = get_geometry_coordinates(pGeometry)
    except Exception as e:
        error_msg = f"Failed to extract coordinates: {str(e)}"
        logger.error(error_msg)
        return None

    if len(aCoords_gcs) != 1:
        error_msg = f"Input geometry must be a single point for buffer creation, got {len(aCoords_gcs)} points"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return None

    # Validate coordinate values
    lon, lat = aCoords_gcs[0][0], aCoords_gcs[0][1]
    if not (-180 <= lon <= 180):
        error_msg = f"Longitude must be in range [-180, 180], got {lon}"
        logger.warning(error_msg)

    if not (-90 <= lat <= 90):
        error_msg = f"Latitude must be in range [-90, 90], got {lat}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Create point object
    point = {"dLongitude_degree": lon, "dLatitude_degree": lat}

    try:
        pPoint = pypoint(point)
    except Exception as e:
        error_msg = f"Failed to create point: {str(e)}"
        logger.error(error_msg)
        return None

    # Calculate buffer zone
    try:
        sWkt_buffer_polygon = pPoint.calculate_buffer_zone(dBuffer_distance_in)
    except Exception as e:
        error_msg = f"Failed to calculate buffer zone: {str(e)}"
        logger.error(error_msg)
        return None

    if sWkt_buffer_polygon is None:
        logger.warning("Buffer calculation returned None")

    return sWkt_buffer_polygon


def create_polyline_buffer_zone(sWkt: str, dBuffer_distance_in: float) -> Optional[str]:
    """
    Create a buffer corridor around a polyline in geographic coordinates.

    This function generates a geodesic buffer (corridor) around a polyline
    specified in WKT format. The buffer is calculated using geodesic distances
    on the WGS84 ellipsoid, creating parallel offset lines on both sides of
    the polyline with rounded caps at endpoints.

    Parameters
    ----------
    sWkt : str
        Well-Known Text (WKT) representation of a LINESTRING geometry.
        Must be a valid WKT string in the format:
        "LINESTRING (lon1 lat1, lon2 lat2, lon3 lat3, ...)"

        Examples:
        - "LINESTRING (-122.4 37.8, -122.5 37.9)" for SF route
        - "LINESTRING (0 0, 1 1, 2 0)" for multi-segment line

        Requirements:
        - At least 2 points (start and end)
        - Coordinates in decimal degrees
        - Valid longitude [-180, 180] and latitude [-90, 90]

    dBuffer_distance_in : float
        Buffer distance in meters (half-width of corridor). Must be positive.

        - Typical values: 10-5000 meters for roads/rivers
        - Creates total corridor width of 2 × dBuffer_distance_in
        - Uses geodesic distance for accuracy

    Returns
    -------
    str or None
        WKT representation of the buffer polygon (corridor).
        Returns None if:
        - Input WKT is invalid
        - Geometry is not a LINESTRING
        - Coordinate extraction fails
        - Polyline has < 2 valid points
        - Buffer calculation fails

        Output format: "POLYGON ((lon1 lat1, lon2 lat2, ...))"
        Polygon includes:
        - Parallel offset on both sides of line
        - Rounded caps at endpoints
        - Smooth joins at vertices

    Raises
    ------
    TypeError
        If sWkt is None or not a string.
        If dBuffer_distance_in is None or not numeric.
    ValueError
        If sWkt is empty string.
        If dBuffer_distance_in is not positive.
        If coordinates are outside valid ranges.

    Notes
    -----
    1. **Geodesic Calculation**: Uses WGS84 ellipsoid for accurate distances
    2. **Coordinate Order**: WKT format is "longitude latitude" (x y)
    3. **Buffer Shape**: Creates corridor with rounded caps at line endpoints
    4. **Duplicate Vertices**: Automatically removed during edge construction
    5. **Segment Processing**: Each line segment buffered individually
    6. **Join Style**: Vertices connected with appropriate geodesic curves
    7. **Cap Style**: Rounded caps approximate semicircles at endpoints
    8. **Coordinate System**: Assumes input/output in WGS84 (EPSG:4326)
    9. **Error Handling**: Returns None on errors with logged messages
    10. **Minimum Points**: Requires at least 2 distinct points for valid line

    Algorithm Details
    -----------------
    1. Parse WKT to extract LINESTRING coordinates
    2. Create pypoint objects for each point
    3. Build edge objects between consecutive vertices
    4. Skip duplicate vertices to avoid zero-length edges
    5. Construct polyline from edge list
    6. Calculate geodesic buffer using pypolyline.calculate_buffer_zone()
    7. Return buffer polygon as WKT

    Examples
    --------
    Create 100m buffer corridor along a road:

    >>> wkt_road = "LINESTRING (-122.4194 37.7749, -122.4089 37.7849)"
    >>> buffer_distance = 100  # 100 meters (50m each side)
    >>> wkt_corridor = create_polyline_buffer_zone(wkt_road, buffer_distance)
    >>> # Creates 200m wide corridor along the road

    Create buffer along a river with multiple segments:

    >>> wkt_river = "LINESTRING (-73.9 40.7, -73.8 40.75, -73.7 40.8)"
    >>> river_buffer = create_polyline_buffer_zone(wkt_river, 500)
    >>> # Creates 1km wide corridor along the river

    Create buffer for pipeline route:

    >>> pipeline = "LINESTRING (0.0 51.5, 0.1 51.6, 0.2 51.65, 0.3 51.7)"
    >>> safety_zone = create_polyline_buffer_zone(pipeline, 1000)
    >>> # Creates 2km wide safety corridor

    Handle invalid input:

    >>> invalid_wkt = "POINT (0 0)"
    >>> result = create_polyline_buffer_zone(invalid_wkt, 100)
    Error: Input geometry must be a LINESTRING for buffer creation, got POINT
    >>> print(result)
    None

    Create buffer for transportation corridor:

    >>> # Highway segment
    >>> highway = "LINESTRING (-122.5 37.8, -122.4 37.85, -122.3 37.9)"
    >>> # 50m buffer for impact analysis
    >>> impact_zone = create_polyline_buffer_zone(highway, 50)
    >>> # Total corridor width is 100m

    Process multi-segment utility line:

    >>> # Power transmission line
    >>> transmission = "LINESTRING (-74.0 40.7, -74.05 40.72, -74.1 40.74)"
    >>> # 25m maintenance corridor
    >>> maintenance_zone = create_polyline_buffer_zone(transmission, 25)

    See Also
    --------
    create_point_buffer_zone : Create buffer around point
    create_buffer_zone_polygon_file : Batch buffer creation for polygons
    pypolyline.calculate_buffer_zone : Underlying buffer calculation method
    pyedge : Edge class for line segments

    References
    ----------
    .. [1] Karney, C. F. F. (2013). Algorithms for geodesics.
           Journal of Geodesy, 87(1), 43-55.
    .. [2] OGC Simple Features Specification - Well-Known Text
           https://www.ogc.org/standards/sfa
    .. [3] PostGIS ST_Buffer documentation for buffer styles
           https://postgis.net/docs/ST_Buffer.html
    """
    # Validate input parameters
    if sWkt is None:
        error_msg = "WKT string cannot be None"
        logger.error(error_msg)
        raise TypeError(error_msg)

    if not isinstance(sWkt, str):
        error_msg = f"WKT must be a string, got {type(sWkt).__name__}"
        logger.error(error_msg)
        raise TypeError(error_msg)

    if sWkt.strip() == "":
        error_msg = "WKT string cannot be empty"
        logger.error(error_msg)
        raise ValueError(error_msg)

    if dBuffer_distance_in is None:
        error_msg = "Buffer distance cannot be None"
        logger.error(error_msg)
        raise TypeError(error_msg)

    try:
        dBuffer_distance_in = float(dBuffer_distance_in)
    except (ValueError, TypeError):
        error_msg = (
            f"Buffer distance must be numeric, got {type(dBuffer_distance_in).__name__}"
        )
        logger.error(error_msg)
        raise TypeError(error_msg)

    if dBuffer_distance_in <= 0:
        error_msg = f"Buffer distance must be positive, got {dBuffer_distance_in}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Create polyline geometry from WKT
    try:
        pGeometry = ogr.CreateGeometryFromWkt(sWkt)
    except Exception as e:
        error_msg = f"Failed to parse WKT string: {str(e)}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    if pGeometry is None:
        error_msg = f"Invalid WKT input for polyline buffer creation: {sWkt}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return None

    # Validate geometry type
    sGeometry_type = pGeometry.GetGeometryName()
    if sGeometry_type != "LINESTRING":
        error_msg = f"Input geometry must be a LINESTRING for buffer creation, got {sGeometry_type}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return None

    # Extract coordinates
    try:
        aCoords_gcs = get_geometry_coordinates(pGeometry)
    except Exception as e:
        error_msg = f"Failed to extract coordinates: {str(e)}"
        logger.error(error_msg)
        return None

    nPoint = len(aCoords_gcs)

    if nPoint < 2:
        error_msg = f"LINESTRING must have at least 2 points, got {nPoint}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return None

    # Create point objects from coordinates
    aPoint = []
    for i in range(nPoint):
        lon, lat = aCoords_gcs[i][0], aCoords_gcs[i][1]

        # Validate coordinate ranges
        if not (-180 <= lon <= 180):
            error_msg = (
                f"Longitude at point {i} must be in range [-180, 180], got {lon}"
            )
            logger.warning(error_msg)

        if not (-90 <= lat <= 90):
            error_msg = f"Latitude at point {i} must be in range [-90, 90], got {lat}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        point = {"dLongitude_degree": lon, "dLatitude_degree": lat}

        try:
            pPoint = pypoint(point)
            aPoint.append(pPoint)
        except Exception as e:
            error_msg = f"Failed to create point at point {i}: {str(e)}"
            logger.error(error_msg)
            return None

    nPoint = len(aPoint)

    # Build edges between consecutive points, skipping duplicates
    aEdge = []
    for i in range(nPoint - 1):
        if aPoint[i] != aPoint[i + 1]:
            try:
                pEdge = pyline(aPoint[i], aPoint[i + 1])
                aEdge.append(pEdge)
            except Exception as e:
                error_msg = (
                    f"Failed to create edge between points {i} and {i+1}: {str(e)}"
                )
                logger.warning(error_msg)
                # Continue processing other edges
        else:
            logger.debug(f"Skipping duplicate point at position {i}")

    if len(aEdge) == 0:
        error_msg = "No valid edges created (all points may be duplicates)"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return None

    # Create polyline from edges
    try:
        pPolyline = pypolyline(aEdge)
    except Exception as e:
        error_msg = f"Failed to create polyline: {str(e)}"
        logger.error(error_msg)
        return None

    # Calculate buffer zone
    try:
        sWkt_buffer_polygon, _, _ = pPolyline.calculate_buffer_zone_polygon(
            dBuffer_distance_in
        )
    except Exception as e:
        error_msg = f"Failed to calculate buffer zone: {str(e)}"
        logger.error(error_msg)
        return None

    if sWkt_buffer_polygon is None:
        logger.warning("Buffer calculation returned None")

    return sWkt_buffer_polygon


def create_buffer_zone_polygon_file(
    sFilename_polygon_in: str,
    sFilename_polygon_out: str,
    dThreshold_in: Optional[float] = 1.0e9,
    dBuffer_distance_in: float = 5000.0,
    verbose: bool = True,
) -> None:
    """
    Create geodesic buffer zones for all polygons in a vector file with batch processing.

    This function reads polygons from an input vector file, creates geodesic buffer
    zones around each polygon using WGS84 ellipsoid calculations, and writes the
    results to an output vector file. Supports area-based filtering to exclude
    small polygons and provides progress tracking for large datasets.

    Parameters
    ----------
    sFilename_polygon_in : str
        Path to input vector file containing polygons.

        Supported formats:
        - Shapefile (.shp)
        - GeoJSON (.geojson, .json)
        - GeoPackage (.gpkg)
        - GML (.gml)
        - KML (.kml)
        - And other OGR-supported formats

        Must contain POLYGON or MULTIPOLYGON geometries.

    sFilename_polygon_out : str
        Path to output vector file for buffer zones.

        Format determined by file extension:
        - .shp → ESRI Shapefile
        - .geojson/.json → GeoJSON
        - .gpkg → GeoPackage
        - .gml → GML
        - .kml → KML

        If file exists, it will be overwritten.

    dThreshold_in : float or None, optional
        Minimum polygon area threshold in square meters.
        Polygons with area < threshold are excluded from buffering.

        Default: 1.0E9 (1,000,000,000 m² ≈ 1,000 km²)

        Special values:
        - None: Process all polygons regardless of size
        - 0: Process all polygons (same as None)
        - Positive value: Filter by area threshold

    dBuffer_distance_in : float, optional
        Buffer distance in meters for geodesic buffer calculation.

        Default: 5000.0 (5 kilometers)

        Typical ranges:
        - Urban analysis: 100-2000 meters
        - Regional planning: 1000-10000 meters
        - Large-scale studies: 10000+ meters

    verbose : bool, optional
        Enable progress reporting and diagnostic messages.

        Default: True

        When True, prints:
        - Format detection
        - Processing progress
        - Feature counts
        - Warnings for skipped features

    Returns
    -------
    None
        Results written to sFilename_polygon_out.

    Raises
    ------
    TypeError
        If input parameters have incorrect types.
    ValueError
        If file paths are invalid or buffer distance is non-positive.
    FileNotFoundError
        If input file does not exist.
    RuntimeError
        If OGR driver is not available or file operations fail.

    Notes
    -----
    1. **Geodesic Buffers**: Uses WGS84 ellipsoid for accurate distances
    2. **Area Calculation**: Uses geodesic area algorithm (iFlag_algorithm=2)
    3. **Multi-Part Support**: Processes each part of MULTIPOLYGON separately
    4. **Degenerate Polygons**: Skips polygons with < 3 vertices
    5. **Duplicate Vertices**: Automatically removed during edge construction
    6. **Coordinate System**: Assumes input in WGS84 (EPSG:4326)
    7. **Progress Tracking**: Reports every 100 features
    8. **Output Fields**:
       - id: Sequential feature identifier
       - orig_area: Original polygon area in m²
       - buffer_dist: Buffer distance applied in meters
    9. **Memory Management**: Properly releases GDAL/OGR resources
    10. **Error Handling**: Continues processing on per-feature failures

    Output Schema
    -------------
    The output file contains these fields:
    - **id** (Integer): Sequential feature ID (1, 2, 3, ...)
    - **orig_area** (Real): Area of original polygon in m²
    - **buffer_dist** (Real): Buffer distance applied in meters
    - **geometry** (Polygon): Buffer zone polygon

    Processing Workflow
    -------------------
    1. Validate input file exists
    2. Detect input/output formats from extensions
    3. Open input dataset and layer
    4. Create output dataset with appropriate driver
    5. For each input feature:
       a. Extract polygon geometry
       b. Calculate geodesic area
       c. Filter by area threshold
       d. Create pypoint/line objects
       e. Calculate geodesic buffer
       f. Write buffer polygon to output
    6. Report processing statistics
    7. Clean up resources

    Examples
    --------
    Process all large polygons with default 5km buffer:

    >>> create_buffer_zone_polygon_file(
    ...     'watersheds.shp',
    ...     'watershed_buffers.geojson'
    ... )
    Processing 1523 features...
    Processed 1523/1523 features, created 234 buffers

    Create 1km buffers around medium-sized protected areas:

    >>> create_buffer_zone_polygon_file(
    ...     'protected_areas.gpkg',
    ...     'protection_buffers.shp',
    ...     dThreshold_in=1.0E6,  # 1 km²
    ...     dBuffer_distance_in=1000.0  # 1 km
    ... )

    Process all polygons regardless of size (silent mode):

    >>> create_buffer_zone_polygon_file(
    ...     'parcels.geojson',
    ...     'parcel_buffers.gpkg',
    ...     dThreshold_in=None,  # No filtering
    ...     dBuffer_distance_in=100.0,  # 100m
    ...     verbose=False
    ... )

    Large-scale regional analysis:

    >>> # Create 10km buffers around cities > 100 km²
    >>> create_buffer_zone_polygon_file(
    ...     'cities.shp',
    ...     'city_influence_zones.geojson',
    ...     dThreshold_in=1.0E8,  # 100 km²
    ...     dBuffer_distance_in=10000.0  # 10 km
    ... )
    Input format: .shp
    Output format: .geojson -> GeoJSON
    Buffer distance: 10000.0 meters
    Area threshold: 100000000.0 m²
    Processing 456 features...
    Processed 456/456 features, created 89 buffers

    See Also
    --------
    create_point_buffer_zone : Create buffer around single point
    create_polyline_buffer_zone : Create buffer around polyline
    pypolygon.calculate_buffer_zone : Underlying polygon buffer method
    calculate_polygon_area : Geodesic area calculation

    References
    ----------
    .. [1] Karney, C. F. F. (2013). Algorithms for geodesics.
           Journal of Geodesy, 87(1), 43-55.
    .. [2] GDAL/OGR Vector Data Model
           https://gdal.org/user/vector_data_model.html
    """
    # Validate input parameters
    if sFilename_polygon_in is None or not isinstance(sFilename_polygon_in, str):
        error_msg = f"Input filename must be a string, got {type(sFilename_polygon_in).__name__}"
        logger.error(error_msg)
        raise TypeError(error_msg)

    if sFilename_polygon_out is None or not isinstance(sFilename_polygon_out, str):
        error_msg = f"Output filename must be a string, got {type(sFilename_polygon_out).__name__}"
        logger.error(error_msg)
        raise TypeError(error_msg)

    if dBuffer_distance_in <= 0:
        error_msg = f"Buffer distance must be positive, got {dBuffer_distance_in}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    if dThreshold_in is not None and dThreshold_in < 0:
        error_msg = f"Area threshold must be non-negative or None, got {dThreshold_in}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Check input file exists
    if not os.path.exists(sFilename_polygon_in):
        error_msg = f"Input file does not exist: {sFilename_polygon_in}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        raise FileNotFoundError(error_msg)

    # Auto-detect formats from file extensions
    input_ext = os.path.splitext(sFilename_polygon_in)[1].lower()
    output_ext = os.path.splitext(sFilename_polygon_out)[1].lower()

    # Get output driver name from extension
    try:
        output_driver_name = get_vector_format_from_filename(sFilename_polygon_out)
    except ValueError as e:
        # Fall back to GeoJSON if format is not supported
        logger.warning(
            f"Unsupported output format {output_ext}, using GeoJSON: {str(e)}"
        )
        output_driver_name = "GeoJSON"
        output_ext = ".geojson"

    if verbose:
        print(f"Input format: {input_ext}")
        print(f"Output format: {output_ext} -> {output_driver_name}")
        print(f"Buffer distance: {dBuffer_distance_in} meters")
        if dThreshold_in is not None and dThreshold_in > 0:
            print(f"Area threshold: {dThreshold_in} m²")
        else:
            print("Area threshold: None (processing all polygons)")

    # Open input dataset (auto-detect format)
    try:
        pDataSource = ogr.Open(sFilename_polygon_in, 0)
    except Exception as e:
        error_msg = f"Failed to open input file: {str(e)}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        raise RuntimeError(error_msg)

    if pDataSource is None:
        error_msg = f"Could not open input file: {sFilename_polygon_in}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        raise RuntimeError(error_msg)

    # Get layer
    pLayer = pDataSource.GetLayer()
    if pLayer is None:
        error_msg = "No layer found in input file"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        pDataSource = None
        raise RuntimeError(error_msg)

    # Get total feature count for progress tracking
    nTotal_features = pLayer.GetFeatureCount()
    if verbose:
        print(f"Processing {nTotal_features} features...")

    # Get output driver
    pDriver_out = ogr.GetDriverByName(output_driver_name)
    if pDriver_out is None:
        error_msg = f"Driver {output_driver_name} not available!"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        print(f"Hint: Check GDAL installation or try a different output format")
        pDataSource = None
        raise RuntimeError(error_msg)

    # Prepare output (overwrite if exists)
    if os.path.exists(sFilename_polygon_out):
        if verbose:
            print(f"Removing existing output file: {sFilename_polygon_out}")
        try:
            pDriver_out.DeleteDataSource(sFilename_polygon_out)
        except Exception as e:
            logger.warning(f"Could not delete existing file: {str(e)}")

    # Create output datasource
    try:
        pOutDataSource = pDriver_out.CreateDataSource(sFilename_polygon_out)
    except Exception as e:
        error_msg = f"Failed to create output file: {str(e)}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        pDataSource = None
        raise RuntimeError(error_msg)

    if pOutDataSource is None:
        error_msg = f"Could not create output file: {sFilename_polygon_out}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        pDataSource = None
        raise RuntimeError(error_msg)

    # Handle layer creation based on output format
    try:
        if output_driver_name == "ESRI Shapefile":
            # Shapefile requires simpler layer name
            layer_name = os.path.splitext(os.path.basename(sFilename_polygon_out))[0]
            pOutLayer = pOutDataSource.CreateLayer(layer_name, geom_type=ogr.wkbPolygon)
        else:
            pOutLayer = pOutDataSource.CreateLayer("buffer", geom_type=ogr.wkbPolygon)
    except Exception as e:
        error_msg = f"Failed to create output layer: {str(e)}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        pDataSource = None
        pOutDataSource = None
        raise RuntimeError(error_msg)

    if pOutLayer is None:
        error_msg = "Could not create output layer"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        pDataSource = None
        pOutDataSource = None
        raise RuntimeError(error_msg)

    # Add fields for statistics
    try:
        pFieldDefn = ogr.FieldDefn("id", ogr.OFTInteger)
        pOutLayer.CreateField(pFieldDefn)
        pFieldDefn = ogr.FieldDefn("orig_area", ogr.OFTReal)
        pOutLayer.CreateField(pFieldDefn)
        pFieldDefn = ogr.FieldDefn("buffer_dist", ogr.OFTReal)
        pOutLayer.CreateField(pFieldDefn)
    except Exception as e:
        error_msg = f"Failed to create output fields: {str(e)}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        pDataSource = None
        pOutDataSource = None
        raise RuntimeError(error_msg)

    # Process features
    pFeature = pLayer.GetNextFeature()
    nProcessed = 0
    nBuffered = 0
    nSkipped_area = 0
    nSkipped_degenerate = 0
    nSkipped_type = 0
    nFailed = 0
    feature_id = 1

    while pFeature:
        nProcessed += 1
        if verbose and nProcessed % 100 == 0:
            print(
                f"Processed {nProcessed}/{nTotal_features} features, created {nBuffered} buffers"
            )

        pGeometry = pFeature.GetGeometryRef()
        if pGeometry is None:
            logger.warning(f"Feature {nProcessed} has no geometry")
            pFeature = pLayer.GetNextFeature()
            nSkipped_degenerate += 1
            continue

        sGeometry_type = pGeometry.GetGeometryName()

        # Handle POLYGON and MULTIPOLYGON
        if sGeometry_type == "MULTIPOLYGON":
            try:
                aaCoords_gcs = get_geometry_coordinates(pGeometry)
                nPart = len(aaCoords_gcs)
            except Exception as e:
                logger.warning(
                    f"Failed to extract coordinates from feature {nProcessed}: {str(e)}"
                )
                pFeature = pLayer.GetNextFeature()
                nFailed += 1
                continue
        elif sGeometry_type == "POLYGON":
            try:
                aCoords_gcs = get_geometry_coordinates(pGeometry)
                aaCoords_gcs = [aCoords_gcs]
                nPart = 1
            except Exception as e:
                logger.warning(
                    f"Failed to extract coordinates from feature {nProcessed}: {str(e)}"
                )
                pFeature = pLayer.GetNextFeature()
                nFailed += 1
                continue
        else:
            if verbose:
                logger.debug(
                    f"Skipping non-polygon geometry type {sGeometry_type} in feature {nProcessed}"
                )
            pFeature = pLayer.GetNextFeature()
            nSkipped_type += 1
            continue

        # Process each polygon part
        for i in range(nPart):
            aCoords_gcs = aaCoords_gcs[i]

            # Skip degenerate polygons
            if len(aCoords_gcs) < 3:
                logger.debug(
                    f"Skipping degenerate polygon with {len(aCoords_gcs)} points"
                )
                nSkipped_degenerate += 1
                continue

            aCoords_gcs = np.array(aCoords_gcs)

            # Calculate the geodesic area of the polygon
            try:
                dArea = calculate_polygon_area(
                    aCoords_gcs[:, 0],
                    aCoords_gcs[:, 1],
                    iFlag_algorithm=2,  # Geodesic algorithm
                )
            except Exception as e:
                logger.warning(
                    f"Failed to calculate area for feature {nProcessed} part {i}: {str(e)}"
                )
                nFailed += 1
                continue

            # Filter by area threshold
            if (
                dThreshold_in is not None
                and dThreshold_in > 0
                and dArea < dThreshold_in
            ):
                logger.debug(
                    f"Skipping polygon with area {dArea:.2f} m² < threshold {dThreshold_in} m²"
                )
                nSkipped_area += 1
                continue

            # Remove the last point which is the same as the first point
            nPoint = len(aCoords_gcs) - 1
            aPoint = []

            # Create point objects
            for j in range(nPoint):
                point = {
                    "dLongitude_degree": aCoords_gcs[j, 0],
                    "dLatitude_degree": aCoords_gcs[j, 1],
                }
                try:
                    pPoint = pypoint(point)
                    aPoint.append(pPoint)
                except Exception as e:
                    logger.warning(
                        f"Failed to create point {j} for feature {nProcessed}: {str(e)}"
                    )

            if len(aPoint) < 3:
                logger.debug(f"Not enough valid points ({len(aPoint)}) for polygon")
                nSkipped_degenerate += 1
                continue

            nPoint2 = len(aPoint)
            aEdge = []

            # Build edges, skipping duplicates
            for j in range(nPoint2 - 1):
                if aPoint[j] != aPoint[j + 1]:
                    try:
                        pEdge = pyline(aPoint[j], aPoint[j + 1])
                        aEdge.append(pEdge)
                    except Exception as e:
                        logger.warning(
                            f"Failed to create edge {j} for feature {nProcessed}: {str(e)}"
                        )

            # Close the polygon
            if nPoint2 > 0:
                try:
                    pEdge = pyline(aPoint[nPoint2 - 1], aPoint[0])
                    aEdge.append(pEdge)
                except Exception as e:
                    logger.warning(
                        f"Failed to close polygon for feature {nProcessed}: {str(e)}"
                    )

                if len(aEdge) < 3:
                    logger.debug(f"Not enough valid edges ({len(aEdge)}) for polygon")
                    nSkipped_degenerate += 1
                    continue

                try:
                    # Create polygon and calculate buffer
                    pPolygon = pypolygon(aEdge)
                    sWkt_buffer = pPolygon.calculate_buffer_zone(dBuffer_distance_in)

                    if sWkt_buffer:
                        buffer_geom = ogr.CreateGeometryFromWkt(sWkt_buffer)
                        if buffer_geom:
                            # Create output feature
                            outFeature = ogr.Feature(pOutLayer.GetLayerDefn())
                            outFeature.SetGeometry(buffer_geom)
                            outFeature.SetField("id", feature_id)
                            outFeature.SetField("orig_area", float(dArea))
                            outFeature.SetField(
                                "buffer_dist", float(dBuffer_distance_in)
                            )
                            pOutLayer.CreateFeature(outFeature)
                            outFeature = None  # Free feature
                            nBuffered += 1
                            feature_id += 1
                        else:
                            logger.warning(
                                f"Failed to create geometry from WKT for feature {nProcessed}"
                            )
                            nFailed += 1
                    else:
                        logger.warning(
                            f"Buffer calculation returned None for feature {nProcessed}"
                        )
                        nFailed += 1
                except Exception as e:
                    if verbose:
                        logger.warning(
                            f"Failed to create buffer for polygon part {i} in feature {nProcessed}: {str(e)}"
                        )
                    nFailed += 1

        pFeature = pLayer.GetNextFeature()

    # Cleanup
    pDataSource = None
    pOutDataSource = None

    # Report final statistics
    if verbose:
        print(f'\n{"="*60}')
        print(f"Processing complete!")
        print(f'{"="*60}')
        print(f"Total features processed: {nProcessed}")
        print(f"Buffer zones created: {nBuffered}")
        print(f"Skipped (below area threshold): {nSkipped_area}")
        print(f"Skipped (degenerate geometry): {nSkipped_degenerate}")
        print(f"Skipped (wrong geometry type): {nSkipped_type}")
        print(f"Failed (errors): {nFailed}")
        print(f"Output saved to: {sFilename_polygon_out}")
        print(f'{"="*60}')

    logger.info(
        f"Buffer zone processing complete: {nBuffered} buffers created from {nProcessed} features"
    )

    return
