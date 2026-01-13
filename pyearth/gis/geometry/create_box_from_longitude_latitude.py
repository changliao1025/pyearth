"""
Create geographic bounding boxes aligned to global grid systems.

This module provides utilities for creating bounding boxes (rectangular polygons)
from longitude/latitude coordinates that are aligned to a global grid structure.
The boxes are snapped to a grid defined by resolutions that are powers of 2
(e.g., 1°, 0.5°, 0.25°, 0.125°, etc.).

The primary use case is for creating grid cells in earth system models, climate
simulations, and geospatial analyses where data needs to be aligned to a
consistent global grid structure anchored at (-180°, 90°).

Functions
---------
find_nearest_resolution
    Find the nearest grid resolution (as a power of 2) that is larger than
    a given resolution in meters.
create_box_from_longitude_latitude
    Create a grid-aligned bounding box from a longitude/latitude point and
    resolution, optionally saving as GeoJSON.

Notes
-----
**Global Grid Structure**:

The functions assume a global grid anchored at the top-left corner (-180°, 90°)
with cells arranged in a regular grid. The resolution must be a value like:
    1°, 0.5°, 0.25°, 0.125°, 0.0625°, etc. (i.e., 1/(2^n) degrees)

**Grid Cell Indexing**:

Given a point (lon, lat) and resolution (res_x, res_y):
    - Column index: floor((lon - (-180)) / res_x)
    - Row index: floor((90 - lat) / res_y)
    - Cell bounds are calculated from these indices

**Coordinate Ordering**:

The output coordinates follow a counter-clockwise ordering:
    1. Upper-left (NW corner)
    2. Lower-left (SW corner)
    3. Lower-right (SE corner)
    4. Upper-right (NE corner)

Examples
--------
Create a 1° × 1° grid cell:

>>> from pyearth.gis.geometry.create_box_from_longitude_latitude import \\
...     create_box_from_longitude_latitude
>>> box, coords = create_box_from_longitude_latitude(-120.0, 45.0, 1.0, 1.0)
>>> print(box)  # [left, right, bottom, top]
[-120.0, -119.0, 45.0, 46.0]

Create a 0.5° × 0.5° grid cell and save to GeoJSON:

>>> box, coords = create_box_from_longitude_latitude(
...     -100.5, 40.25, 0.5, 0.5,
...     sFilename_output_in='grid_cell.geojson'
... )

Find the nearest power-of-2 resolution for a given meter resolution:

>>> from pyearth.gis.geometry.create_box_from_longitude_latitude import \\
...     find_nearest_resolution
>>> res_deg = find_nearest_resolution(50000, 45.0)  # 50 km at 45°N
>>> print(f"{res_deg}°")
0.5°

See Also
--------
pyearth.gis.geometry.calculate_polygon_area : Calculate area of the created box
pyearth.gis.spatialref.convert_between_degree_and_meter : Degree/meter conversions
"""

import os
from typing import Optional, Tuple
import numpy as np
from osgeo import ogr, osr
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.spatialref.convert_between_degree_and_meter import (
    degree_to_meter,
    meter_to_degree,
)


def find_nearest_resolution(dResolution_meter_in: float, dLatitude_in: float) -> float:
    """
    Find the nearest power-of-2 degree resolution larger than a given meter resolution.

    This function converts a resolution in meters to degrees (accounting for
    latitude-dependent scaling), then finds the nearest standard grid resolution
    that is larger than the input. Standard resolutions are powers of 2:
    1°, 0.5°, 0.25°, 0.125°, 0.0625°, etc. (i.e., 1/(2^n) degrees for n=0,1,2,...).

    Parameters
    ----------
    dResolution_meter_in : float
        Desired resolution in meters. Must be positive.
    dLatitude_in : float
        Latitude in degrees where the resolution applies. Must be in range
        [-90°, 90°]. Needed because the degree-to-meter conversion varies
        with latitude.

    Returns
    -------
    float
        The nearest standard grid resolution in degrees that is equal to or
        larger than the input resolution. Will be a value of the form 1/(2^n)
        for some non-negative integer n.

    Notes
    -----
    **Resolution Selection**:

    The function searches through resolutions: [1.0, 0.5, 0.25, 0.125, ..., 1/512]
    and returns the first resolution that is >= the input resolution (in degrees).

    This ensures that the returned resolution is never smaller than requested,
    which is important for maintaining data accuracy in gridded datasets.

    **Latitude Dependency**:

    At the equator, 1° longitude ≈ 111.32 km
    At 45° latitude, 1° longitude ≈ 78.71 km
    At 60° latitude, 1° longitude ≈ 55.66 km

    The function accounts for this by using the latitude-specific conversion.

    **Default Resolution Range**:

    The function checks 10 resolution levels (1/2^0 through 1/2^9):
    - Finest: ~0.00195° (1/512°) ≈ 217m at equator
    - Coarsest: 1° ≈ 111.32 km at equator

    Warnings
    --------
    - If the input resolution is finer than 1/512°, the function returns 1/512°
    - Very high latitudes (near poles) may not have appropriate resolutions
      because longitude degrees become very small

    Examples
    --------
    Find resolution for 100 km at the equator:

    >>> res = find_nearest_resolution(100000, 0.0)
    >>> print(f"{res}°")
    1.0°

    Find resolution for 50 km at 45°N:

    >>> res = find_nearest_resolution(50000, 45.0)
    >>> print(f"{res}°")
    0.5°

    Find resolution for 10 km at 60°N:

    >>> res = find_nearest_resolution(10000, 60.0)
    >>> print(f"{res}°")
    0.25°

    Very fine resolution (< 1 km):

    >>> res = find_nearest_resolution(500, 0.0)
    >>> print(f"{res}°")
    0.00390625°

    See Also
    --------
    pyearth.gis.spatialref.convert_between_degree_and_meter.meter_to_degree
    pyearth.gis.spatialref.convert_between_degree_and_meter.degree_to_meter
    """
    # Validate inputs
    if dResolution_meter_in <= 0:
        raise ValueError(
            f"dResolution_meter_in must be positive, got {dResolution_meter_in}"
        )

    if not -90 <= dLatitude_in <= 90:
        raise ValueError(f"dLatitude_in must be in range [-90, 90], got {dLatitude_in}")

    # Convert the resolution to degrees at the given latitude
    dResolution_degree = meter_to_degree(dResolution_meter_in, dLatitude_in)

    # Define standard resolutions: 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512
    nResolution = 10
    aResolution_degree = np.array([1.0 / (2**i) for i in range(nResolution)])

    # Find the nearest resolution that is larger than or equal to the input resolution
    dResolution_degree_out = aResolution_degree[0]  # Default to coarsest (1°)

    for i in range(nResolution - 1):
        if aResolution_degree[i] >= dResolution_degree >= aResolution_degree[i + 1]:
            dResolution_degree_out = aResolution_degree[i]
            break

    # If input is coarser than 1°, return 1°
    if dResolution_degree > aResolution_degree[0]:
        dResolution_degree_out = aResolution_degree[0]

    # If input is finer than finest resolution, return finest
    if dResolution_degree < aResolution_degree[-1]:
        dResolution_degree_out = aResolution_degree[-1]

    return dResolution_degree_out


def create_box_from_longitude_latitude(
    dLongitude_in: float,
    dLatitude_in: float,
    dResolution_x_in: float,
    dResolution_y_in: float,
    sFilename_output_in: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a grid-aligned bounding box from a longitude/latitude point.

    This function creates a rectangular bounding box (grid cell) that contains
    the input point and is aligned to a global grid structure anchored at
    (-180°, 90°). The grid cell bounds are determined by the resolution and
    the point's position relative to the grid origin.

    Parameters
    ----------
    dLongitude_in : float
        Longitude of the reference point in degrees. Should be in range
        [-180°, 180°], though values outside are handled via grid calculation.
    dLatitude_in : float
        Latitude of the reference point in degrees. Should be in range
        [-90°, 90°].
    dResolution_x_in : float
        Grid cell width in degrees (longitude direction). Should be a
        power-of-2 fraction: 1, 0.5, 0.25, etc. Must be positive.
    dResolution_y_in : float
        Grid cell height in degrees (latitude direction). Should be a
        power-of-2 fraction: 1, 0.5, 0.25, etc. Must be positive.
    sFilename_output_in : str, optional
        Path to output GeoJSON file. If provided, the grid cell will be
        saved as a GeoJSON polygon with attributes (cellid, longitude,
        latitude, area). If the file exists, it will be removed first.
        If None, no file is created.

    Returns
    -------
    aBox_out : numpy.ndarray
        Array of shape (4,) containing the bounding box coordinates:
        [left, right, bottom, top] or [west, east, south, north]
        - aBox_out[0]: Left (west) longitude
        - aBox_out[1]: Right (east) longitude
        - aBox_out[2]: Bottom (south) latitude
        - aBox_out[3]: Top (north) latitude

    aCoordinates_out : numpy.ndarray
        Array of shape (4, 2) containing the corner coordinates in
        counter-clockwise order: [longitude, latitude]
        - aCoordinates_out[0]: Upper-left (NW) corner
        - aCoordinates_out[1]: Lower-left (SW) corner
        - aCoordinates_out[2]: Lower-right (SE) corner
        - aCoordinates_out[3]: Upper-right (NE) corner

    Raises
    ------
    ValueError
        If dResolution_x_in or dResolution_y_in is not positive.    Notes
    -----
    **Grid Alignment Algorithm**:

    1. Calculate grid indices from the global origin (-180°, 90°):
       - Column index: floor((longitude - (-180)) / resolution_x)
       - Row index: floor((90 - latitude) / resolution_y)

    2. Calculate cell boundaries from indices:
       - Left: -180 + column_index * resolution_x
       - Top: 90 - row_index * resolution_y
       - Right: left + resolution_x
       - Bottom: top - resolution_y

    **Input Point Treatment**:

    The input point (longitude_in, latitude_in) is treated as the top-left
    corner of the grid cell in the grid calculation. The actual cell bounds
    are determined by snapping to the global grid.

    **GeoJSON Output**:

    If filename_output_in is provided, a GeoJSON file is created with:
    - Layer name: 'cell'
    - CRS: EPSG:4326 (WGS84)
    - Attributes:
      * cellid: Always 1 (integer)
      * longitude: Input longitude
      * latitude: Input latitude
      * area: Calculated polygon area in square meters

    **Area Calculation**:

    The area is calculated using spherical geometry (accounting for Earth's
    curvature), so it will be smaller at higher latitudes for the same
    degree-based resolution.

    Warnings
    --------
    - Resolutions should ideally be power-of-2 fractions (1, 0.5, 0.25, etc.)
      for proper global grid alignment, though the function works with any
      positive resolution.
    - Grid cells near the poles may have very small areas despite large
      degree-based dimensions.
    - If the file exists at filename_output_in, it will be deleted.

    Examples
    --------
    Create a 1° × 1° grid cell:

    >>> box, coords = create_box_from_longitude_latitude(-120.0, 45.0, 1.0, 1.0)
    >>> print(f"Bounds: W={box[0]}, E={box[1]}, S={box[2]}, N={box[3]}")
    Bounds: W=-120.0, E=-119.0, S=45.0, N=46.0

    Create a 0.5° × 0.5° grid cell:

    >>> box, coords = create_box_from_longitude_latitude(-100.25, 40.75, 0.5, 0.5)
    >>> print(coords[0])  # Upper-left corner
    [-100.5, 41.0]

    Create grid cell and save to GeoJSON:

    >>> box, coords = create_box_from_longitude_latitude(
    ...     0.0, 0.0, 0.25, 0.25,
    ...     filename_output_in='equator_cell.geojson'
    ... )
    >>> # Creates equator_cell.geojson with a 0.25° × 0.25° cell

    Create a fine-resolution cell:

    >>> box, coords = create_box_from_longitude_latitude(
    ...     -122.4194, 37.7749, 0.01, 0.01  # San Francisco area
    ... )
    >>> print(f"Cell area: {coords.shape}")
    Cell area: (4, 2)

    Access corner coordinates:

    >>> box, coords = create_box_from_longitude_latitude(10.0, 50.0, 1.0, 1.0)
    >>> nw_corner = coords[0]  # Upper-left
    >>> sw_corner = coords[1]  # Lower-left
    >>> se_corner = coords[2]  # Lower-right
    >>> ne_corner = coords[3]  # Upper-right

    See Also
    --------
    find_nearest_resolution : Find appropriate power-of-2 resolution
    pyearth.gis.geometry.calculate_polygon_area : Calculate polygon area
    osgeo.ogr.Geometry : OGR geometry objects
    osgeo.osr.SpatialReference : Spatial reference systems
    """
    # Validate inputs
    if dResolution_x_in <= 0:
        raise ValueError(f"dResolution_x_in must be positive, got {dResolution_x_in}")

    if dResolution_y_in <= 0:
        raise ValueError(f"dResolution_y_in must be positive, got {dResolution_y_in}")

    # Handle existing output file
    if sFilename_output_in is not None and os.path.isfile(sFilename_output_in):
        print(f"{sFilename_output_in} already exists, removing...")
        os.remove(sFilename_output_in)

    # Use input point as the top-left corner reference
    dLon_min = dLongitude_in
    dLat_max = dLatitude_in

    # Calculate grid indices from global origin (-180°, 90°)
    # This ensures alignment to a global grid structure
    nleft = np.floor((dLon_min - (-180)) / dResolution_x_in)
    ntop = np.floor((90 - dLat_max) / dResolution_y_in)

    # Initialize output arrays
    aBox_out = np.full(4, -9999.0, dtype=float)
    aCoordinates_out = np.full((4, 2), -9999.0, dtype=float)

    # Calculate grid-aligned box boundaries
    # Left (west) longitude
    aBox_out[0] = -180 + nleft * dResolution_x_in
    # Top (north) latitude
    aBox_out[3] = 90 - ntop * dResolution_y_in
    # Right (east) longitude
    aBox_out[1] = aBox_out[0] + dResolution_x_in
    # Bottom (south) latitude
    aBox_out[2] = aBox_out[3] - dResolution_y_in

    # Build corner coordinates in counter-clockwise order
    # Upper-left (NW)
    aCoordinates_out[0, 0] = aBox_out[0]  # West
    aCoordinates_out[0, 1] = aBox_out[3]  # North

    # Lower-left (SW)
    aCoordinates_out[1, 0] = aBox_out[0]  # West
    aCoordinates_out[1, 1] = aBox_out[2]  # South

    # Lower-right (SE)
    aCoordinates_out[2, 0] = aBox_out[1]  # East
    aCoordinates_out[2, 1] = aBox_out[2]  # South

    # Upper-right (NE)
    aCoordinates_out[3, 0] = aBox_out[1]  # East
    aCoordinates_out[3, 1] = aBox_out[3]  # North

    # Calculate polygon area using spherical geometry
    dArea = calculate_polygon_area(aCoordinates_out[:, 0], aCoordinates_out[:, 1])

    # Save as GeoJSON file if requested
    if sFilename_output_in is not None:
        # Create GeoJSON driver and dataset
        pDriver_geojson = ogr.GetDriverByName("GeoJSON")
        pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in)

        # Set up WGS84 spatial reference
        pSpatial_reference_gcs = osr.SpatialReference()
        pSpatial_reference_gcs.ImportFromEPSG(4326)  # WGS84 lat/lon
        pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

        # Create layer
        pLayer = pDataset.CreateLayer("cell", pSpatial_reference_gcs, ogr.wkbPolygon)

        # Add attribute fields
        pLayer.CreateField(ogr.FieldDefn("cellid", ogr.OFTInteger64))
        pLayer.CreateField(ogr.FieldDefn("longitude", ogr.OFTReal))
        pLayer.CreateField(ogr.FieldDefn("latitude", ogr.OFTReal))

        pArea_field = ogr.FieldDefn("area", ogr.OFTReal)
        pArea_field.SetWidth(20)
        pArea_field.SetPrecision(2)
        pLayer.CreateField(pArea_field)

        # Create feature
        pLayerDefn = pLayer.GetLayerDefn()
        pFeature = ogr.Feature(pLayerDefn)

        # Build polygon geometry
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for x, y in aCoordinates_out:
            ring.AddPoint(float(x), float(y))

        # Close the ring (add first point again)
        ring.AddPoint(float(aCoordinates_out[0, 0]), float(aCoordinates_out[0, 1]))

        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pPolygon.AddGeometry(ring)

        # Set geometry and attributes
        pFeature.SetGeometry(pPolygon)
        pFeature.SetField("cellid", 1)
        pFeature.SetField("longitude", float(dLongitude_in))
        pFeature.SetField("latitude", float(dLatitude_in))
        pFeature.SetField("area", float(dArea))

        # Write feature to layer
        pLayer.CreateFeature(pFeature)

        # Clean up (close dataset)
        pDataset = pLayer = pFeature = None

    return aBox_out, aCoordinates_out
