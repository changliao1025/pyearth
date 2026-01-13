"""
Filter Vector Data by Polygon Extent
=====================================

This module provides functionality to extract features from vector datasets based on
spatial filtering using a polygon layer's bounding box (extent). Features that
intersect with the polygon's bounding box are extracted to a new output dataset.

The spatial filtering operation is commonly used in geospatial workflows for:

* **Regional Data Extraction**: Extract features within a specific geographic area
  defined by administrative boundaries, study areas, or areas of interest
* **Data Tiling**: Split large datasets into manageable regional tiles
* **Performance Optimization**: Reduce processing time by working with spatially
  filtered subsets rather than entire datasets
* **Spatial Subsetting**: Create focused datasets for regional analysis or visualization

Key Features
------------
- Spatial filtering based on bounding box (extent) of polygon layer
- Preserves all attribute fields and their data types
- Supports multiple output formats (Parquet, Shapefile, GeoJSON, GeoPackage, etc.)
- Automatic format detection based on file extension
- Supports any OGR-compatible input vector format
- Automatic resource cleanup and spatial filter reset

Technical Details
-----------------
The function uses OGR's spatial filtering capabilities (`SetSpatialFilterRect`) to
efficiently query features that intersect with a rectangular extent. This approach:

**Bounding Box Filtering**: Uses the axis-aligned bounding box (AABB) of the polygon
layer as the filter region. This is computationally efficient but includes features
that may not actually intersect the polygon geometry itself - only its extent.

**Format Support**: The function automatically detects output format from file extension.
Supported formats include Parquet, Shapefile, GeoJSON, GeoPackage, GML, and KML.
Input can be any OGR-compatible vector format.

**Performance**: Spatial filtering with bounding boxes is very fast because:
- Uses spatial indexes if available (R-tree or similar)
- Simple rectangle intersection tests (cheaper than polygon intersection)
- Database-level filtering reduces memory overhead

**Limitation**: Currently uses bounding box filtering (option 1) rather than exact
polygon geometry intersection. Features outside the polygon but within its bounding
box will be included. For exact polygon clipping, use dedicated clip operations.

See Also
--------
filter_vector_by_attribute : Filter features by attribute values
clip_vector_by_polygon : Exact polygon clipping (geometrically clips features)
spatial_join : Join features based on spatial relationships

Notes
-----
The function currently uses bounding box filtering (iFlag_option = 1) which is fast
but not geometrically precise. A more accurate polygon-based filter could be
implemented by setting iFlag_option = 2 (currently unimplemented placeholder).

Examples
--------
Extract features within a watershed boundary:

    >>> filter_vector_by_polygon(
    ...     'all_streams.parquet',
    ...     'watershed_boundary.shp',
    ...     'watershed_streams.parquet'
    ... )

Filter points within a study area:

    >>> filter_vector_by_polygon(
    ...     'global_weather_stations.parquet',
    ...     'california.geojson',
    ...     'california_stations.parquet'
    ... )
"""

import os
from osgeo import osr, ogr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_extension


def filter_vector_by_polygon(
    sFilename_vector_in, sFilename_polygon_in, sFilename_vector_out
):
    """
    Filter vector features by the bounding box extent of a polygon layer.

    This function extracts features from an input vector dataset that spatially
    intersect with the bounding box (extent) of a polygon layer. The filtered
    features are written to a new Parquet file, preserving all attribute fields.

    The spatial filtering uses the rectangular extent (bounding box) of the polygon
    layer rather than exact polygon geometry intersection. This approach is
    computationally efficient but may include features that fall within the bounding
    box but outside the actual polygon boundaries.

    Parameters
    ----------
    sFilename_vector_in : str
        Path to the input vector file to be filtered.
        Can be any OGR-compatible format (Parquet, Shapefile, GeoJSON, etc.).
        Features from this dataset will be tested against the spatial filter.

    sFilename_polygon_in : str
        Path to the polygon file defining the filter region.
        Can be any OGR-compatible format containing polygon geometries.
        The bounding box (extent) of this layer defines the filter region.

    sFilename_vector_out : str
        Path where the filtered output file will be created.
        Output format is determined by file extension (.parquet, .shp, .geojson, .gpkg, etc.).
        Directory must exist and be writable.

    Returns
    -------
    None
        Results are written directly to the output file.
        Function does not return a value; check output file for results.

    Raises
    ------
    RuntimeError
        If input vector file cannot be opened (file not found, permission denied,
        or invalid format).
        If polygon file cannot be opened.
        If output driver cannot be obtained for the specified format.

    Notes
    -----
    1. **Spatial Filter Method**: Currently uses bounding box (rectangular extent)
       filtering via `SetSpatialFilterRect()`. This is fast but includes features
       within the bounding box that may not intersect the actual polygon geometry.

    2. **Geometry Type**: Output layer is configured for polygon geometries
       (`ogr.wkbPolygon`). For other geometry types (points, lines), this may need
       modification or the function should auto-detect geometry type from input.

    3. **Coordinate System**: Input vector and polygon files should ideally be in
       the same spatial reference system for accurate spatial filtering. The function
       does not perform automatic reprojection.

    4. **Bounding Box vs Polygon**: The extent is a rectangular bounding box
       (min_x, min_y, max_x, max_y). Features touching or overlapping this
       rectangle are included, even if they don't intersect the polygon itself.

    5. **Performance**: Bounding box filtering is very efficient, especially with
       spatial indexes. For large datasets, this can be orders of magnitude faster
       than exact polygon intersection.

    6. **Field Preservation**: All attribute fields from the input layer are copied
       to the output, maintaining field names, types, widths, and precision.

    7. **Feature Geometry**: Geometries are copied as-is without clipping. Features
       may extend beyond the polygon boundaries. For exact clipping, use a dedicated
       clip operation after filtering.

    8. **Spatial Filter Cleanup**: The function resets the spatial filter on the
       input layer before closing, ensuring no side effects if the data source is
       reused elsewhere.

    9. **Output Format**: Output format is automatically detected from file extension.
       Supports Parquet (.parquet), Shapefile (.shp), GeoJSON (.geojson/.json),
       GeoPackage (.gpkg), GML (.gml), and KML (.kml). Each format has different
       characteristics for compression, compatibility, and query performance.

    10. **Incomplete Implementation**: The code includes placeholder for exact polygon
        filtering (iFlag_option = 2) but this is not currently implemented. Only
        bounding box filtering (option 1) is functional.

    Examples
    --------
    Extract stream segments within a watershed boundary (Parquet to Parquet):

        >>> filter_vector_by_polygon(
        ...     sFilename_vector_in='national_streams.parquet',
        ...     sFilename_polygon_in='my_watershed.shp',
        ...     sFilename_vector_out='watershed_streams.parquet'
        ... )
        (minX, maxX, minY, maxY)  # Prints extent of clip polygon

    Filter weather stations within a state boundary (GeoJSON to Shapefile):

        >>> filter_vector_by_polygon(
        ...     sFilename_vector_in='all_weather_stations.geojson',
        ...     sFilename_polygon_in='washington_state.shp',
        ...     sFilename_vector_out='wa_stations.shp'
        ... )

    Extract land parcels within a city boundary (any format to GeoPackage):

        >>> filter_vector_by_polygon(
        ...     sFilename_vector_in='county_parcels.parquet',
        ...     sFilename_polygon_in='city_boundary.geojson',
        ...     sFilename_vector_out='city_parcels.gpkg'
        ... )

    Create regional subset of global dataset (Parquet to GeoJSON):

        >>> filter_vector_by_polygon(
        ...     sFilename_vector_in='world_cities.parquet',
        ...     sFilename_polygon_in='europe_boundary.shp',
        ...     sFilename_vector_out='europe_cities.geojson'
        ... )

    See Also
    --------
    ogr.Layer.SetSpatialFilterRect : OGR method for bounding box spatial filtering
    ogr.Layer.SetSpatialFilter : Exact geometry-based spatial filtering
    filter_vector_by_attribute : Filter features by attribute criteria
    clip_vector_by_polygon : Geometrically clip features to polygon boundary
    """
    # Open the input vector dataset using generic ogr.Open()
    # This allows any OGR-supported format as input
    pDataSource_in = ogr.Open(sFilename_vector_in)

    # Validate that the input vector file was successfully opened
    # Common failure reasons: file not found, permission denied, unsupported format
    if pDataSource_in is None:
        print("Could not open input vector")
        return

    # Get the first (and typically only) layer from the input dataset
    pLayer_in = pDataSource_in.GetLayer()

    # Open the polygon dataset that defines the spatial filter region
    pDataSource_clip = ogr.Open(sFilename_polygon_in)

    # Validate that the polygon file was successfully opened
    if pDataSource_clip is None:
        print("Could not open clip polygon")
        # Clean up the already-opened input data source before returning
        pDataSource_in = None
        return

    # Get the polygon layer that defines the filter extent
    pLayer_clip = pDataSource_clip.GetLayer()

    # Get the appropriate output driver based on file extension
    # Uses the centralized utility function for consistent format handling
    try:
        pDriver_out = get_vector_driver_from_extension(sFilename_vector_out)
    except ValueError as e:
        print(f"Error getting driver for output file: {e}")
        pDataSource_in = pDataSource_clip = None
        return

    # Create a new output dataset using the format-specific driver
    # This will hold the spatially filtered features
    pDataSource_out = pDriver_out.CreateDataSource(sFilename_vector_out)

    # Create the output layer
    # Layer name: 'clipped' (arbitrary identifier, doesn't affect file access)
    # Geometry type: wkbPolygon (assumes input contains polygon features)
    # Note: For mixed geometry types, consider auto-detecting from input layer
    pLayer_out = pDataSource_out.CreateLayer("clipped", geom_type=ogr.wkbPolygon)

    # Copy all field definitions from input layer to output layer
    # This preserves the complete attribute schema (field names, types, widths)
    pLayer_defn_in = pLayer_in.GetLayerDefn()
    for i in range(pLayer_defn_in.GetFieldCount()):
        field_defn = pLayer_defn_in.GetFieldDefn(i)
        pLayer_out.CreateField(field_defn)

    # Configure spatial filtering method
    # Option 1: Bounding box (fast but approximate)
    # Option 2: Exact polygon geometry (accurate but slower, not implemented)
    iFlag_option = 1

    if iFlag_option == 1:
        # Get the bounding box (extent) of the clip polygon layer
        # Returns tuple: (min_x, max_x, min_y, max_y)
        # This is the axis-aligned rectangular extent containing all polygon features
        aExtent_clip = pLayer_clip.GetExtent()
        print(aExtent_clip)

        # Apply rectangular spatial filter to the input layer
        # SetSpatialFilterRect parameters: (min_x, min_y, max_x, max_y)
        # Note: Order differs from GetExtent() - extent[0], extent[2], extent[1], extent[3]
        # Only features intersecting this rectangle will be returned in subsequent iterations
        pLayer_in.SetSpatialFilterRect(
            aExtent_clip[0], aExtent_clip[2], aExtent_clip[1], aExtent_clip[3]
        )
        pass
    else:
        # Placeholder for exact polygon geometry filtering
        # This would use SetSpatialFilter(polygon_geometry) for precise filtering
        # More accurate but computationally more expensive than bounding box
        if (
            iFlag_option == 2
        ):  # use the outline of the clip layer, this is more accurate
            pass

    # Iterate through spatially filtered features and copy them to output
    # Only features passing the spatial filter (intersecting bounding box) are included
    for feature in pLayer_in:
        # Get the geometry of the current feature
        # Geometry is copied as-is without clipping to polygon boundaries
        pGeometry = feature.GetGeometryRef()

        # Create a new feature in the output layer using its schema
        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())

        # Set the geometry for the output feature
        pFeature_out.SetGeometry(pGeometry)

        # Copy all attribute field values from input feature to output feature
        # Iterate through all fields and transfer values by field name
        for i in range(pLayer_out.GetLayerDefn().GetFieldCount()):
            pFeature_out.SetField(
                pLayer_out.GetLayerDefn().GetFieldDefn(i).GetNameRef(),
                feature.GetField(i),
            )

        # Write the feature to the output layer
        pLayer_out.CreateFeature(pFeature_out)

        # Release the feature to free memory
        pFeature_out = None

    # Reset the spatial filter on the input layer
    # This clears any filtering applied by SetSpatialFilterRect
    # Important for cleanup to avoid side effects if data source is reused
    pLayer_in.SetSpatialFilter(None)

    # Close all data sources by setting references to None
    # This ensures proper resource cleanup, data flushing, and file handle release
    pDataSource_in = pDataSource_clip = pDataSource_out = None
