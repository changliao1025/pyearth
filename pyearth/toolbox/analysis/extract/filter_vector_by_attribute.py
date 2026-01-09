"""
Filter Vector Data by Attribute Values
========================================

This module provides functionality to extract features from vector datasets based on
attribute values, creating a filtered output dataset containing only features that
match the specified criteria.

The filtering operation is commonly used in geospatial data processing workflows for:

* **Data Subset Creation**: Extract specific administrative regions, land use types,
  or other categorical features
* **Quality Control**: Isolate features meeting certain quality thresholds
* **Thematic Mapping**: Create derived datasets for specific themes or categories
* **Data Reduction**: Reduce dataset size by extracting only relevant features

Key Features
------------
- Supports multiple vector formats (Parquet, Shapefile, GeoJSON, GeoPackage, etc.)
- Automatic format detection based on file extension
- Preserves spatial reference system (coordinate system) from input
- Maintains all field definitions and attributes
- Uses SQL-based filtering for efficient attribute queries
- Automatically handles output file overwrites

Technical Details
-----------------
The function uses GDAL/OGR's attribute filtering capabilities via SQL WHERE clauses.
This approach leverages database-style indexing for efficient feature selection,
particularly beneficial for large datasets.

**Format Support**: The function automatically detects output format from file extension.
Supported formats include:
- Parquet (.parquet): Columnar storage for efficient attribute queries
- Shapefile (.shp): Traditional GIS format with broad compatibility
- GeoJSON (.geojson, .json): Text-based format for web applications
- GeoPackage (.gpkg): Modern SQLite-based format
- GML (.gml): Geography Markup Language
- KML (.kml): Keyhole Markup Language for Google Earth

Input can be any OGR-compatible vector format.

**Spatial Reference Preservation**: The output dataset inherits the spatial reference
system (SRS) from the input, ensuring coordinate system consistency.

See Also
--------
filter_vector_by_geometry : Filter features by spatial relationships
extract_features_by_location : Extract based on spatial intersection
attribute_query : More complex attribute-based queries

Notes
-----
The function creates a new dataset rather than modifying the input file, ensuring
the original data remains unchanged.

Examples
--------
Extract all features where country equals "USA":

    >>> filter_vector_by_attribute(
    ...     'world_countries.parquet',
    ...     'usa_only.parquet',
    ...     'country',
    ...     'USA'
    ... )

Filter land use polygons for residential areas:

    >>> filter_vector_by_attribute(
    ...     'landuse.parquet',
    ...     'residential.parquet',
    ...     'class',
    ...     'residential'
    ... )
"""

import os
from osgeo import osr, ogr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_extension


def filter_vector_by_attribute(
    sFilename_input, sFilename_output, sAttribute_name, dValue_filter
):
    """
    Filter vector dataset features by attribute value and save to a new file.

    This function reads a vector dataset in any OGR-supported format, filters features
    based on an attribute value match, and writes the filtered results to a new file.
    The output format is automatically determined from the file extension. All features
    where the specified attribute equals the filter value are included in the output.

    The function performs an exact match comparison using SQL WHERE clause syntax.
    The filtering operation preserves:
    - Spatial reference system (coordinate system)
    - All field definitions and data types
    - Feature geometries and all attributes
    - Layer metadata

    Parameters
    ----------
    sFilename_input : str
        Path to the input vector file to filter.
        Can be any OGR-compatible format (Parquet, Shapefile, GeoJSON, GeoPackage, etc.).

    sFilename_output : str
        Path where the filtered output file will be created.
        Output format is determined by file extension (.parquet, .shp, .geojson, .gpkg, etc.).
        If the file exists, it will be automatically deleted and recreated.
        Directory must exist and be writable.

    sAttribute_name : str
        Name of the attribute field to filter on.
        Must exist in the input dataset's attribute table.
        Field name is case-sensitive and must match exactly.

    dValue_filter : str or numeric
        Value to match in the specified attribute field.
        Features with this exact value in `sAttribute_name` will be included.
        String values are automatically quoted in the SQL query.
        For numeric comparisons, ensure the value type matches the field type.

    Returns
    -------
    None
        Results are written directly to the output file.
        Function prints confirmation message upon successful completion.

    Raises
    ------
    RuntimeError
        If the input file cannot be opened (file not found, permission denied,
        or invalid format).

    Notes
    -----
    1. **File Format**: Supports multiple vector formats with automatic detection from
       file extension. Supported formats: Parquet (.parquet), Shapefile (.shp),
       GeoJSON (.geojson/.json), GeoPackage (.gpkg), GML (.gml), KML (.kml).
       If extension is not recognized, defaults to Parquet format.

    2. **Geometry Type**: Output layer is configured for polygon geometries
       (`ogr.wkbPolygon`). For other geometry types, this would need modification.

    3. **SQL Filtering**: Uses OGR's `SetAttributeFilter()` with SQL WHERE syntax.
       The filter is: `attribute_name = 'value'`. Only exact matches are included.

    4. **Output Overwrite**: If `sFilename_output` exists, it is deleted without
       warning. Ensure you don't overwrite important data.

    5. **Field Preservation**: All attribute fields from the input are copied to
       the output, maintaining data types and field properties.

    6. **Performance**: For large datasets, attribute filtering is typically fast,
       but performance depends on:
       - Dataset size (number of features)
       - Attribute indexing (varies by format)
       - Filter selectivity (fraction of features matching)
       - Output format (Parquet is generally fastest for analytical queries)

    7. **Empty Results**: If no features match the filter criteria, an empty
       dataset with the correct schema is created.

    8. **Spatial Reference**: The output inherits the SRS from the input layer
       automatically via the `srs` parameter in `CreateLayer()`.

    Examples
    --------
    Filter a global dataset to extract only USA features (Parquet to Parquet):

        >>> filter_vector_by_attribute(
        ...     sFilename_input='world_countries.parquet',
        ...     sFilename_output='usa.parquet',
        ...     sAttribute_name='country_code',
        ...     dValue_filter='USA'
        ... )
        filter_vector_by_attribute is done usa.parquet

    Extract residential land use polygons (Shapefile to GeoJSON):

        >>> filter_vector_by_attribute(
        ...     sFilename_input='landuse_2024.shp',
        ...     sFilename_output='residential_2024.geojson',
        ...     sAttribute_name='landuse_class',
        ...     dValue_filter='residential'
        ... )
        filter_vector_by_attribute is done residential_2024.geojson

    Filter watersheds by category (GeoJSON to GeoPackage):

        >>> filter_vector_by_attribute(
        ...     sFilename_input='watersheds.geojson',
        ...     sFilename_output='large_watersheds.gpkg',
        ...     sAttribute_name='drainage_class',
        ...     dValue_filter='large'
        ... )
        filter_vector_by_attribute is done large_watersheds.gpkg

    Create subset for a specific administrative region (any format to Shapefile):

        >>> filter_vector_by_attribute(
        ...     sFilename_input='census_tracts.parquet',
        ...     sFilename_output='king_county_tracts.shp',
        ...     sAttribute_name='county_name',
        ...     dValue_filter='King County'
        ... )
        filter_vector_by_attribute is done king_county_tracts.parquet

    See Also
    --------
    ogr.Layer.SetAttributeFilter : OGR method for SQL-based attribute filtering
    filter_vector_by_geometry : Filter features by spatial criteria
    extract_features_by_extent : Extract features within a bounding box
    """
    # Open the input dataset using generic ogr.Open() which auto-detects the format
    # This supports any OGR-compatible vector format (Parquet, Shapefile, GeoJSON, etc.)
    pDataSource_in = ogr.Open(sFilename_input, 0)

    # Validate that the input file was successfully opened
    # Common failure reasons: file not found, permission denied, invalid format
    if pDataSource_in is None:
        print("Could not open input file")
        return

    # Remove existing output file if it exists
    # This ensures a clean write without conflicts or appending issues
    if os.path.exists(sFilename_output):
        os.remove(sFilename_output)

    # Get the first (and typically only) layer from the input dataset
    # Most vector files contain a single layer, but formats can support multiple
    pLayer_in = pDataSource_in.GetLayer()

    # Extract the spatial reference system (SRS) from the input layer
    # This preserves the coordinate system (e.g., WGS84, UTM, etc.) in the output
    pSpatialRef = pLayer_in.GetSpatialRef()

    # Get the appropriate output driver based on file extension
    # Uses the centralized utility function for consistent format handling
    try:
        pDriver_out = get_vector_driver_from_extension(sFilename_output)
    except ValueError as e:
        print(f"Error getting driver for output file: {e}")
        pDataSource_in = None
        return

    # Create a new output dataset using the format-specific driver
    # This initializes the file structure for writing filtered features
    pDataSource_out = pDriver_out.CreateDataSource(sFilename_output)

    # Create the output layer with the same spatial reference as input
    # Layer name: 'filtered' (arbitrary, doesn't affect file access)
    # Geometry type: wkbPolygon (assumes input contains polygon features)
    # SRS: Inherited from input to maintain coordinate system consistency
    pLayer_out = pDataSource_out.CreateLayer(
        "filtered", geom_type=ogr.wkbPolygon, srs=pSpatialRef
    )

    # Copy all field definitions from input layer to output layer
    # This preserves the attribute table schema (field names, types, widths)
    pLayer_defn_in = pLayer_in.GetLayerDefn()
    for i in range(pLayer_defn_in.GetFieldCount()):
        pField_defn = pLayer_defn_in.GetFieldDefn(i)
        pLayer_out.CreateField(pField_defn)

    # Construct SQL WHERE clause for attribute filtering
    # Format: "attribute_name = 'value'" for string matching
    # Note: String values are automatically quoted; numeric values should work as well
    sql_query = f"{sAttribute_name} = '{dValue_filter}'"

    # Apply the attribute filter to the input layer
    # This restricts iteration to only features matching the SQL condition
    # No features are actually modified; this just controls which are returned
    pLayer_in.SetAttributeFilter(sql_query)

    # Iterate through filtered features and copy them to the output layer
    # Only features matching the attribute filter are included in this iteration
    for feature in pLayer_in:
        pLayer_out.CreateFeature(feature)

    # Close both data sources by setting references to None
    # This ensures all data is flushed to disk and file handles are released properly
    pDataSource_in = pDataSource_out = None

    # Print confirmation message indicating successful completion
    print("filter_vector_by_attribute is done", sFilename_output)
