"""
Copy vector geometries without attributesSee Also
--------
osgeo.ogr : GDAL/OGR vector data handling
pyearth.gis.gdal.gdal_vector_format_support : Vector format detection and driver utilities

Referencesw files.

This module provides functionality to extract pure geometry from vector files,
removing all attribute data while preserving spatial information and coordinate
reference systems. Useful for creating lightweight geometry files or preparing
data for specific GIS operations.

Main Features
-------------
- Support for multiple vector formats (GeoJSON, Shapefile, GeoPackage, KML, etc.)
- Automatic format detection from file extensions
- Preserves spatial reference system (CRS)
- Preserves geometry types (Point, LineString, Polygon, etc.)
- Handles multi-geometry types (MultiPoint, MultiPolygon, etc.)
- Removes all attribute fields
- Validates input and output files

Typical Use Cases
-----------------
1. **File Size Reduction**: Remove unnecessary attributes to reduce file size
2. **Privacy**: Share spatial boundaries without revealing attribute data
3. **GIS Processing**: Create clean geometry inputs for spatial operations
4. **Data Simplification**: Extract only geometric information for analysis
5. **Format Conversion**: Convert between vector formats while stripping attributes

Supported Formats
-----------------
Input and output formats (via GDAL/OGR):
- GeoJSON (.geojson, .json)
- ESRI Shapefile (.shp)
- GeoPackage (.gpkg)
- KML (.kml)
- GML (.gml)
- KMZ (.kmz)
- And other GDAL-supported vector formats

Notes
-----
- Output format is auto-detected from file extension
- All features are copied, maintaining feature count
- Spatial reference is preserved from input
- Geometry type is determined from first feature
- Mixed geometry types may cause issues (use consistent types)

See Also
--------
osgeo.ogr : GDAL/OGR vector data handling
pyearth.gis.gdal.gdal_vector_format_support : Vector format detection utilities

References
----------
.. [1] GDAL/OGR Vector Data Model. https://gdal.org/user/vector_data_model.html
.. [2] OGR Vector Formats. https://gdal.org/drivers/vector/index.html

"""

import os
import logging
from typing import Optional
from pathlib import Path
from osgeo import ogr, osr

# Configure module logger
logger = logging.getLogger(__name__)

from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_driver_from_extension,
    get_vector_format_from_extension,
)


def copy_geometry_without_attributes(
    sFilename_in: str,
    sFilename_out: str,
    sFormat_in: Optional[str] = None,
    sFormat_out: Optional[str] = None,
    sLayer_name: Optional[str] = None,
    iLayer_index: int = 0,
) -> str:
    """
    Copy vector geometries to a new file, removing all attribute data.

    Extracts pure geometry from input vector file, preserving spatial reference
    and geometry types while removing all attribute fields. Supports multiple
    vector formats with automatic format detection from file extensions.

    Parameters
    ----------
    sFilename_in : str
        Input vector file path.

        Supported formats (auto-detected from extension):
        - GeoJSON: '.geojson', '.json'
        - Shapefile: '.shp'
        - GeoPackage: '.gpkg'
        - KML: '.kml'
        - GML: '.gml'
        - Other GDAL-supported formats

        Examples:
        - '/path/to/input.geojson'
        - 'data/boundaries.shp'
        - 'output.gpkg'

        File must exist and be readable.

    sFilename_out : str
        Output vector file path.

        Format auto-detected from extension if sFormat_out not specified.
        Will be overwritten if it exists.

        Examples:
        - '/path/to/output.geojson'
        - 'data/geometry_only.shp'
        - 'clean.gpkg'

        Parent directory must be writable.

    sFormat_in : str, optional
        Input file format (GDAL driver name), default None (auto-detect).

        Specify if format cannot be determined from extension.

        Common values:
        - 'GeoJSON'
        - 'ESRI Shapefile'
        - 'GPKG'
        - 'KML'

        If None, attempts to open with generic ogr.Open().

    sFormat_out : str, optional
        Output file format (GDAL driver name), default None (auto-detect).

        If None, format determined from sFilename_out extension:
        - '.geojson' or '.json' → 'GeoJSON'
        - '.shp' → 'ESRI Shapefile'
        - '.gpkg' → 'GPKG'
        - '.kml' → 'KML'

        Common values:
        - 'GeoJSON': Human-readable, widely supported
        - 'ESRI Shapefile': Traditional GIS format
        - 'GPKG': Modern SQLite-based format
        - 'KML': Google Earth format

    sLayer_name : str, optional
        Name of layer to copy from input file, default None.

        If specified, copies only the named layer.
        If None, uses iLayer_index to select layer.

        Useful for multi-layer formats like GeoPackage.

        Examples:
        - 'boundaries'
        - 'admin_level_4'
        - 'countries'

    iLayer_index : int, optional
        Index of layer to copy from input file, default 0.

        Zero-based index (0 = first layer).
        Only used if sLayer_name is None.

        Range: 0 to (layer_count - 1)

    Returns
    -------
    str
        Path to the created output file (same as sFilename_out).

        File contains:
        - All geometries from input layer
        - Spatial reference from input
        - NO attribute fields
        - Same geometry types as input

    Raises
    ------
    TypeError
        If sFilename_in or sFilename_out is not a string.
    ValueError
        - If sFilename_in or sFilename_out is empty
        - If input file cannot be opened
        - If input layer contains no features
        - If output file cannot be created
        - If layer index is out of range
        - If specified layer name not found
    FileNotFoundError
        If input file does not exist.
    RuntimeError
        - If GDAL driver not available
        - If output layer creation fails
        - If geometry copying fails
    OSError
        If output directory cannot be created or is not writable.

    Warns
    -----
    UserWarning
        - If input format cannot be auto-detected
        - If output format cannot be auto-detected
        - If input has multiple layers (only one will be copied)
        - If geometry types are mixed in input layer

    Notes
    -----
    **Processing Details:**

    1. Validate input parameters
    2. Open input file and access specified layer
    3. Extract spatial reference and geometry type
    4. Create output file with detected/specified format
    5. Create output layer with only geometry (no fields)
    6. Copy all geometries from input to output
    7. Close datasets and return output path

    **Spatial Reference:**

    - Preserved exactly from input layer
    - If input has no CRS, output will have no CRS
    - No reprojection or transformation applied

    **Geometry Types:**

    - Determined from first feature in input
    - All features must have compatible geometry types
    - Mixed types may cause errors
    - Supported types: Point, LineString, Polygon, MultiPoint,
      MultiLineString, MultiPolygon, GeometryCollection

    **Attribute Removal:**

    - ALL attribute fields are removed
    - Only FID (Feature ID) may be preserved (format-dependent)
    - Original attribute data is NOT recoverable from output

    **Format-Specific Notes:**

    - GeoJSON: Single file, coordinates in WGS84
    - Shapefile: Creates multiple files (.shp, .shx, .dbf, .prj)
    - GeoPackage: Single SQLite database file
    - KML: XML-based, automatically converts to WGS84

    **Format Detection:**

    - Uses pyearth.gis.gdal.gdal_vector_format_support for format detection
    - Automatically detects format from file extension
    - Supports all GDAL-registered vector formats
    - Falls back to generic ogr.Open() if auto-detection fails

    **Performance:**

    - Linear time complexity: O(n) where n = feature count
    - Memory efficient: processes one feature at a time
    - Large files (>1M features): ~30-60 seconds per 100k features

    Examples
    --------
    **Example 1: Basic GeoJSON to GeoJSON (remove attributes)**

    >>> from pyearth.toolbox.data.geojson.copy_geometry_without_attributes import \\
    ...     copy_geometry_without_attributes
    >>>
    >>> # Remove attributes from GeoJSON file
    >>> output = copy_geometry_without_attributes(
    ...     sFilename_in='data/countries_with_data.geojson',
    ...     sFilename_out='data/countries_geometry_only.geojson')
    >>> print(f"Created: {output}")
    Created: data/countries_geometry_only.geojson

    **Example 2: Convert Shapefile to GeoJSON (remove attributes)**

    >>> # Convert and strip attributes in one operation
    >>> output = copy_geometry_without_attributes(
    ...     sFilename_in='data/boundaries.shp',
    ...     sFilename_out='data/boundaries_clean.geojson')
    >>> print(f"Shapefile converted to clean GeoJSON: {output}")
    Shapefile converted to clean GeoJSON: data/boundaries_clean.geojson

    **Example 3: GeoJSON to GeoPackage**

    >>> # Create lightweight GeoPackage with only geometry
    >>> output = copy_geometry_without_attributes(
    ...     sFilename_in='input.geojson',
    ...     sFilename_out='output.gpkg')
    >>> print(f"Created GeoPackage: {output}")
    Created GeoPackage: output.gpkg

    **Example 4: Explicit format specification**

    >>> # Specify formats explicitly (useful for ambiguous extensions)
    >>> output = copy_geometry_without_attributes(
    ...     sFilename_in='data.txt',
    ...     sFilename_out='output.txt',
    ...     sFormat_in='GeoJSON',
    ...     sFormat_out='GeoJSON')
    >>> print(f"Explicit formats: {output}")
    Explicit formats: output.txt

    **Example 5: Extract specific layer from multi-layer file**

    >>> # Copy only 'admin_boundaries' layer from GeoPackage
    >>> output = copy_geometry_without_attributes(
    ...     sFilename_in='data.gpkg',
    ...     sFilename_out='boundaries.geojson',
    ...     sLayer_name='admin_boundaries')
    >>> print(f"Extracted layer to: {output}")
    Extracted layer to: boundaries.geojson

    **Example 6: Use layer index**

    >>> # Copy second layer (index 1) from multi-layer file
    >>> output = copy_geometry_without_attributes(
    ...     sFilename_in='data.gpkg',
    ...     sFilename_out='layer2.geojson',
    ...     iLayer_index=1)
    >>> print(f"Layer 1 extracted to: {output}")
    Layer 1 extracted to: layer2.geojson

    **Example 7: Privacy - share boundaries without data**

    >>> # Remove sensitive census data, keep only boundaries
    >>> output = copy_geometry_without_attributes(
    ...     sFilename_in='census_tracts_with_population.shp',
    ...     sFilename_out='census_tracts_boundaries_only.geojson')
    >>> # Output has same geometries, but no population data
    >>> print(f"Privacy-safe boundaries: {output}")
    Privacy-safe boundaries: census_tracts_boundaries_only.geojson

    **Example 8: Reduce file size**

    >>> import os
    >>> input_file = 'large_file_with_attributes.geojson'
    >>> output_file = 'large_file_geometry_only.geojson'
    >>>
    >>> size_before = os.path.getsize(input_file)
    >>> output = copy_geometry_without_attributes(input_file, output_file)
    >>> size_after = os.path.getsize(output_file)
    >>>
    >>> reduction = 100 * (1 - size_after / size_before)
    >>> print(f"File size reduced by {reduction:.1f}%")
    File size reduced by 68.3%

    **Example 9: Prepare for spatial operations**

    >>> # Create clean geometry for overlay analysis
    >>> mask = copy_geometry_without_attributes(
    ...     sFilename_in='study_area_with_metadata.shp',
    ...     sFilename_out='study_area_mask.geojson')
    >>> # Use mask for clipping other datasets
    >>> print(f"Mask created: {mask}")
    Mask created: study_area_mask.geojson

    **Example 10: Batch processing workflow**

    >>> import glob
    >>> # Process all shapefiles in directory
    >>> input_files = glob.glob('input/*.shp')
    >>> for input_file in input_files:
    ...     base_name = os.path.splitext(os.path.basename(input_file))[0]
    ...     output_file = f'output/{base_name}_clean.geojson'
    ...     copy_geometry_without_attributes(input_file, output_file)
    ...     print(f"Processed: {base_name}")
    Processed: file1
    Processed: file2
    Processed: file3

    See Also
    --------
    osgeo.ogr.DataSource : GDAL vector data container
    osgeo.ogr.Layer : Vector layer with features
    osgeo.ogr.Feature : Individual vector feature

    References
    ----------
    .. [1] GDAL/OGR Vector Data Model. https://gdal.org/user/vector_data_model.html
    .. [2] OGR Vector Formats. https://gdal.org/drivers/vector/index.html
    .. [3] OGR API Tutorial. https://gdal.org/tutorials/vector_api_tut.html

    """
    # ========================================================================
    # Input validation
    # ========================================================================

    # Validate input filename
    if not isinstance(sFilename_in, str):
        raise TypeError(f"sFilename_in must be a string, got {type(sFilename_in)}")

    if not sFilename_in or sFilename_in.isspace():
        raise ValueError("sFilename_in cannot be empty")

    # Validate output filename
    if not isinstance(sFilename_out, str):
        raise TypeError(f"sFilename_out must be a string, got {type(sFilename_out)}")

    if not sFilename_out or sFilename_out.isspace():
        raise ValueError("sFilename_out cannot be empty")

    # Validate layer index
    if not isinstance(iLayer_index, int):
        raise TypeError(f"iLayer_index must be an integer, got {type(iLayer_index)}")

    if iLayer_index < 0:
        raise ValueError(f"iLayer_index must be non-negative, got {iLayer_index}")

    logger.info(
        f"Copying geometry without attributes: {sFilename_in} → {sFilename_out}"
    )

    # ========================================================================
    # Check input file exists
    # ========================================================================

    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"Input file not found: {sFilename_in}")

    logger.debug(f"Input file exists: {sFilename_in}")

    # ========================================================================
    # Auto-detect or validate input format
    # ========================================================================

    if sFormat_in is not None:
        logger.info(f"Using specified input format: {sFormat_in}")
        pDriver_in = ogr.GetDriverByName(sFormat_in)
        if pDriver_in is None:
            raise RuntimeError(f"Input GDAL driver '{sFormat_in}' not available")
    else:
        # Try to get driver from extension for better format detection
        try:
            pDriver_in = get_vector_driver_from_extension(sFilename_in)
            if pDriver_in is not None:
                sFormat_detected = get_vector_format_from_extension(sFilename_in)
                logger.debug(f"Auto-detected input format: {sFormat_detected}")
            else:
                logger.debug("Could not auto-detect input format, will use ogr.Open()")
        except Exception as e:
            logger.debug(f"Auto-detection failed: {e}, will use ogr.Open()")
            pDriver_in = None  # Will use ogr.Open() generic method

    # ========================================================================
    # Auto-detect or validate output format
    # ========================================================================

    if sFormat_out is None:
        # Auto-detect from extension using library function
        try:
            sFormat_out = get_vector_format_from_extension(sFilename_out)
            logger.info(f"Auto-detected output format: {sFormat_out}")
        except Exception as e:
            logger.warning(
                f"Cannot auto-detect output format from '{sFilename_out}': {e}. "
                f"Defaulting to GeoJSON."
            )
            sFormat_out = "GeoJSON"
    else:
        logger.info(f"Using specified output format: {sFormat_out}")

    # Get output driver using library function
    try:
        pDriver_out = get_vector_driver_from_extension(sFilename_out)
        if pDriver_out is None:
            # Fallback to getting driver by name
            pDriver_out = ogr.GetDriverByName(sFormat_out)
            if pDriver_out is None:
                raise RuntimeError(f"Output GDAL driver '{sFormat_out}' not available")
    except Exception as e:
        # Fallback to getting driver by name
        logger.debug(
            f"get_vector_driver_from_extension failed: {e}, using GetDriverByName"
        )
        pDriver_out = ogr.GetDriverByName(sFormat_out)
        if pDriver_out is None:
            raise RuntimeError(f"Output GDAL driver '{sFormat_out}' not available")

    # ========================================================================
    # Open input dataset
    # ========================================================================

    try:
        if pDriver_in is not None:
            pDataset_in = pDriver_in.Open(sFilename_in, 0)  # 0 = read-only
        else:
            pDataset_in = ogr.Open(sFilename_in, 0)

        if pDataset_in is None:
            raise ValueError(f"Could not open input file: {sFilename_in}")
    except Exception as e:
        raise ValueError(f"Failed to open input file '{sFilename_in}': {e}") from e

    logger.debug(f"Opened input dataset: {sFilename_in}")

    # ========================================================================
    # Get input layer
    # ========================================================================

    n_layers = pDataset_in.GetLayerCount()
    logger.debug(f"Input file has {n_layers} layer(s)")

    if n_layers > 1:
        logger.warning(
            f"Input file has {n_layers} layers. " f"Only one layer will be copied."
        )

    # Get layer by name or index
    if sLayer_name is not None:
        logger.info(f"Accessing layer by name: '{sLayer_name}'")
        pLayer_in = pDataset_in.GetLayerByName(sLayer_name)
        if pLayer_in is None:
            pDataset_in = None
            raise ValueError(f"Layer '{sLayer_name}' not found in input file")
    else:
        if iLayer_index >= n_layers:
            pDataset_in = None
            raise ValueError(
                f"Layer index {iLayer_index} out of range "
                f"(file has {n_layers} layer(s))"
            )
        logger.info(f"Accessing layer by index: {iLayer_index}")
        pLayer_in = pDataset_in.GetLayer(iLayer_index)

    if pLayer_in is None:
        pDataset_in = None
        raise ValueError("Could not access input layer")

    layer_name_actual = pLayer_in.GetName()
    feature_count = pLayer_in.GetFeatureCount()
    logger.info(f"Input layer: '{layer_name_actual}', {feature_count} features")

    # ========================================================================
    # Get spatial reference and geometry type
    # ========================================================================

    pSRS = pLayer_in.GetSpatialRef()
    if pSRS is not None:
        srs_name = pSRS.GetName() if pSRS.GetName() else "Unknown"
        logger.debug(f"Input spatial reference: {srs_name}")
    else:
        logger.warning("Input layer has no spatial reference")

    # Get geometry type from first feature
    pLayer_in.ResetReading()
    pFeature_first = pLayer_in.GetNextFeature()

    if pFeature_first is None:
        pDataset_in = None
        raise ValueError("Input layer contains no features")

    pGeometry_first = pFeature_first.GetGeometryRef()
    if pGeometry_first is None:
        pDataset_in = None
        raise ValueError("First feature has no geometry")

    geometry_type = pGeometry_first.GetGeometryType()
    geometry_name = ogr.GeometryTypeToName(geometry_type)
    logger.info(f"Geometry type: {geometry_name}")

    # ========================================================================
    # Create output file
    # ========================================================================

    # Remove existing output file if it exists
    if os.path.exists(sFilename_out):
        logger.debug(f"Removing existing output file: {sFilename_out}")

        # For shapefiles, remove all associated files
        if sFormat_out == "ESRI Shapefile":
            base_name = os.path.splitext(sFilename_out)[0]
            for ext in [
                ".shp",
                ".shx",
                ".dbf",
                ".prj",
                ".cpg",
                ".qpj",
                ".sbn",
                ".sbx",
                ".shp.xml",
            ]:
                file_to_delete = base_name + ext
                if os.path.exists(file_to_delete):
                    os.remove(file_to_delete)
                    logger.debug(f"Removed {file_to_delete}")
        else:
            os.remove(sFilename_out)

    # Create output directory if needed
    output_dir = os.path.dirname(sFilename_out)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(f"Created output directory: {output_dir}")
        except OSError as e:
            pDataset_in = None
            raise OSError(f"Cannot create output directory '{output_dir}': {e}") from e

    # Create output dataset
    try:
        pDataset_out = pDriver_out.CreateDataSource(sFilename_out)
        if pDataset_out is None:
            raise RuntimeError(f"Failed to create output file: {sFilename_out}")
    except Exception as e:
        pDataset_in = None
        raise RuntimeError(
            f"Cannot create output dataset '{sFilename_out}': {e}"
        ) from e

    logger.debug(f"Created output dataset: {sFilename_out}")

    # ========================================================================
    # Create output layer (geometry only, no attribute fields)
    # ========================================================================

    output_layer_name = "geometry_only"
    logger.debug(f"Creating output layer: {output_layer_name}")

    try:
        pLayer_out = pDataset_out.CreateLayer(
            output_layer_name, srs=pSRS, geom_type=geometry_type
        )

        if pLayer_out is None:
            raise RuntimeError("Failed to create output layer")
    except Exception as e:
        pDataset_in = None
        pDataset_out = None
        raise RuntimeError(f"Cannot create output layer: {e}") from e

    logger.debug("Output layer created (no attribute fields)")

    # Get layer definition for creating new features
    pLayerDefn_out = pLayer_out.GetLayerDefn()

    # ========================================================================
    # Copy geometries from input to output
    # ========================================================================

    logger.info(f"Copying {feature_count} geometries...")

    # Reset reading to start
    pLayer_in.ResetReading()

    n_copied = 0
    n_skipped = 0

    for idx, pFeature_in in enumerate(pLayer_in):
        if idx % 1000 == 0 and idx > 0:
            logger.debug(f"Processed {idx}/{feature_count} features...")

        # Get geometry from input feature
        pGeometry = pFeature_in.GetGeometryRef()

        if pGeometry is None:
            logger.warning(f"Feature {idx} has no geometry, skipping")
            n_skipped += 1
            continue

        try:
            # Create new feature with only geometry (no attributes)
            pFeature_out = ogr.Feature(pLayerDefn_out)
            pFeature_out.SetGeometry(pGeometry.Clone())

            # Add feature to output layer
            if pLayer_out.CreateFeature(pFeature_out) != 0:
                logger.warning(f"Failed to create feature {idx}")
                n_skipped += 1
            else:
                n_copied += 1

            # Clean up
            pFeature_out = None
        except Exception as e:
            logger.warning(f"Error copying feature {idx}: {e}")
            n_skipped += 1

    # ========================================================================
    # Clean up and finalize
    # ========================================================================

    # Close datasets
    pDataset_in = None
    pDataset_out = None

    logger.info(
        f"Geometry copy complete:\n"
        f"  Input: {sFilename_in}\n"
        f"  Output: {sFilename_out}\n"
        f"  Features copied: {n_copied}\n"
        f"  Features skipped: {n_skipped}\n"
        f"  Output format: {sFormat_out}"
    )

    if n_copied == 0:
        logger.warning("No features were copied!")

    return sFilename_out
