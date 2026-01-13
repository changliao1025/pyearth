"""
Raster to Vector Conversion (Polygonization)

This module provides functionality for converting raster datasets into vector polygons,
a process also known as polygonization or vectorization. It extracts contiguous regions
of equal pixel values and converts them into polygon features with associated attributes.

Main Functions
--------------
vectorize_raster : Convert raster data to vector polygons with dominant value extraction

Key Features
------------
- Automatic polygonization of raster cells into contiguous regions
- NoData value masking to exclude background pixels
- Dominant value calculation for each polygon from original raster
- Custom field name and type specification
- Sequential polygon ID assignment
- Multi-format vector output (Shapefile, GeoJSON, GeoPackage, etc.)
- Support for integer and floating-point rasters
- Memory-efficient processing using GDAL's Polygonize function

Use Cases
---------
1. **Classification Maps**: Convert classified raster to polygon land use zones
2. **Contour Extraction**: Generate polygon boundaries from continuous surfaces
3. **Mask Vectorization**: Convert binary masks to polygon features
4. **Region Extraction**: Identify and extract homogeneous regions
5. **Segmentation Output**: Convert image segmentation results to vectors
6. **Change Detection**: Polygonize areas of change from raster analysis
7. **Model Output**: Convert gridded model results to GIS-compatible polygons

Technical Details
-----------------
The module uses GDAL's Polygonize function to convert raster cells into polygon
features. The process involves:

1. **Masking**: Create binary mask excluding NoData pixels
2. **Polygonization**: Group contiguous cells with same value into polygons
3. **Value Extraction**: For each polygon, extract dominant (most frequent) value
   from original raster using clipping and statistical analysis
4. **Attribute Assignment**: Store polygon ID and dominant value in output fields

Dominant Value Calculation:
For each polygon, the original raster is clipped to the polygon's extent and the
most frequent pixel value is computed using bincount (integers) or unique value
counting (floats). This handles cases where polygonization merges slightly
different values or when resampling occurs.

Performance Characteristics
---------------------------
- Time Complexity: O(N*P) where N=pixels, P=polygons (value extraction per polygon)
- Space Complexity: O(P) where P=number of polygons
- Memory Usage: Mask array + polygon clipping requires ~2x raster size in memory
- Large Datasets: May generate thousands of small polygons requiring cleanup

Optimization Tips:
- Use appropriate NoData values to exclude background
- Consider raster smoothing before vectorization to reduce polygon count
- Apply minimum area threshold post-processing to remove tiny polygons
- Use integer rasters when possible (faster dominant value calculation)

Dependencies
------------
- GDAL/OGR: Rasterization and vector I/O
- NumPy: Array operations and statistical calculations
- OSR: Spatial reference system handling

See Also
--------
- rasterize_vector: Inverse operation (vector to raster)
- gdal.Polygonize: GDAL's core polygonization function
- numpy_to_gdal_type: Data type conversion utility
"""

import os
import logging
from typing import Optional
import numpy as np
from osgeo import gdal, ogr, osr

from pyearth.gis.gdal.gdal_to_numpy_datatype import numpy_to_gdal_type
from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_driver_from_extension,
    get_vector_format_from_extension,
)

# Configure logging
logger = logging.getLogger(__name__)


def vectorize_raster(
    sFilename_raster_in: str,
    sFilename_vector_out: str,
    sFieldname: str = "value",
    sFieldtype: int = ogr.OFTInteger64,
) -> bool:
    """
    Convert raster dataset to vector polygons with dominant value extraction.

    This function polygonizes a raster dataset by grouping contiguous pixels with
    similar values into polygon features. Each polygon is assigned a unique ID and
    a value field containing the dominant (most frequent) pixel value from the
    original raster within that polygon's extent.

    The function creates a mask to exclude NoData pixels, performs polygonization
    using GDAL's Polygonize algorithm, then extracts the most common value within
    each polygon through statistical analysis of the clipped raster data.

    Parameters
    ----------
    sFilename_raster_in : str
        Absolute path to input raster file. Supports all GDAL-readable formats
        (GeoTIFF, IMG, HDF, etc.). Must have defined NoData value for proper
        masking. Single-band rasters only.
    sFilename_vector_out : str
        Absolute path for output vector file. Format automatically detected from
        file extension. Supports multiple formats:
        - Shapefile: .shp (includes .shx, .dbf, .prj components)
        - GeoJSON: .geojson, .json
        - GeoPackage: .gpkg
        - GML: .gml
        - KML: .kml
        Existing files will be deleted and recreated.
    sFieldname : str, optional
        Name of the attribute field to store dominant raster values. Default is
        'value'. Field name must comply with Shapefile restrictions (max 10 chars,
        no spaces or special characters except underscore).
    sFieldtype : int, optional
        OGR field type for the value attribute. Default is ogr.OFTInteger64.
        Common options:
        - ogr.OFTInteger: 32-bit integer
        - ogr.OFTInteger64: 64-bit integer (default)
        - ogr.OFTReal: Floating-point (for continuous rasters)
        Must match the raster data type to avoid precision loss.

    Returns
    -------
    bool
        True if vectorization completed successfully, False if errors occurred
        (e.g., file not found, invalid raster format, polygonization failure).

    Raises
    ------
    FileNotFoundError
        If input raster file does not exist (implicit via GDAL).
    RuntimeError
        If raster cannot be opened or polygonization fails (implicit via GDAL).

    Notes
    -----
    1. **NoData Handling**: The function creates a binary mask where pixels equal
       to the raster's NoData value are excluded from polygonization. Only valid
       data pixels are converted to polygons. Rasters without NoData value may
       include unwanted background polygons.

    2. **Polygonization Algorithm**: Uses GDAL's Polygonize function which groups
       4-connected or 8-connected pixels (implementation dependent) into polygons.
       Small variations in pixel values may create many small polygons. Consider
       reclassification or smoothing before vectorization.

    3. **Dominant Value Extraction**: For each polygon, the function:
       - Clips the original raster to the polygon's bounding box using Warp
       - Masks pixels outside the polygon geometry
       - Computes most frequent value using bincount (integers) or unique counts (floats)
       This handles edge cases where polygon boundaries don't align perfectly with
       pixels or when resampling creates slight value variations.

    4. **Output Schema**: Each polygon feature contains two attributes:
       - 'id': Sequential integer starting from 1 (unique polygon identifier)
       - <sFieldname>: Dominant pixel value (type specified by sFieldtype)

    5. **Spatial Reference**: Output inherits the spatial reference system from
       the input raster. No coordinate transformation is performed.

    6. **Memory Usage**: The function creates:
       - Full raster mask array in memory (~size of input raster)
       - Temporary clip for each polygon (typically small)
       For very large rasters (>4GB), may require substantial memory.

    7. **Performance Considerations**: Processing time grows with polygon count.
       For rasters producing thousands of tiny polygons, consider:
       - Increasing minimum mapping unit through raster preprocessing
       - Post-processing to remove small polygons
       - Using higher resolution for initial classification

    8. **Data Type Compatibility**: Ensure sFieldtype matches raster data:
       - Integer rasters → ogr.OFTInteger or ogr.OFTInteger64
       - Float rasters → ogr.OFTReal
       Mismatches may cause precision loss or errors.

    9. **Output Format Limitations**: Different formats have different constraints:
       - Shapefile: Field names max 10 characters, ~2GB file size limit
       - GeoJSON: No field name limits, but large files can be slow to parse
       - GeoPackage: No practical limits, good for large datasets
       - Consider format choice based on data size and field naming needs.

    10. **Edge Effects**: Polygon boundaries follow pixel edges, creating
        "staircase" effects. For smooth boundaries, consider:
        - Higher resolution input rasters
        - Post-processing simplification
        - Alternative smoothing algorithms

    Examples
    --------
    Basic land cover classification vectorization (Shapefile):

    >>> success = vectorize_raster(
    ...     sFilename_raster_in='/data/landcover_classified.tif',
    ...     sFilename_vector_out='/output/landcover_polygons.shp',
    ...     sFieldname='lc_code',
    ...     sFieldtype=ogr.OFTInteger
    ... )
    >>> print(f"Success: {success}")
    # Output: Shapefile with polygons, each with 'id' and 'lc_code' fields
    # 'lc_code' contains the land cover classification value

    Vectorize continuous elevation data to GeoJSON:

    >>> success = vectorize_raster(
    ...     sFilename_raster_in='/data/elevation.tif',
    ...     sFilename_vector_out='/output/elevation_zones.geojson',
    ...     sFieldname='elev_m',
    ...     sFieldtype=ogr.OFTReal
    ... )
    # Output: GeoJSON with polygons containing average elevation values

    Convert binary mask to GeoPackage polygons:

    >>> success = vectorize_raster(
    ...     sFilename_raster_in='/data/water_mask.tif',
    ...     sFilename_vector_out='/output/water_bodies.gpkg',
    ...     sFieldname='water',
    ...     sFieldtype=ogr.OFTInteger
    ... )
    # Output: GeoPackage with water body polygons, background (NoData) excluded

    Expected processing workflow:
    ```
    INFO: Processing raster: /data/landcover_classified.tif
    INFO: Raster dimensions: 5000 x 4000 pixels
    INFO: NoData value: -9999
    INFO: Creating mask for valid pixels
    INFO: Polygonizing raster...
    INFO: Created 1247 preliminary polygons
    INFO: Extracting dominant values for each polygon...
    INFO: Processing polygon 100/1247
    INFO: Processing polygon 200/1247
    ...
    INFO: Processing polygon 1200/1247
    INFO: Vectorization completed successfully
    INFO: Output: /output/landcover_polygons.shp (1247 features)
    ```

    See Also
    --------
    rasterize_vector : Inverse operation (convert vector to raster)
    gdal.Polygonize : GDAL's polygonization function
    numpy_to_gdal_type : Convert NumPy values to GDAL-compatible types
    """
    # Validate input file exists
    if not os.path.exists(sFilename_raster_in):
        logger.error(f"Input raster file not found: {sFilename_raster_in}")
        return False

    logger.info(f"Processing raster: {sFilename_raster_in}")

    # Get vector driver based on output file extension
    try:
        pDriver_vector = get_vector_driver_from_extension(sFilename_vector_out)
        sFormat_out = get_vector_format_from_extension(sFilename_vector_out)
        logger.info(f"Output format: {sFormat_out}")
    except ValueError as e:
        logger.error(f"Unsupported output format: {e}")
        return False

    if pDriver_vector is None:
        logger.error(f"Could not load driver for output format: {sFormat_out}")
        return False

    # Open input raster dataset
    pDataset_raster = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)
    if pDataset_raster is None:
        logger.error(f"Could not open raster file: {sFilename_raster_in}")
        return False

    # Log raster information
    ncolumn = pDataset_raster.RasterXSize
    nrow = pDataset_raster.RasterYSize
    logger.info(f"Raster dimensions: {ncolumn} x {nrow} pixels")

    # Delete existing output file if present
    if os.path.exists(sFilename_vector_out):
        # Handle format-specific cleanup
        if sFormat_out == "ESRI Shapefile":
            # Delete all shapefile component files
            base_name = os.path.splitext(sFilename_vector_out)[0]
            for ext in [".shp", ".shx", ".dbf", ".prj", ".cpg", ".shp.xml"]:
                component_file = base_name + ext
                if os.path.exists(component_file):
                    os.remove(component_file)
            logger.info(
                f"Removed existing shapefile components: {sFilename_vector_out}"
            )
        else:
            # For other formats, use driver's delete method
            pDriver_vector.DeleteDataSource(sFilename_vector_out)
            logger.info(f"Removed existing output file: {sFilename_vector_out}")

    # Extract geotransform parameters
    pGeotransform = pDataset_raster.GetGeoTransform()
    dOriginX = pGeotransform[0]
    dOriginY = pGeotransform[3]
    dPixelWidth = pGeotransform[1]
    dPixelHeight = pGeotransform[5]

    # Get first raster band
    pBand = pDataset_raster.GetRasterBand(1)
    if pBand is None:
        logger.error("Could not access raster band")
        pDataset_raster = None
        return False

    # Get spatial reference from raster
    pProjection_source = pDataset_raster.GetProjection()
    spatialRef = osr.SpatialReference()
    if pProjection_source:
        spatialRef.ImportFromWkt(pProjection_source)
    else:
        logger.warning("Input raster has no spatial reference defined")

    # Get NoData value and data type
    nodata_value = pBand.GetNoDataValue()
    if nodata_value is None:
        logger.warning(
            "Raster has no NoData value defined - all pixels will be polygonized"
        )
        nodata_value = -9999  # Default fallback

    logger.info(f"NoData value: {nodata_value}")
    eType_out = pBand.DataType

    # Create in-memory mask dataset
    logger.info("Creating mask for valid pixels")
    mask_ds = gdal.GetDriverByName("MEM").Create(
        "", pDataset_raster.RasterXSize, pDataset_raster.RasterYSize, 1, gdal.GDT_Byte
    )
    if mask_ds is None:
        logger.error("Could not create memory dataset for mask")
        pDataset_raster = None
        return False

    mask_band = mask_ds.GetRasterBand(1)

    # Create binary mask: 1 for valid pixels, 0 for NoData
    try:
        raster_array = pBand.ReadAsArray()
        if raster_array is None:
            logger.error("Could not read raster data into array")
            pDataset_raster = None
            mask_ds = None
            return False

        mask_array = np.not_equal(raster_array, nodata_value).astype(np.uint8)
        mask_band.WriteArray(mask_array)
        mask_band.FlushCache()

        valid_pixel_count = np.sum(mask_array)
        total_pixels = mask_array.size
        logger.info(
            f"Valid pixels: {valid_pixel_count}/{total_pixels} ({100*valid_pixel_count/total_pixels:.1f}%)"
        )
    except Exception as e:
        logger.error(f"Error creating mask array: {e}")
        pDataset_raster = None
        mask_ds = None
        return False

    # Create output vector dataset
    logger.info(f"Creating output {sFormat_out} dataset")
    pDataset_vector = pDriver_vector.CreateDataSource(sFilename_vector_out)
    if pDataset_vector is None:
        logger.error(f"Could not create output file: {sFilename_vector_out}")
        pDataset_raster = None
        mask_ds = None
        return False

    pLayer_out = pDataset_vector.CreateLayer("polygonized", spatialRef, ogr.wkbPolygon)
    if pLayer_out is None:
        logger.error("Could not create output layer")
        pDataset_raster = None
        mask_ds = None
        pDataset_vector = None
        return False

    # Create ID field
    pField_defn_id = ogr.FieldDefn("id", ogr.OFTInteger64)
    pLayer_out.CreateField(pField_defn_id)

    # Create value field with specified type
    pField_defn_value = ogr.FieldDefn(sFieldname, sFieldtype)
    pLayer_out.CreateField(pField_defn_value)

    logger.info(
        f"Output fields created: 'id' (Integer64), '{sFieldname}' ({sFieldtype})"
    )

    # Perform polygonization
    logger.info("Polygonizing raster...")
    field_index = 0  # First field for polygon values (temporary, will be replaced)

    def polygonize_callback(progress, message, callback_data):
        """Callback function for polygonization progress monitoring."""
        return 1  # Return 1 to continue, 0 to abort

    try:
        result = gdal.Polygonize(
            pBand, mask_band, pLayer_out, field_index, [], callback=polygonize_callback
        )

        if result != 0:
            logger.error("Polygonization failed")
            pDataset_raster = None
            mask_ds = None
            pDataset_vector = None
            return False
    except Exception as e:
        logger.error(f"Error during polygonization: {e}")
        pDataset_raster = None
        mask_ds = None
        pDataset_vector = None
        return False

    # Get polygon count
    nFeature_count = pLayer_out.GetFeatureCount()
    logger.info(f"Created {nFeature_count} preliminary polygons")

    # Extract dominant values for each polygon
    logger.info("Extracting dominant values for each polygon...")
    pLayer_out.ResetReading()
    pFeature = pLayer_out.GetNextFeature()
    lID = 1

    while pFeature:
        # Log progress periodically
        if lID % 100 == 0:
            logger.info(f"Processing polygon {lID}/{nFeature_count}")

        # Get polygon geometry
        geom = pFeature.GetGeometryRef()
        if geom is None:
            logger.warning(f"Polygon {lID} has no geometry, skipping")
            pFeature = pLayer_out.GetNextFeature()
            lID += 1
            continue

        pPolygonWKT = geom.ExportToWkt()
        aPolygon_extent = geom.GetEnvelope()
        minX, maxX, minY, maxY = aPolygon_extent
        aPolygon_extent = [minX, minY, maxX, maxY]

        # Clip raster to polygon extent
        try:
            pWrapOption = gdal.WarpOptions(
                cropToCutline=False,
                cutlineWKT=pPolygonWKT,
                xRes=dPixelWidth,
                yRes=abs(dPixelHeight),
                outputBounds=aPolygon_extent,
                dstSRS=spatialRef,
                format="MEM",
                resampleAlg="near",
                dstNodata=nodata_value,
                outputType=eType_out,
            )

            pDataset_clip_warped = gdal.Warp(
                "", sFilename_raster_in, options=pWrapOption
            )
            if pDataset_clip_warped is None:
                logger.warning(f"Could not clip raster for polygon {lID}, using NoData")
                dominant_value = nodata_value
            else:
                # Convert clipped raster to array
                aData_clip = pDataset_clip_warped.ReadAsArray()

                # Calculate dominant (most frequent) value
                if aData_clip is not None and aData_clip.size > 0:
                    # Flatten array
                    valid_data = aData_clip.flatten()

                    # Remove NoData values
                    if nodata_value is not None:
                        valid_data = valid_data[valid_data != nodata_value]

                    # Find most frequent value
                    if valid_data.size > 0:
                        if np.issubdtype(valid_data.dtype, np.integer):
                            # For integer arrays, use bincount (faster)
                            dominant_value = np.bincount(valid_data).argmax()
                        else:
                            # For float arrays, use unique value counting
                            unique_values, counts = np.unique(
                                valid_data, return_counts=True
                            )
                            dominant_value = unique_values[np.argmax(counts)]
                    else:
                        dominant_value = nodata_value
                else:
                    dominant_value = nodata_value

                pDataset_clip_warped = None
        except Exception as e:
            logger.warning(f"Error extracting value for polygon {lID}: {e}")
            dominant_value = nodata_value

        # Set polygon ID
        pFeature.SetField("id", lID)

        # Convert and set dominant value
        try:
            dominant_value_gdal = numpy_to_gdal_type(dominant_value, sFieldtype)
            pFeature.SetField(sFieldname, dominant_value_gdal)
        except Exception as e:
            logger.warning(f"Error setting field value for polygon {lID}: {e}")

        # Update feature in layer
        pLayer_out.SetFeature(pFeature)

        # Move to next polygon
        lID += 1
        pFeature = pLayer_out.GetNextFeature()

    # Cleanup resources
    pDataset_vector = None
    pDataset_raster = None
    mask_ds = None

    logger.info("Vectorization completed successfully")
    logger.info(f"Output: {sFilename_vector_out} ({nFeature_count} features)")

    return True


if __name__ == "__main__":
    pass
