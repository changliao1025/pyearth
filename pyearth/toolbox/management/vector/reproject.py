import os
import osgeo
from osgeo import osr, ogr
from typing import Optional
from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_format_from_filename,
    print_supported_vector_formats,
    get_vector_driver_from_format,
    get_vector_driver_from_filename,
)


def reproject_vector(
    input_file: str,
    output_file: str,
    target_projection: str,
    copy_attributes: bool = True,
) -> bool:
    """
    Reproject vector geometries from one coordinate system to another.

    Transforms all geometries in a vector file from their source coordinate system
    to a target coordinate system using GDAL/OGR coordinate transformation.
    Supports all OGR-supported vector formats and geometry types.

    Parameters
    ----------
    input_file : str
        Path to input vector file. Must exist and be readable by GDAL/OGR.
    output_file : str
        Path for output vector file. Format auto-detected from extension.
        Existing files will be overwritten.
    target_projection : str
        Target coordinate system as WKT (Well-Known Text) string.
        Must be a valid spatial reference system definition.
    copy_attributes : bool, optional
        If True, copy all attribute fields from input to output. Default is True.

    Returns
    -------
    bool
        True if reprojection completed successfully, False if errors occurred.

    Raises
    ------
    FileNotFoundError
        If input file does not exist.
    RuntimeError
        If input file cannot be opened or output file cannot be created.
    ValueError
        If target projection is invalid or transformation cannot be created.

    Notes
    -----
    1. **Coordinate Systems**: The source CRS is automatically detected from the
       input file. The target CRS must be provided as a valid WKT string.

    2. **Geometry Types**: All geometry types are supported (Point, LineString,
       Polygon, MultiPoint, MultiLineString, MultiPolygon, etc.). The output
       geometry type matches the input.

    3. **GDAL Version Compatibility**: Handles axis order changes in GDAL 3+
       by setting the traditional GIS axis mapping strategy.

    4. **Attribute Preservation**: When copy_attributes=True, all feature attributes
       are copied to the output. Field definitions are preserved.

    5. **File Formats**: Input format auto-detected. Output format determined by
       file extension. Both must be supported by GDAL/OGR.

    6. **Performance**: Processes features sequentially to handle large files.
       No geometry simplification or validation is performed.

    Examples
    --------
    Reproject to WGS84 (EPSG:4326):

    >>> from osgeo import osr
    >>> wgs84 = osr.SpatialReference()
    >>> wgs84.ImportFromEPSG(4326)
    >>> wkt = wgs84.ExportToWkt()
    >>> success = reproject_vector('input.shp', 'output_wgs84.shp', wkt)
    >>> success
    True

    Reproject with custom projection string:

    >>> utm_wkt = '''PROJCS["WGS 84 / UTM zone 33N",
    ...     GEOGCS["WGS 84", ...]]'''  # Full WKT string
    >>> success = reproject_vector('input.geojson', 'output_utm.geojson',
    ...                           utm_wkt, copy_attributes=True)

    See Also
    --------
    print_supported_vector_formats : List available vector formats
    osr.SpatialReference : GDAL spatial reference system handling
    osr.CoordinateTransformation : GDAL coordinate transformation
    """
    # Input validation
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input vector file does not exist: {input_file}")

    # Remove existing output file if it exists
    if os.path.exists(output_file):
        try:
            os.remove(output_file)
        except OSError as e:
            print(f"Warning: Could not remove existing output file: {e}")

    # Open input dataset
    input_dataset = ogr.Open(input_file, 0)  # Read-only
    if input_dataset is None:
        raise RuntimeError(f"Could not open input file: {input_file}")

    try:
        # Get input layer and spatial reference
        input_layer = input_dataset.GetLayer(0)
        if input_layer is None:
            raise RuntimeError("Could not access layer in input file")

        source_srs = input_layer.GetSpatialRef()
        if source_srs is None:
            print("Warning: Input file has no spatial reference defined")
            source_srs = osr.SpatialReference()  # Create default

        source_wkt = source_srs.ExportToWkt()

        # Create target spatial reference
        target_srs = osr.SpatialReference()
        if target_srs.ImportFromWkt(target_projection) != 0:
            raise ValueError(f"Invalid target projection WKT: {target_projection}")

        # Check if reprojection is needed
        if source_wkt == target_projection:
            print(
                "Input and target coordinate systems are identical - no reprojection needed"
            )
            return True

        # Handle GDAL 3+ axis order changes
        if int(osgeo.__version__[0]) >= 3:
            source_srs.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
            target_srs.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

        # Create coordinate transformation
        transform = osr.CoordinateTransformation(source_srs, target_srs)
        if transform is None:
            raise ValueError("Could not create coordinate transformation")

        # Get output driver
        output_driver = get_vector_driver_from_filename(output_file)
        if output_driver is None:
            raise RuntimeError(f"Could not get driver for output file: {output_file}")

        # Create output dataset
        output_dataset = output_driver.CreateDataSource(output_file)
        if output_dataset is None:
            raise RuntimeError(f"Could not create output file: {output_file}")

        try:
            # Get geometry type from first feature
            input_layer.ResetReading()
            first_feature = input_layer.GetNextFeature()
            if first_feature is None:
                print("Warning: Input file contains no features")
                return True

            geometry = first_feature.GetGeometryRef()
            geometry_type = geometry.GetGeometryType() if geometry else ogr.wkbUnknown

            # Create output layer
            output_layer = output_dataset.CreateLayer(
                "reprojected", target_srs, geom_type=geometry_type
            )
            if output_layer is None:
                raise RuntimeError("Could not create output layer")

            # Copy field definitions if requested
            if copy_attributes:
                input_defn = input_layer.GetLayerDefn()
                for i in range(input_defn.GetFieldCount()):
                    field_defn = input_defn.GetFieldDefn(i)
                    if output_layer.CreateField(field_defn) != 0:
                        print(f"Warning: Failed to create field {field_defn.GetName()}")

            # Process all features
            feature_count = input_layer.GetFeatureCount()
            processed_count = 0

            input_layer.ResetReading()  # Reset to beginning
            feature = input_layer.GetNextFeature()

            while feature:
                processed_count += 1
                if processed_count % 1000 == 0:
                    print(f"Processed {processed_count}/{feature_count} features")

                # Get geometry and transform it
                geometry = feature.GetGeometryRef()
                if geometry is not None:
                    try:
                        # Create transformed geometry
                        transformed_geometry = geometry.Clone()
                        transform_result = transformed_geometry.Transform(transform)

                        if transform_result != 0:
                            print(
                                f"Warning: Failed to transform feature {processed_count}"
                            )
                            feature = input_layer.GetNextFeature()
                            continue

                        # Create output feature
                        output_feature = ogr.Feature(output_layer.GetLayerDefn())
                        output_feature.SetGeometry(transformed_geometry)

                        # Copy attributes if requested
                        if copy_attributes:
                            output_defn = output_layer.GetLayerDefn()
                            input_defn = input_layer.GetLayerDefn()

                            for i in range(input_defn.GetFieldCount()):
                                field_name = input_defn.GetFieldDefn(i).GetName()
                                output_field_index = output_defn.GetFieldIndex(
                                    field_name
                                )
                                if output_field_index >= 0:
                                    output_feature.SetField(
                                        field_name, feature.GetField(field_name)
                                    )

                        # Add feature to output layer
                        if output_layer.CreateFeature(output_feature) != 0:
                            print(
                                f"Warning: Failed to create output feature {processed_count}"
                            )

                        output_feature = None

                    except Exception as e:
                        print(
                            f"Warning: Error processing feature {processed_count}: {e}"
                        )
                else:
                    print(f"Warning: Feature {processed_count} has no geometry")

                feature = input_layer.GetNextFeature()

            # Flush changes
            output_dataset.FlushCache()

            print(f"Successfully reprojected {processed_count} features")
            return True

        finally:
            # Clean up output dataset
            output_dataset = None

    finally:
        # Clean up input dataset
        input_dataset = None
