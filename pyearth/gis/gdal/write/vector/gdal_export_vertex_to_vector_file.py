import os
from typing import List, Optional, Any

from osgeo import ogr, osr
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_format_from_extension, get_vector_driver_from_extension


def export_vertex_to_vector_file(
    aVertex_in: List[Any],
    sFilename_vector_in: str,
    iFlag_projected_in: Optional[int] = None,
    pSpatial_reference_in: Optional[osr.SpatialReference] = None,
    aAttribute_data: Optional[List[int]] = None,
    overwrite: bool = True
) -> None:
    """Export vertices as points to a vector file.

    Parameters
    ----------
    aVertex_in : list
        List of vertex objects with coordinate information. Each vertex must have:
        - For geographic (iFlag_projected_in=0): dLongitude_degree, dLatitude_degree
        - For projected (iFlag_projected_in=1): dx, dy
    sFilename_vector_in : str
        Output filename with extension. Format is determined from extension.
        Supports GeoJSON, Shapefile, GeoPackage, KML, etc.
    iFlag_projected_in : int, optional
        Coordinate system flag: 0 for geographic, 1 for projected.
        Default is 0 (geographic).
    pSpatial_reference_in : osr.SpatialReference, optional
        Spatial reference system. If None, defaults to WGS84 (EPSG:4326).
    aAttribute_data : list of int, optional
        Optional connectivity attribute data for each vertex.
    overwrite : bool, optional
        If True, overwrite existing file. Default is True.

    Raises
    ------
    ValueError
        If aVertex_in is empty or if file extension is not supported.
    FileExistsError
        If the output file exists and overwrite is False.
    RuntimeError
        If dataset, layer, or feature creation fails.

    Notes
    -----
    - Creates a point layer with 'pointid' field (1-based index)
    - If aAttribute_data is provided, adds 'connectivity' field
    - Automatically determines output format from file extension
    """

    if not aVertex_in:
        raise ValueError("Vertex list cannot be empty.")

    # Handle overwrite
    if os.path.exists(sFilename_vector_in):
        if not overwrite:
            raise FileExistsError(f"Output file already exists: {sFilename_vector_in}")
        os.remove(sFilename_vector_in)

    # Set defaults
    iFlag_projected = 0 if iFlag_projected_in is None else iFlag_projected_in

    if pSpatial_reference_in is None:
        spatial_ref = osr.SpatialReference()
        spatial_ref.ImportFromEPSG(4326)  # WGS84 lat/lon
    else:
        spatial_ref = pSpatial_reference_in

    # Get driver from file extension
    driver = get_vector_driver_from_extension(sFilename_vector_in)

    iFlag_attribute = aAttribute_data is not None
    attribute_data = aAttribute_data if aAttribute_data is not None else []

    dataset = None
    layer = None
    feature = None

    try:
        # Create dataset and layer
        dataset = driver.CreateDataSource(sFilename_vector_in)
        if dataset is None:
            raise RuntimeError(f"Could not create output dataset: {sFilename_vector_in}")

        layer = dataset.CreateLayer('vertex', spatial_ref, ogr.wkbPoint)
        if layer is None:
            raise RuntimeError(f"Could not create layer in output dataset: {sFilename_vector_in}")

        # Create fields
        layer.CreateField(ogr.FieldDefn('pointid', ogr.OFTInteger64))
        if iFlag_attribute:
            layer.CreateField(ogr.FieldDefn('connectivity', ogr.OFTInteger64))

        layer_defn = layer.GetLayerDefn()
        feature = ogr.Feature(layer_defn)

        # Add vertices as features
        for vertex_id, vertex in enumerate(aVertex_in):
            point = ogr.Geometry(ogr.wkbPoint)

            if iFlag_projected == 1:
                if not hasattr(vertex, 'dx') or not hasattr(vertex, 'dy'):
                    raise ValueError(f"Vertex {vertex_id} missing dx or dy attributes for projected coordinates.")
                point.AddPoint(vertex.dx, vertex.dy)
            else:
                if not hasattr(vertex, 'dLongitude_degree') or not hasattr(vertex, 'dLatitude_degree'):
                    raise ValueError(f"Vertex {vertex_id} missing coordinate attributes.")
                point.AddPoint(vertex.dLongitude_degree, vertex.dLatitude_degree)

            geometry = ogr.CreateGeometryFromWkb(point.ExportToWkb())
            if geometry is None:
                raise RuntimeError(f"Could not create geometry for vertex {vertex_id}.")

            feature.SetGeometry(geometry)
            feature.SetField("pointid", vertex_id + 1)

            if iFlag_attribute:
                if vertex_id >= len(attribute_data):
                    raise ValueError(f"Attribute data missing for vertex {vertex_id}.")
                feature.SetField("connectivity", int(attribute_data[vertex_id]))

            result = layer.CreateFeature(feature)
            if result != ogr.OGRERR_NONE:
                raise RuntimeError(f"Could not create feature for vertex {vertex_id}.")

        dataset.FlushCache()

    finally:
        feature = None
        layer = None
        dataset = None


def export_vertex_as_polygon_file(
    aVertex_in: List[Any],
    sFilename_out: str,
    pSpatial_reference_in: Optional[osr.SpatialReference] = None,
    overwrite: bool = True
) -> None:
    """Export vertices as a closed polygon to a vector file.

    Parameters
    ----------
    aVertex_in : list
        List of vertex objects with dLongitude_degree and dLatitude_degree attributes.
    sFilename_out : str
        Output filename with extension. Format is determined from extension.
        Supports GeoJSON, Shapefile, GeoPackage, KML, etc.
    pSpatial_reference_in : osr.SpatialReference, optional
        Spatial reference system. If None, defaults to WGS84 (EPSG:4326).
    overwrite : bool, optional
        If True, overwrite existing file. Default is True.

    Raises
    ------
    ValueError
        If aVertex_in has fewer than 3 vertices or if file extension is not supported.
    FileExistsError
        If the output file exists and overwrite is False.
    RuntimeError
        If dataset, layer, or feature creation fails.

    Notes
    -----
    - Automatically closes the polygon by connecting the last vertex to the first
    - Adds an 'area' field containing the calculated polygon area
    - Automatically determines output format from file extension
    """

    if not aVertex_in or len(aVertex_in) < 3:
        raise ValueError("Polygon requires at least 3 vertices.")

    # Handle overwrite
    if os.path.exists(sFilename_out):
        if not overwrite:
            raise FileExistsError(f"Output file already exists: {sFilename_out}")
        os.remove(sFilename_out)

    # Get driver from file extension
    driver = get_vector_driver_from_extension(sFilename_out)

    if pSpatial_reference_in is None:
        spatial_ref = osr.SpatialReference()
        spatial_ref.ImportFromEPSG(4326)
    else:
        spatial_ref = pSpatial_reference_in

    dataset = None
    layer = None
    feature = None

    try:
        # Create dataset and layer
        dataset = driver.CreateDataSource(sFilename_out)
        if dataset is None:
            raise RuntimeError(f"Could not create output dataset: {sFilename_out}")

        layer = dataset.CreateLayer('polygon', spatial_ref, geom_type=ogr.wkbPolygon)
        if layer is None:
            raise RuntimeError(f"Could not create layer in output dataset: {sFilename_out}")

        # Create area field
        area_field = ogr.FieldDefn('area', ogr.OFTReal)
        area_field.SetWidth(20)
        area_field.SetPrecision(2)
        layer.CreateField(area_field)

        # Build polygon ring
        ring = ogr.Geometry(ogr.wkbLinearRing)
        lon_list = []
        lat_list = []

        for vertex_id, vertex in enumerate(aVertex_in):
            if not hasattr(vertex, 'dLongitude_degree') or not hasattr(vertex, 'dLatitude_degree'):
                raise ValueError(f"Vertex {vertex_id} missing coordinate attributes.")

            ring.AddPoint(vertex.dLongitude_degree, vertex.dLatitude_degree)
            lon_list.append(vertex.dLongitude_degree)
            lat_list.append(vertex.dLatitude_degree)

        # Close the ring
        ring.AddPoint(aVertex_in[0].dLongitude_degree, aVertex_in[0].dLatitude_degree)

        # Calculate area
        area = calculate_polygon_area(lon_list, lat_list)

        # Create polygon geometry
        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)

        # Create and populate feature
        layer_defn = layer.GetLayerDefn()
        feature = ogr.Feature(layer_defn)
        feature.SetGeometry(polygon)
        feature.SetField("area", area)

        result = layer.CreateFeature(feature)
        if result != ogr.OGRERR_NONE:
            raise RuntimeError(f"Could not create polygon feature in {sFilename_out}.")

        dataset.FlushCache()

    finally:
        feature = None
        layer = None
        dataset = None




