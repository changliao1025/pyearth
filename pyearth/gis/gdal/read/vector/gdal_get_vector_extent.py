from osgeo import gdal, osr, ogr
def gdal_get_vector_extent(sFilename_in, iFlag_return_union_geometry=False):
    """
    Get the extent of a vector file.

    This function can also return a single geometry that is a collection of all
    feature geometries in the source file.

    Args:
        sFilename_in (str): The input vector file.
        return_union_geometry (bool, optional): If True, return a collection of
            all geometries. Defaults to False.

    Returns:
        tuple or None:
            If return_union_geometry is False, returns the extent as a
            tuple (min_x, max_x, min_y, max_y), or None on error.
            If return_union_geometry is True, returns a tuple (extent, union_geometry),
            or (None, None) on error. The union_geometry will be a single OGR
            Geometry object (e.g., MultiPolygon) containing all geometries from
            the input layer.
    """
    pDataset = None
    try:
        pDataset = ogr.Open(sFilename_in)
        if pDataset is None:
            print(f"Could not open {sFilename_in}")
            return (None, None) if iFlag_return_union_geometry else None

        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            print(f"Could not get layer from {sFilename_in}")
            return (None, None) if iFlag_return_union_geometry else None

        aExtent = pLayer.GetExtent()

        if not iFlag_return_union_geometry:
            return aExtent

        # The following part creates a collection of geometries.
        pLayer_geom_type = pLayer.GetGeomType()

        # Determine the target collection type using the flattened geometry type
        #base_geom_type = ogr.wkbFlatten(pLayer_geom_type)
        base_geom_type = pLayer_geom_type & 0xff
        if base_geom_type == ogr.wkbPoint:
            union_geom = ogr.Geometry(ogr.wkbMultiPoint)
        elif base_geom_type == ogr.wkbLineString:
            union_geom = ogr.Geometry(ogr.wkbMultiLineString)
        elif base_geom_type == ogr.wkbPolygon:
            union_geom = ogr.Geometry(ogr.wkbMultiPolygon)
        else:
            sGeomType = ogr.GeometryTypeToName(pLayer_geom_type)
            print(f"Geometry type '{sGeomType}' not supported for geometry collection.")
            return aExtent, None

        pLayer.ResetReading()
        for pFeature in pLayer:
            geometry = pFeature.GetGeometryRef()
            if geometry is not None:
                geom_clone = geometry.Clone()
                geom_clone.FlattenTo2D()

                flat_geom_type = geom_clone.GetGeometryType()

                if flat_geom_type in [ogr.wkbMultiPoint, ogr.wkbMultiLineString, ogr.wkbMultiPolygon]:
                    for i in range(geom_clone.GetGeometryCount()):
                        part = geom_clone.GetGeometryRef(i)
                        union_geom.AddGeometry(part)
                elif flat_geom_type in [ogr.wkbPoint, ogr.wkbLineString, ogr.wkbPolygon]:
                    union_geom.AddGeometry(geom_clone)

        return aExtent, union_geom

    except Exception as e:
        print(f"Error in gdal_get_vector_extent: {e}")
        return (None, None) if iFlag_return_union_geometry else None
    finally:
        if pDataset is not None:
            pDataset = None