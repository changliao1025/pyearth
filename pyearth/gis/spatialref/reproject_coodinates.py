
import osgeo
from osgeo import  osr

#both function use wkt instead of srs object

def reproject_coordinates(dx_in, dy_in, pProjection_source_in, pProjection_target_in=None):
    """ Reproject a pair of x,y coordinates.

    Args:
        dx_in (float): X Coordinate of point
        dy_in (float): Y Coordinate of
        pSpatial_reference_source_in (osr): The source spatial reference of point
        pSpatial_reference_target_in (osr, optional): The target spatial reference of point. Defaults to None.

    Returns:
        Tuple: dx_new, dy_new
    """

    if pProjection_target_in is not None:
        #convert the projection to spatial reference
        pSpatial_reference_target = osr.SpatialReference()
        pSpatial_reference_target.ImportFromWkt(pProjection_target_in)
    else:
        pSpatial_reference_target = osr.SpatialReference()
        pSpatial_reference_target.ImportFromEPSG(4326)
        pass

    pSpatial_reference_source = osr.SpatialReference()
    pSpatial_reference_source.ImportFromWkt(pProjection_source_in)
    if int(osgeo.__version__[0]) >= 3:
    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546

        pSpatial_reference_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        pSpatial_reference_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)


    pTransform = osr.CoordinateTransformation( pSpatial_reference_source, pSpatial_reference_target)

    dx_out, dy_out, z = pTransform.TransformPoint( dx_in,dy_in)

    return dx_out, dy_out

def reproject_coordinates_batch(aX_in, aY_in, pProjection_source_in, pProjection_target_in=None):
    """ Reproject a list of x, y coordinates.

    Args:
        aX_in (list): A list of X Coordinate of points
        aY_in (list): A list of Y Coordinate of points
        pSpatial_reference_source_in (osr): The source spatial reference of point
        pSpatial_reference_target_in (osr, optional): The target spatial reference of point. Defaults to None.

    Returns:
        Tuple: aX_out, aY_out
    """

    if pProjection_target_in is not None:
        pSpatial_reference_target = osr.SpatialReference()
        pSpatial_reference_target.ImportFromWkt(pProjection_target_in)
    else:
        pSpatial_reference_target = osr.SpatialReference()
        pSpatial_reference_target.ImportFromEPSG(4326)

    pSpatial_reference_source = osr.SpatialReference()
    pSpatial_reference_source.ImportFromWkt(pProjection_source_in)
    if int(osgeo.__version__[0]) >= 3:
    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546

        pSpatial_reference_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        pSpatial_reference_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)


    pTransform = osr.CoordinateTransformation( pSpatial_reference_source, pSpatial_reference_target)

    npoint = len(aX_in)
    aX_out=list()
    aY_out=list()
    for i in range(npoint):
        x0 = aX_in[i]
        y0 = aY_in[i]

        x1,y1, z = pTransform.TransformPoint( x0,y0)

        aX_out.append(x1)
        aY_out.append(y1)

    return aX_out,aY_out

