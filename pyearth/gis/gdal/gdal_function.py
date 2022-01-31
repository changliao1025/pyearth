import os, sys
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal, gdalconst

gdal.UseExceptions()    # Enable exceptions


def reproject_coordinates(x, y, spatial_reference_source, spatial_reference_target=None):
    """ Reproject a pair of x,y coordinates. 

    Args:
        x ([type]): [description]
        y ([type]): [description]
        spatial_reference_source ([type]): [description]
        spatial_reference_target ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """

    if spatial_reference_target is not None:

        pass
    else:
        spatial_reference_target = osr.SpatialReference()
        spatial_reference_target.ImportFromEPSG(4326)
        
        pass

    
    if int(osgeo.__version__[0]) >= 3:
    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
                    
        spatial_reference_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        spatial_reference_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    
    pTransform = osr.CoordinateTransformation( spatial_reference_source, spatial_reference_target)
   
    x_new,y_new, z = pTransform.TransformPoint( x,y)
    
    return x_new,y_new

def reproject_coordinates_batch(x, y, spatial_reference_source, spatial_reference_target=None):
    """ Reproject a list of x, y coordinates.

    Args:
        x (list): list of x coordinates
        y (list): list of y coordinates
        spatial_reference_source ([type]): [description]
        spatial_reference_target ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """

    if spatial_reference_target is not None:

        pass
    else:
        spatial_reference_target = osr.SpatialReference()
        spatial_reference_target.ImportFromEPSG(4326)
        
        pass

    
    if int(osgeo.__version__[0]) >= 3:
    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
                    
        spatial_reference_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        spatial_reference_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    
    pTransform = osr.CoordinateTransformation( spatial_reference_source, spatial_reference_target)

    npoint = len(x)
    x_new=list()
    y_new=list()
    for i in range(npoint):
        x0 = x[i]
        y0 = y[i]
   
        x1,y1, z = pTransform.TransformPoint( x0,y0)

        x_new.append(x1)
        y_new.append(y1)
    
    return x_new,y_new

def obtain_raster_metadata(sFilename_geotiff):
    """retrieve the metadata of a geotiff file

    Args:
        sFilename_geotiff ([type]): [description]

    Returns:
        [type]: [description]
    """
    
    #pDriver = gdal.GetDriverByName('GTiff')
   
    pDataset = gdal.Open(sFilename_geotiff, gdal.GA_ReadOnly)

    if pDataset is None:
        print("Couldn't open this file: " + sFilename_geotiff)
        sys.exit("Try again!")
    else: 
        pProjection = pDataset.GetProjection()
        pSpatialRef = osr.SpatialReference(wkt=pProjection)
    
    
        ncolumn = pDataset.RasterXSize
        nrow = pDataset.RasterYSize
        #nband = pDataset.RasterCount

        pGeotransform = pDataset.GetGeoTransform()
        dOriginX = pGeotransform[0]
        dOriginY = pGeotransform[3]
        dPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]       
        
        print( dPixelWidth, dOriginX, dOriginY, nrow, ncolumn)
        return dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, pSpatialRef, pProjection, pGeotransform

def obtain_shapefile_metadata(sFilename_shapefile):
    """[summary]

    Args:
        sFilename_shapefile ([type]): [description]

    Returns:
        [type]: [description]
    """
    if os.path.exists(sFilename_shapefile):
        pass
    else:
        print('The  shapefile does not exist!')
        return

    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
   
    pDataset = pDriver_shapefile.Open(sFilename_shapefile, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    #pSrs = pLayer.GetSpatialRef()

    if pDataset is None:
        print("Couldn't open this file: " + sFilename_shapefile)
        sys.exit("Try again!")
    else:    
        
        iFlag_first=1
    
        for feature in pLayer:
            pGeometry = feature.GetGeometryRef()
            pEnvelope = pGeometry.GetEnvelope()

            #"minX: %d, minY: %d, maxX: %d, maxY: %d" %(env[0],env[2],env[1],env[3])

            if iFlag_first ==1:
                left_min = pEnvelope[0]
                right_max =  pEnvelope[1]
                bot_min =  pEnvelope[2]
                top_max =  pEnvelope[3]
                iFlag_first = 0
            else:
                left_min = np.min([left_min,  pEnvelope[0]])
                right_max = np.max([right_max,  pEnvelope[1]])
                bot_min = np.min([bot_min,  pEnvelope[2]])
                top_max = np.max([top_max,  pEnvelope[3]])

        

            print( pEnvelope )
      
        return left_min, right_max, bot_min, top_max

