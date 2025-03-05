import os
import numpy as np
import osgeo
from osgeo import gdal, osr
from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.gdal.gdal_to_numpy_datatype import gdal_to_numpy_datatype
def reproject_raster(sFilename_raster_in, sFilename_raster_out, pProjection_target, iFlag_force_resolution_in = 0):
    # Open the input raster file
    pDriver_tiff = gdal.GetDriverByName('GTiff')
    pDataset_raster = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)
    if pDataset_raster is None:
        raise FileNotFoundError(f"Could not open {sFilename_raster_in}")

    # Get the source projection
    pProjection_source = pDataset_raster.GetProjection()
    pSpatial_reference_source = osr.SpatialReference()
    pSpatial_reference_source.ImportFromWkt(pProjection_source)

    # Create the target projection
    pSpatial_reference_target = osr.SpatialReference()
    pSpatial_reference_target.ImportFromWkt(pProjection_target)

    # Get the geotransform and dimensions of the input raster
    geotransform = pDataset_raster.GetGeoTransform()
    width = pDataset_raster.RasterXSize
    height = pDataset_raster.RasterYSize

    # Create the output raster file
    driver = gdal.GetDriverByName('GTiff')
    pDataset_transform = driver.Create(sFilename_raster_out, width, height, pDataset_raster.RasterCount, gdal.GDT_Float32)
    pDataset_transform.SetGeoTransform(geotransform)
    pDataset_transform.SetProjection(pSpatial_reference_target.ExportToWkt())

    # Reproject each band
    for i in range(1, pDataset_raster.RasterCount + 1):
        band = pDataset_raster.GetRasterBand(i)
        data = band.ReadAsArray()
        gdal.ReprojectImage(pDataset_raster, pDataset_transform, pProjection_source, pProjection_target, gdal.GRA_Bilinear)

    # Flush the cache and close the datasets
    pDataset_transform.FlushCache()
    pDataset_raster = None
    pDataset_transform = None

    return

def reproject_raster_gdalwarp(sFilename_raster_in, sFilename_raster_out, pProjection_target,
                              xRes=None, yRes=None,
                               sResampleAlg = 'near',
                                 iFlag_force_resolution_in = 0):
    pDriver_tiff = gdal.GetDriverByName('GTiff')
    # Open the input raster file
    pSpatialRef_default = osr.SpatialReference()
    pSpatialRef_default.ImportFromEPSG(4326)  # WGS84
    sProjection_default = pSpatialRef_default.ExportToWkt()

    pDataset_raster = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)
    if pDataset_raster is None:
        raise FileNotFoundError(f"Could not open {sFilename_raster_in}")

    #get the nodata value
    band = pDataset_raster.GetRasterBand(1)
    dMissing_value_source = band.GetNoDataValue()
    dMissing_value_target = dMissing_value_source
    #get the min and max value
    dMin = band.GetMinimum()
    dMax = band.GetMaximum()

    #get the datatype
    iData_type = band.DataType
    #get projection
    pProjection_source = pDataset_raster.GetProjection()
    pSpatialRef_source = osr.SpatialReference()
    pSpatialRef_source.ImportFromWkt(pProjection_source)
    if pProjection_source == sProjection_default:
        iFlag_wgs84_source = 1
    else:
        iFlag_wgs84_source = 0

    if pProjection_target == sProjection_default:
        iFlag_wgs84_target = 1
    else:
        iFlag_wgs84_target = 0

    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(pProjection_target)
    if int(osgeo.__version__[0]) >= 3:
        # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
        pSpatialRef_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        pSpatialRef_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    if xRes is None or yRes is None:
        xRes = pDataset_raster.GetGeoTransform()[1]
        yRes = pDataset_raster.GetGeoTransform()[5]
    else:
        dResolution_x = xRes
        dResolution_y = abs(yRes)
        iFlag_force_resolution_in = 1

    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']
    if iFlag_force_resolution_in == 1:
        pWrapOption = gdal.WarpOptions( cropToCutline=False, #could be true if vector file is provided
                                    xRes=dResolution_x,
                                   yRes=dResolution_y,
                                        dstSRS=pSpatialRef_target , format = 'MEM',
                                        resampleAlg=sResampleAlg ) #this resample algorithm may be provided as an input argument
        pDataset_clip_warped = gdal.Warp('', pDataset_raster, options=pWrapOption)#gdal.Warp(sFilename_out, pDataset_in, options=pWrapOption)
        newGeoTransform = pDataset_clip_warped.GetGeoTransform()
        aData_clip = pDataset_clip_warped.ReadAsArray()
        iNewWidth = aData_clip.shape[1]
        iNewHeigh = aData_clip.shape[0]
        pDataset_clip = pDriver_tiff.Create(sFilename_raster_out, iNewWidth, iNewHeigh, 1, iData_type, options= options) #this data type may be provided as an input argument
        pDataset_clip.GetRasterBand(1).SetNoDataValue(dMissing_value_target)
        pDataset_clip.SetGeoTransform( newGeoTransform )
        pDataset_clip.SetProjection( pProjection_target)
        pass
    else:
        if iFlag_wgs84_source == 1:       #from wgs84
            dLon_min, dLon_max, dLat_min, dLat_max = gdal_get_raster_extent(sFilename_raster_in)
            if iFlag_wgs84_target == 1: #to wgs84
                #if the target is wgs84, we can use the same approach
                nleft  = np.floor(  (dLon_min - (-180)) /(dResolution_x)  )
                nright = np.ceil(  (dLon_max - (-180)) /(dResolution_x)  )
                ntop  = np.floor(  (90 - dLat_max) /(dResolution_y)  )
                nbot = np.ceil(  (90 - dLat_min) /(dResolution_y)  )
                nrow = max(1, int(nbot - ntop))
                ncolumn = max(1, int(nright - nleft))
                dMin_x = -180 + nleft*dResolution_x
                dMax_y = 90 - ntop*dResolution_y
            else: #to something else
                #if the target is not wgs84, we need to reproject the extent
                if int(osgeo.__version__[0]) >= 3:
                    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
                    pSpatialRef_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                    pSpatialRef_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                pTransform = osr.CoordinateTransformation(pSpatialRef_source, pSpatialRef_target)
                dX_min, dY_min, _ = pTransform.TransformPoint(dLon_min, dLat_min)
                dX_max, dY_max, _ = pTransform.TransformPoint(dLon_max, dLat_max)
                nrow = max(1, int((dY_max - dY_min) / dResolution_y))
                ncolumn = max(1, int((dX_max - dX_min) / dResolution_x))
                dMin_x = dX_min
                dMax_y = dY_max
        else:
            #if it is not wgs84, it has a projection already
            #get the extent of raster geotiff file
            dX_min, dX_max, dY_min, dY_max = gdal_get_raster_extent(sFilename_raster_in)
            #if the projection is not the same, we need to reproject the extent
            if iFlag_wgs84_target == 1: #to wgs84
                if int(osgeo.__version__[0]) >= 3:
                    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
                    pSpatialRef_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                    pSpatialRef_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                pTransform = osr.CoordinateTransformation(pSpatialRef_source, pSpatialRef_target)
                dLon_min, dLat_min, _ = pTransform.TransformPoint(dX_min, dY_min)
                dLon_max, dLat_max, _ = pTransform.TransformPoint(dX_max, dY_max)
                nleft  = np.floor(  (dLon_min - (-180)) /(dResolution_x)  )
                nright = np.ceil(  (dLon_max - (-180)) /(dResolution_x)  )
                ntop  = np.floor(  (90 - dLat_max) /(dResolution_y)  )
                nbot = np.ceil(  (90 - dLat_min) /(dResolution_y)  )
                nrow = int(nbot-ntop)
                ncolumn = int(nright - nleft)
                dMin_x = -180 + nleft*dResolution_x
                dMax_y = 90 - ntop*dResolution_y
            else: #to something else, but could be the same or different
                if pProjection_source == pProjection_target:
                    nrow = int((dY_max - dY_min)/dResolution_y)
                    ncolumn = int((dX_max - dX_min)/dResolution_x)
                    dMin_x = dX_min
                    dMax_y = dY_max
                else: #not the same
                    if int(osgeo.__version__[0]) >= 3:
                        # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
                        pSpatialRef_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                        pSpatialRef_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                    pTransform = osr.CoordinateTransformation(pSpatialRef_source, pSpatialRef_target)
                    dX_min, dY_min, _ = pTransform.TransformPoint(dX_min, dY_min)
                    dX_max, dY_max, _ = pTransform.TransformPoint(dX_max, dY_max)
                    nrow = int((dY_max - dY_min)/dResolution_y)
                    ncolumn = int((dX_max - dX_min)/dResolution_x)
                    dMin_x = dX_min
                    dMax_y = dY_max

        iNewWidth = ncolumn
        iNewHeigh = nrow
        newGeoTransform = (dMin_x, dResolution_x, 0, dMax_y, 0, -dResolution_y)
        pDataset_clip = pDriver_tiff.Create(sFilename_raster_out, iNewWidth, iNewHeigh, 1, iData_type, options= options) #this data type may be provided as an input argument
        pDataset_clip.GetRasterBand(1).SetNoDataValue(dMissing_value_target)
        pDataset_clip.SetGeoTransform( newGeoTransform )
        pDataset_clip.SetProjection( pProjection_target)
        pWrapOption = gdal.WarpOptions( cropToCutline=False, #could be true if vector file is provided
                                width=iNewWidth,
                                    height=iNewHeigh,
                                        dstSRS=pSpatialRef_target , format = 'MEM',
                                        resampleAlg=sResampleAlg ) #this resample algorithm may be provided as an input argument

        pDataset_clip_warped = gdal.Warp('', pDataset_raster, options=pWrapOption)#gdal.Warp(sFilename_out, pDataset_in, options=pWrapOption)
        #convert the warped dataset to an array
        aData_clip = pDataset_clip_warped.ReadAsArray()

    #change the gdal data type to numpy data type
    iData_type_numpy = gdal_to_numpy_datatype(iData_type)
    aData_clip = aData_clip.astype(iData_type_numpy)
    aData_clip[aData_clip == dMissing_value_source] = dMissing_value_target

    dummy_index = np.where(aData_clip < dMin)
    aData_clip[dummy_index] = dMissing_value_target
    dummy_index = np.where(aData_clip > dMax)
    aData_clip[dummy_index] = dMissing_value_target
    pDataset_clip.GetRasterBand(1).WriteArray(aData_clip)
    pDataset_clip.GetRasterBand(1).FlushCache()  # Corrected method name to FlushCache()
    pDataset_clip.FlushCache()

    # Close the input dataset
    pDataset_raster = None
    pSpatialRef_source = None
    pSpatialRef_target = None
    pSpatialRef_default = None
    pSpatial_reference_source = None

    return