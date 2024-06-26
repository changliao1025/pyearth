import os
from osgeo import gdal


def gdal_write_geotiff_file(sFilename_in,
                            aData_in,
                            dPixelWidth_in,
                            dPixelHeight_in,
                            dOriginX_in,
                            dOriginY_in,
                            dMissing_value_in,
                            pProjection_in,
                            datatype=gdal.GDT_Float32):
    """
    Write a Geotiff standard format raster file

    Args:
        sFilename_in (string): The filename
        aData_in (numpy.array): The data
        dPixelWidth_in (float): The resolution
        dOriginX_in (float): The location of origin x
        dOriginY_in (float): The location of origin y
        dMissing_value_in (float): The missinge value
        pProjection_in (osr): The spatial reference

    Returns:
        Tuple: pGeotransform_out, pProjection_out
    """
    if os.path.exists(sFilename_in):
        os.remove(sFilename_in)
        pass

    sDriverName = 'GTiff'
    pDriver = gdal.GetDriverByName(sDriverName)

    if pDriver is None:
        print("%s pDriver not available.\n" % sDriverName)
    else:
        print("%s pDriver IS available.\n" % sDriverName)

    nrow, ncolumn = aData_in.shape
    nband = 1

    pDataset = pDriver.Create(
        sFilename_in,
        ncolumn,
        nrow,
        nband,
        datatype)

    pDataset.SetGeoTransform([
        dOriginX_in,    # 0
        dPixelWidth_in,  # 1
        0,                      # 2
        dOriginY_in,    # 3
        0,                      # 4
        dPixelHeight_in])

    pDataset.SetProjection(pProjection_in)

    pBand = pDataset.GetRasterBand(1)
    pBand.WriteArray(aData_in)
    pBand.SetNoDataValue(dMissing_value_in)

    pGeotransform_out = pDataset.GetGeoTransform()
    pProjection_out = pDataset.GetProjection()

    pDataset.FlushCache()  # Write to disk.
    pDriver = None
    pDataset = None
    pBand = None
    return sFilename_in, pGeotransform_out, pProjection_out


def gdal_write_geotiff_file_multiple_band(sFilename_in,
                                          aData_in,
                                          dPixelWidth_in,
                                          dPixelHeight_in,
                                          dOriginX_in,
                                          dOriginY_in,
                                          dMissing_value_in,
                                          pProjection_in,
                                          datatype=gdal.GDT_Float32):
    """
    Write a multi-band geotiff raster file

    Args:
        sFilename_in (string): The filename
        aData_in (numpy.array): The data
        dPixelWidth_in (float): The resolution
        dOriginX_in (float): The location of origin x
        dOriginY_in (float): The location of origin y
        dMissing_value_in (float): The missinge value
        pProjection_in (osr): The spatial reference

    Returns:
        Tuple: pGeotransform_out, pProjection_out
    """
    if os.path.exists(sFilename_in):
        os.remove(sFilename_in)
        pass

    sDriverName = 'GTiff'
    pDriver = gdal.GetDriverByName(sDriverName)

    if pDriver is None:
        print("%s pDriver not available.\n" % sDriverName)
    else:
        print("%s pDriver IS available.\n" % sDriverName)

    nband, nrow, ncolumn = aData_in.shape

    # Creates a new raster data source
    pDataset = pDriver.Create(sFilename_in,
                              ncolumn,
                              nrow,
                              nband,
                              datatype)

    # Write metadata

    pDataset.SetGeoTransform([dOriginX_in,
                              dPixelWidth_in,
                              0.0,
                              dOriginY_in,
                              0.0,
                              dPixelHeight_in])


    pDataset.SetProjection(pProjection_in)

    # Write raster datasets
    for iBand in range(nband):
        pBand = pDataset.GetRasterBand(iBand + 1)
        pBand.WriteArray(aData_in[iBand, :, :])
        pBand.SetNoDataValue(dMissing_value_in)

    pGeotransform_out = pDataset.GetGeoTransform()
    pProjection_out = pDataset.GetProjection()

    pDataset.FlushCache()  # Write to disk.
    pDriver = None
    pDataset = None
    pBand = None

    return sFilename_in, pGeotransform_out, pProjection_out
