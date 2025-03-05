import os
import numpy as np
from osgeo import osr, gdal, ogr

def merge_vector_polygon_files(iFiletype_in,
                              aFilename_in,
                              sFilename_out,
                              aVariable_in=None,
                              sVariable_out=None,
                              dMissing_value_in=None,
                              dData_max_in=None,
                              dData_min_in=None,
                              ):
    """
    Map a vector polyline data

    Args:

    """
    if aVariable_in is not None:

        sVar = sVariable_out[0:4].lower()



    if iFiletype_in == 1:  # geojson
        pDriver = ogr.GetDriverByName('GeoJSON')
    else:
        if iFiletype_in == 2:  # shapefile
            pDriver = ogr.GetDriverByName('Esri Shapefile')

    nDataset = len(aFilename_in)
    aDataset = []
    for sFilename_in in aFilename_in:
        pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
        aDataset.append(pDataset)

    pLayer = aDataset[0].GetLayer(0)

    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon

    # get the number of features
    nFeature = pLayer.GetFeatureCount()

    if os.path.exists(sFilename_out):
        os.remove(sFilename_out)

    pDataset_out = pDriver.CreateDataSource(sFilename_out)
    pLayer_out = pDataset_out.CreateLayer(
        'cell', pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    # long type for high resolution
    pLayer_out.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64))
    # long type for high resolution
    if aVariable_in is not None:
        pLayer_out.CreateField(ogr.FieldDefn(sVar, ogr.OFTReal))
    pLayerDefn = pLayer_out.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayerDefn)

    lID = 1
    for iFeature in range(nFeature):
        if aVariable_in is not None:
            aValue = list()
        for iData in range(nDataset):
            pLayer = aDataset[iData].GetLayer(0)
            pFeature = pLayer.GetFeature(iFeature)
            if iData == 0:  # use the first data as geometry
                pGeometry_in = pFeature.GetGeometryRef()
                pGeometry_out = pGeometry_in.Clone()

            if aVariable_in is not None:
                sVar_dummy = aVariable_in[iData].lower()
                sVar_dummy = sVar_dummy[0:4]
                dValue = float(pFeature.GetField(sVar_dummy))
                aValue.append(dValue)

        if aVariable_in is not None:
            aValue = np.array(aValue)
            # take the sum
            dValue = np.sum(aValue)

        pFeature_out.SetGeometry(pGeometry_out)
        pFeature_out.SetField('id', lID)
        if aVariable_in is not None:
            pFeature_out.SetField(sVar, float(dValue))
            
        pLayer_out.CreateFeature(pFeature_out)
        lID = lID + 1

    pDataset = pLayer = pFeature = None  # save, close

    print("finished")
