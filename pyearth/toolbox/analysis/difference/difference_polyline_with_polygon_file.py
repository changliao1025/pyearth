#import gdal
import os
from osgeo import ogr, osr
from datetime import datetime
import importlib

def calculate_polyline_polygon_difference(sFilename_base,
                                              sFilename_new,
                                               sFilename_difference_out  ):

    iFlag_cython = importlib.util.find_spec("cython")
    iFlag_rtree = importlib.util.find_spec("tinyr")
    if iFlag_cython is not None:
        from tinyr import RTree
    else:
        import rtree #regular rtree?


    start_time = datetime.now()

    #read the base file
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    #delete if file exist
    if os.path.exists(sFilename_difference_out):
        os.remove(sFilename_difference_out)

    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_difference_out)
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon

    pLayerOut = pDataset_out.CreateLayer('diff', pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayerOut.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution

    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)
    #read the base file

    pDataset_base = pDriver_geojson.Open(sFilename_base, 0)
    #find the numebr of geometries in the base file
    pLayer_base = pDataset_base.GetLayer()
    nFeature_base = pLayer_base.GetFeatureCount()
    #check whether the base file has the attribute ID
    pLayerDefn_base = pLayer_base.GetLayerDefn()

    interleaved = True
    index_base = RTree(interleaved=interleaved, max_cap=5, min_cap=2)

    #for i in range(nFeature_base):
    for idx, pFeature_base in enumerate(pLayer_base):
        ID = idx + 1
        fid = pFeature_base.GetFID()
        pGeometry_base = pFeature_base.GetGeometryRef()
        left, right, bottom, top= pGeometry_base.GetEnvelope()
        pBound= (left, bottom, right, top)
        index_base.insert(fid, pBound)  #

    #read the new file
    pDataset_new = pDriver_geojson.Open(sFilename_new, 0)
    pLayer_new = pDataset_new.GetLayer()
    nFeature_new =  pLayer_new.GetFeatureCount()
    lID_polygon = 1

    for idx, pFeature_new in enumerate(pLayer_new):
        ID = idx + 1
        pGeometry_new = pFeature_new.GetGeometryRef()
        left, right, bottom, top= pGeometry_new.GetEnvelope()
        pBound= (left, bottom, right, top)
        aIntersect = list(index_base.search(pBound))
        for k in aIntersect:
            ID2 = k + 1
            pFeature_base = pLayer_base.GetFeature(k)
            pGeometry_base = pFeature_base.GetGeometryRef()
            iFlag_intersect = pGeometry_new.Intersects( pGeometry_base )
            if( iFlag_intersect == True):
                pGeometry_diff = pGeometry_base.Difference(pGeometry_new)
                pGeometrytype_intersect = pGeometry_diff.GetGeometryName()
                if pGeometrytype_intersect == 'POLYGON':
                    pFeatureOut.SetGeometry(pGeometry_diff)
                    pFeatureOut.SetField("id", lID_polygon)
                    pLayerOut.CreateFeature(pFeatureOut)
                    lID_polygon = lID_polygon + 1
                else:
                    for i in range(pGeometry_diff.GetGeometryCount()):
                        pGeometry_intersect2 = pGeometry_diff.GetGeometryRef(i)
                        pFeatureOut.SetGeometry(pGeometry_intersect2)
                        pFeatureOut.SetField("id", lID_polygon)
                        pLayerOut.CreateFeature(pFeatureOut)
                        lID_polygon = lID_polygon + 1

    #close files
    pDataset_base = None
    pDataset_new = None
    end_time = datetime.now()
    # get difference
    delta = end_time - start_time

    sec = delta.total_seconds()
    print('difference in seconds:', sec)
    min = sec / 60
    print('difference in minutes:', min)

    return None

if __name__ == '__main__':
    pass