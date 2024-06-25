import os
from osgeo import ogr, gdal

def convert_geojson_to_geoparquet(sFilename_geojson_in, sFilename_geoparqiet_out = None, sLayername_in = 'layer'):
    pDrive_parquet = ogr.GetDriverByName('Parquet')
    #check if the pDrive_parquet is available
    if pDrive_parquet is None:
        print('Parquet pDrive_parquet not available.')
        print("GDAL version:", gdal.__version__)
        return

    pDrive_geojson = ogr.GetDriverByName('GeoJSON')
    pDateset_in = pDrive_geojson.Open(sFilename_geojson_in, 0)
    if pDateset_in is None:
        print(f"Failed to open file: {sFilename_geojson_in}")
        return

    pLayer_in = pDateset_in.GetLayer()

    if sFilename_geoparqiet_out is None:
        sFilename_geoparqiet_out = sFilename_geojson_in.replace('.geojson', '.parquet')

    if os.path.exists(sFilename_geoparqiet_out):
        os.remove(sFilename_geoparqiet_out)

    pDataset_out = pDrive_parquet.CreateDataSource(sFilename_geoparqiet_out)
    pLayer_out = pDataset_out.CopyLayer(pLayer_in, sLayername_in)
    pDateset_in = pDataset_out = None
    pLayer_in = pLayer_out = None

    return