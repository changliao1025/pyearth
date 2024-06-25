import os
from osgeo import ogr, gdal

def convert_geojson_to_geopackage(sFilename_geojson_in, sFilename_geopackage_out = None, sLayername_in = 'layer'):
    pDriver_gpkg = ogr.GetDriverByName('GPKG')
    #check if the pDriver_gpkg is available
    if pDriver_gpkg is None:
        print('GPKG pDriver_gpkg not available.')
        print("GDAL version:", gdal.__version__)
        return

    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    pDataset_in = pDriver_geojson.Open(sFilename_geojson_in, 0)
    if pDataset_in is None:
        print(f"Failed to open file: {sFilename_geojson_in}")
        return
    pLayer_in = pDataset_in.GetLayer()

    if sFilename_geopackage_out is None:
        sFilename_geopackage_out = sFilename_geojson_in.replace('.geojson', '.gpkg')

    if os.path.exists(sFilename_geopackage_out):
        os.remove(sFilename_geopackage_out)

    pDataset_out = pDriver_gpkg.CreateDataSource(sFilename_geopackage_out)
    pLayer_out = pDataset_out.CopyLayer(pLayer_in, sLayername_in)
    pDataset_in = pDataset_out = None
    pLayer_in = pLayer_out = None

if __name__ == '__main__':
    sFilename_geojson_in = '/compyfs/liao313/04model/pyhexwatershed/k34/pyhexwatershed20231001003/hexwatershed/00000001/variable_polygon.geojson'
    sFilename_geopackage_out = '/compyfs/liao313/04model/pyhexwatershed/k34/pyhexwatershed20231001003/hexwatershed/00000001/variable_polygon.parquet'
    sFilename_geojson_in = '/compyfs/liao313/04model/pyhexwatershed/k34/pyhexwatershed20231001003/hexwatershed/00000001/flow_direction.geojson'
    sFilename_geopackage_out = '/compyfs/liao313/04model/pyhexwatershed/k34/pyhexwatershed20231001003/hexwatershed/00000001/flow_direction.parquet'
    sFilename_geojson_in = '/compyfs/liao313/04model/pyhexwatershed/amazon/pyhexwatershed20240101012/pyflowline/dggrid.geojson'
    sFilename_geopackage_out = '/compyfs/liao313/04model/pyhexwatershed/amazon/pyhexwatershed20240101012/pyflowline/dggrid.gpkg'

    convert_geojson_to_geopackage(sFilename_geojson_in, sFilename_geopackage_out=sFilename_geopackage_out)