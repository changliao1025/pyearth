import os
from osgeo import ogr, gdal

def convert_geojson_to_geoparquet(input_geojson, output_geoparquet):
    driver = ogr.GetDriverByName('Parquet')
    #check if the driver is available
    if driver is None:
        print('Parquet driver not available.')
        print("GDAL version:", gdal.__version__)
        return 
    
    driver2 = ogr.GetDriverByName('GeoJSON')
    in_dataset = driver2.Open(input_geojson, 0)
    in_layer = in_dataset.GetLayer()

    if os.path.exists(output_geoparquet):
        os.remove(output_geoparquet)

    out_dataset = driver.CreateDataSource(output_geoparquet)
    out_layer = out_dataset.CopyLayer(in_layer, 'flow')
    in_dataset = out_dataset = None

if __name__ == '__main__':
    input_geojson = '/compyfs/liao313/04model/pyhexwatershed/k34/pyhexwatershed20231001003/hexwatershed/00000001/variable_polygon.geojson'
    output_geoparquet = '/compyfs/liao313/04model/pyhexwatershed/k34/pyhexwatershed20231001003/hexwatershed/00000001/variable_polygon.parquet'
    input_geojson = '/compyfs/liao313/04model/pyhexwatershed/k34/pyhexwatershed20231001003/hexwatershed/00000001/flow_direction.geojson'
    output_geoparquet = '/compyfs/liao313/04model/pyhexwatershed/k34/pyhexwatershed20231001003/hexwatershed/00000001/flow_direction.parquet'
    input_geojson = '/compyfs/liao313/04model/pyhexwatershed/amazon/pyhexwatershed20240101012/pyflowline/dggrid.geojson'
    output_geoparquet = '/compyfs/liao313/04model/pyhexwatershed/amazon/pyhexwatershed20240101012/pyflowline/dggrid.parquet'
   
    convert_geojson_to_geoparquet(input_geojson, output_geoparquet)