import os
from osgeo import ogr

def convert_geojson_to_geoparquet(input_geojson, output_geoparquet):
    driver = ogr.GetDriverByName('Parquet')
    driver2 = ogr.GetDriverByName('GeoJSON')
    in_dataset = driver2.Open(input_geojson, 0)
    in_layer = in_dataset.GetLayer()

    if os.path.exists(output_geoparquet):
        os.remove(output_geoparquet)

    out_dataset = driver.CreateDataSource(output_geoparquet)
    out_layer = out_dataset.CopyLayer(in_layer, 'mesh')
    in_dataset = out_dataset = None