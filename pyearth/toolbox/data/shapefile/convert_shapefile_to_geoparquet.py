import os
from osgeo import ogr, gdal

def convert_shapefile_to_geoparquet(sFilename_shapefile_in, sFilename_geoparquet_out = None, sLayername_in = 'layer'):
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    #check if the pDriver_gpkg is available
    if pDriver_shapefile is None:
        print('shapefile driver not available.')
        print("GDAL version:", gdal.__version__)
        return

    pDriver_parquet = ogr.GetDriverByName('Parquet')
    pDataset_in = pDriver_parquet.Open(sFilename_shapefile_in, 0)
    if pDataset_in is None:
        print(f"Failed to open file: {sFilename_shapefile_in}")
        return
    pLayer_in = pDataset_in.GetLayer()

    if sFilename_geoparquet_out is None:
        sFilename_geoparquet_out = sFilename_shapefile_in.replace('.geojson', '.gpkg')

    if os.path.exists(sFilename_geoparquet_out):
        os.remove(sFilename_geoparquet_out)

    pDataset_out = pDriver_parquet.CreateDataSource(sFilename_geoparquet_out)
    pLayer_out = pDataset_out.CopyLayer(pLayer_in, sLayername_in)
    pDataset_in = pDataset_out = None
    pLayer_in = pLayer_out = None

    return

