import os
from osgeo import osr, ogr
#accepte user input arguments

def filter_vector_by_attribute(sFilename_input, sFilename_output, sAttribute_name, dValue_filter):
    # Open the input shapefile
    pDriver_parquet = ogr.GetDriverByName('Parquet')
    pDataSource_in = pDriver_parquet.Open(sFilename_input, 0)
    if pDataSource_in is None:
        print("Could not open input file")
        return    
    
    if os.path.exists(sFilename_output):
        os.remove(sFilename_output)
    
    pLayer_in = pDataSource_in.GetLayer()
    #get spatial reference
    pSpatialRef = pLayer_in.GetSpatialRef()
    # Create a new shapefile for the filtered data with the same spatial reference

    pDataSource_out = pDriver_parquet.CreateDataSource(sFilename_output)
    pLayer_out = pDataSource_out.CreateLayer('filtered', geom_type=ogr.wkbPolygon, srs=pSpatialRef)

    # Copy the fields from the input layer to the output layer
    pLayer_defn_in = pLayer_in.GetLayerDefn()
    for i in range(pLayer_defn_in.GetFieldCount()):
        pField_defn = pLayer_defn_in.GetFieldDefn(i)
        pLayer_out.CreateField(pField_defn)

    # Create a SQL query to filter by the attribute value
    sql_query = f"{sAttribute_name} = '{dValue_filter}'"

    # Apply the filter and copy features to the output layer
    pLayer_in.SetAttributeFilter(sql_query)
    for feature in pLayer_in:
        pLayer_out.CreateFeature(feature)

    # Close the data sources
    pDataSource_in = pDataSource_out = None