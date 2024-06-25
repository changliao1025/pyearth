import numpy as np
from osgeo import ogr
def retrieve_field_value(sFilename_vector_in, sField_name_in, dMissing_value = None):
    # Open the vector file
    ds = ogr.Open(sFilename_vector_in)
    if ds is None:
        raise ValueError(f"Could not open {sFilename_vector_in}")

    # Get the first layer
    layer = ds.GetLayer(0)

    # Collect all values of the field
    if dMissing_value is None:
        aValues = [feature.GetField(sField_name_in) for feature in layer]
    else:
        aValues = []
        for feature in layer:
            dValue = feature.GetField(sField_name_in)
            if dValue != dMissing_value:
               aValues.append(dValue)


    # convert to numpy
    aValues = np.array(aValues)
    return aValues


