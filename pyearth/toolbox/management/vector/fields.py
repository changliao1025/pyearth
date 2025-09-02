import numpy as np
from osgeo import ogr
def get_field_value(sFilename_vector_in, sField_name_in, dMissing_value = None):
    # Open the vector file
    pDataset = ogr.Open(sFilename_vector_in)
    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_vector_in}")

    # Get the first pLayer
    pLayer = pDataset.GetLayer(0)

    # Collect all values of the field
    if dMissing_value is None:
        aValues = [pFeature.GetField(sField_name_in) for pFeature in pLayer]
    else:
        aValues = []
        for pFeature in pLayer:
            dValue = pFeature.GetField(sField_name_in)
            if dValue != dMissing_value:
               aValues.append(dValue)


    # convert to numpy
    aValues = np.array(aValues)
    return aValues

def get_field_and_value(sFilename_vector_in):
    # Open the vector file
    pDataset = ogr.Open(sFilename_vector_in)
    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_vector_in}")

    # Get the first pLayer
    pLayer = pDataset.GetLayer(0)

    # Collect all field names and values
    aField_names = [field.name for field in pLayer.schema]
    aValues= list()


    for pFeature in pLayer:
        feature_dict = {}
        # Collect all field values for this feature
        for sFieldname in aField_names:
            feature_dict[sFieldname] = pFeature.GetField(sFieldname)

        # Add this feature's dictionary to the list
        aValues.append(feature_dict)

    return aField_names, aValues


def add_field_to_vector_file(sFilename_vector_in, aField, aValue):

    #open the file for updating
    pDataset = ogr.Open(sFilename_vector_in, update=1)
    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_vector_in} for updating")

    pLayer = pDataset.GetLayer(0)

    # Ensure fields and values are both lists or arrays
    if not isinstance(aField, (list, np.ndarray)):
        aField = [aField]

    if not isinstance(aValue, (list, np.ndarray, dict)):
        aValue = [aValue]

    # Get field names that already exist in the layer
    aField_existing = [field.GetName() for field in pLayer.schema]

    # Process each field
    aValue_row = aValue[0]
    for i, sFieldname in enumerate(aField):
        # Check if the field already exists
        if sFieldname not in aField_existing:
            # Determine the field type based on the value
            if isinstance(aValue_row, dict):
                sample_value = aValue_row[sFieldname][0] if len(aValue_row[sFieldname]) > 0 else None

            if sample_value is None:
                sFieldtype = ogr.OFTString
            elif isinstance(sample_value, int):
                sFieldtype = ogr.OFTInteger
            elif isinstance(sample_value, float):
                sFieldtype = ogr.OFTReal
            else:
                sFieldtype = ogr.OFTString

            # Create the field
            field_defn = ogr.FieldDefn(sFieldname, sFieldtype)
            pLayer.CreateField(field_defn)

    # Reset reading
    pLayer.ResetReading()

    # Update each pFeature with the values
    pFeature = pLayer.GetNextFeature()
    nrow = len(aValue)
    if nrow == 1:
        aValue_row = aValue[0]
        while pFeature:
            for j, sFieldname in enumerate(aField):
                if isinstance(aValue_row, dict):
                    pFeature.SetField(sFieldname, aValue_row[sFieldname])
            pLayer.SetFeature(pFeature)
            pFeature = pLayer.GetNextFeature()
    else:
        i = 0
        while pFeature:
            aValue_row = aValue[i]
            for j, sFieldname in enumerate(aField):
                if isinstance(aValue_row, dict):
                    pFeature.SetField(sFieldname, aValue_row[sFieldname])
            i += 1
            pLayer.SetFeature(pFeature)
            pFeature = pLayer.GetNextFeature()

    # Commit changes and close
    pDataset.FlushCache()
    pDataset = None

    return






