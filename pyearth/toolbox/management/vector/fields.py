import numpy as np
from osgeo import ogr, gdal

def check_parquet_support():
    """
    Check if Parquet/Arrow driver is available in OGR
    """
    driver_names = [ogr.GetDriver(i).GetName() for i in range(ogr.GetDriverCount())]
    parquet_drivers = [name for name in driver_names if 'parquet' in name.lower() or 'arrow' in name.lower()]
    return parquet_drivers

def get_supported_formats():
    """
    Get list of all supported OGR formats
    """
    formats = {}
    for i in range(ogr.GetDriverCount()):
        driver = ogr.GetDriver(i)
        formats[driver.GetName()] = driver.GetMetadata()
    return formats

def diagnose_parquet_issue():
    """
    Diagnose parquet support issues and suggest solutions
    """
    print("=== GDAL/OGR Parquet Support Diagnosis ===")
    print(f"GDAL Version: {gdal.VersionInfo()}")

    parquet_drivers = check_parquet_support()
    if parquet_drivers:
        print(f"Available Parquet drivers: {parquet_drivers}")
    else:
        print("No Parquet drivers found!")
        print("\nSolutions:")
        print("1. Install GDAL with Arrow/Parquet support:")
        print("   conda install -c conda-forge gdal>=3.5 pyarrow")
        print("2. Or convert parquet to supported format:")
        print("   Use geopandas: gdf.to_file('file.gpkg', driver='GPKG')")
        print("3. Check available drivers:")

        print("Available drivers:")
        for i in range(min(10, ogr.GetDriverCount())):
            driver = ogr.GetDriver(i)
            print(f"  - {driver.GetName()}")
        print("  ... (showing first 10)")

    return parquet_drivers
def get_field_value(sFilename_vector_in, sField_name_in, dMissing_value = None):
    # Open the vector file
    try:
        pDataset = ogr.Open(sFilename_vector_in)
    except Exception as e:
        if sFilename_vector_in.lower().endswith('.parquet'):
            print("Parquet support may require GDAL with Arrow/Parquet driver.")
            print("Available parquet drivers:", check_parquet_support())
        raise ValueError(f"Could not open {sFilename_vector_in}: {e}")

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
    try:
        pDataset = ogr.Open(sFilename_vector_in)
    except Exception as e:
        if sFilename_vector_in.lower().endswith('.parquet'):
            print("Parquet support may require GDAL with Arrow/Parquet driver.")
            print("Available parquet drivers:", check_parquet_support())
        raise ValueError(f"Could not open {sFilename_vector_in}: {e}")

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
    try:
        pDataset = ogr.Open(sFilename_vector_in, update=1)
    except Exception as e:
        # Try alternative approach for parquet files
        if sFilename_vector_in.lower().endswith('.parquet'):
            print(f"Warning: Direct parquet update not supported. Error: {e}")
            print("Consider converting to a different format for field updates.")
            return
        else:
            raise ValueError(f"Could not open {sFilename_vector_in} for updating: {e}")

    if pDataset is None:
        # Check if it's a parquet file
        if sFilename_vector_in.lower().endswith('.parquet'):
            print("Warning: Parquet format may not support in-place updates.")
            print("Available drivers:", [ogr.GetDriver(i).GetName() for i in range(ogr.GetDriverCount())])
            return
        else:
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






