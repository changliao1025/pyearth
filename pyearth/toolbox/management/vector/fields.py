import numpy as np
from osgeo import ogr, gdal
from pyearth.gis.gdal.gdal_vector_format_support import (
    get_vector_driver_from_extension,
    get_vector_format_from_extension
)
from pyearth.gis.gdal.gdal_vector_format_support import (
    check_parquet_support,
    has_parquet_support,
)

def get_field_value(sFilename_vector_in, sField_name_in, dMissing_value = None):
    """Return values for a field from the first (or named) layer of a vector.

    This helper opens a vector datasource, finds the requested field (case-
    insensitive), and returns the values as a NumPy array. It supports an
    optional missing value filter and will provide a clear error message if
    the field is not found or the file cannot be opened.

    Parameters
    ----------
    sFilename_vector_in : str
        Path to the vector datasource.
    sField_name_in : str
        Field name to extract (case-insensitive).
    dMissing_value : optional
        If provided, filter out values equal to this missing value.

    Returns
    -------
    numpy.ndarray
        Array of field values (dtype inferred). Empty array if no values found.
    """
    # Open the vector file
    pDataset = None
    # Prefer opening with a resolved driver (handles Parquet/Arrow resolution)
    try:
        drv = get_vector_driver_from_extension(sFilename_vector_in)
        if drv is not None:
            pDataset = drv.Open(sFilename_vector_in, 0)
    except Exception:
        # fallback to generic open
        pDataset = None

    if pDataset is None:
        try:
            pDataset = ogr.Open(sFilename_vector_in)
        except Exception as e:
            # try to provide a helpful hint based on extension
            try:
                fmt = get_vector_format_from_extension(sFilename_vector_in)
                if fmt and fmt.lower().startswith('parquet') and not has_parquet_support():
                    raise ValueError(
                        f"Could not open {sFilename_vector_in}: GDAL built without Parquet/Arrow support."
                    )
            except Exception:
                pass
            raise ValueError(f"Could not open {sFilename_vector_in}: {e}")

    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_vector_in}")

    # Get the first pLayer
    pLayer = pDataset.GetLayer(0)

    # Resolve field name case-insensitively
    layer_fields = [f.name for f in pLayer.schema]
    field_map = {f.lower(): f for f in layer_fields}
    key = sField_name_in.lower()
    if key not in field_map:
        raise KeyError(f"Field '{sField_name_in}' not found in layer. Available fields: {layer_fields}")
    field_name_actual = field_map[key]

    # Collect values
    aValues = []
    for pFeature in pLayer:
        val = pFeature.GetField(field_name_actual)
        # Treat empty strings and None as missing if dMissing_value is provided
        if dMissing_value is not None:
            try:
                if val == dMissing_value:
                    continue
            except Exception:
                pass
        aValues.append(val)


    # convert to numpy with inferred dtype
    if len(aValues) == 0:
        return np.array([])
    try:
        arr = np.array(aValues)
    except Exception:
        # fallback to object array
        arr = np.array(aValues, dtype=object)
    return arr

def get_field_and_value(sFilename_vector_in):
    # Open the vector file
    pDataset = None
    try:
        drv = get_vector_driver_from_extension(sFilename_vector_in)
        if drv is not None:
            pDataset = drv.Open(sFilename_vector_in, 0)
    except Exception:
        pDataset = None

    if pDataset is None:
        try:
            pDataset = ogr.Open(sFilename_vector_in)
        except Exception as e:
            # Provide a helpful hint for parquet files
            try:
                fmt = get_vector_format_from_extension(sFilename_vector_in)
                if fmt and fmt.lower().startswith('parquet'):
                    print("Parquet support may require GDAL built with Arrow/Parquet driver.")
                    print("Available parquet drivers:", check_parquet_support())
            except Exception:
                pass
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






