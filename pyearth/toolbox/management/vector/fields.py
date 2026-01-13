import numpy as np
import logging
from typing import List, Tuple, Union, Optional, Any, Dict
from osgeo import ogr, gdal
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_format_from_filename


def get_field_value(
    sFilename_vector_in: str,
    sField_name_in: str,
    dMissing_value: Optional[Any] = None,
    iLayer: Optional[Union[int, str]] = None,
    bVerbose: bool = False,
    bUnique: bool = False,
) -> np.ndarray:
    """Return values for a field from a specified layer of a vector.

    This helper opens a vector datasource, finds the requested field (case-
    insensitive), and returns the values as a NumPy array. It supports optional
    missing value filtering, layer selection, verbose logging, and unique value
    extraction.

    Parameters
    ----------
    sFilename_vector_in : str
        Path to the vector datasource.
    sField_name_in : str
        Field name to extract (case-insensitive).
    dMissing_value : optional
        If provided, filter out values equal to this missing value.
    iLayer : int or str, optional
        Layer index (int) or name (str) to read from. Defaults to first layer (0).
    bVerbose : bool, optional
        If True, log progress information. Defaults to False.
    bUnique : bool, optional
        If True, return only unique values. Defaults to False.

    Returns
    -------
    numpy.ndarray
        Array of field values (dtype inferred). Empty array if no values found.
        If bUnique is True, returns unique values.
    """
    # Open the vector file
    pDataset = None
    if bVerbose:
        logging.info(f"Opening vector file: {sFilename_vector_in}")
    try:
        drv = get_vector_format_from_filename(sFilename_vector_in)
        if drv is not None:
            pDataset = drv.Open(sFilename_vector_in, 0)
    except Exception:
        pDataset = None

    if pDataset is None:
        try:
            pDataset = ogr.Open(sFilename_vector_in)
        except Exception as e:
            raise ValueError(f"Could not open {sFilename_vector_in}: {e}")

    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_vector_in}")

    # Get the specified layer
    if iLayer is not None:
        if isinstance(iLayer, int):
            pLayer = pDataset.GetLayer(iLayer)
        elif isinstance(iLayer, str):
            pLayer = pDataset.GetLayerByName(iLayer)
        else:
            raise ValueError("iLayer must be int or str")
        if pLayer is None:
            raise ValueError(f"Layer {iLayer} not found in {sFilename_vector_in}")
    else:
        pLayer = pDataset.GetLayer(0)

    if bVerbose:
        logging.info(f"Reading layer: {pLayer.GetName()}")

    # Resolve field name case-insensitively
    if bVerbose:
        logging.info(f"Resolving field: {sField_name_in}")
    layer_fields = [f.name for f in pLayer.schema]
    field_map = {f.lower(): f for f in layer_fields}
    key = sField_name_in.lower()
    if key not in field_map:
        raise KeyError(
            f"Field '{sField_name_in}' not found in layer. Available fields: {layer_fields}"
        )
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

    if bVerbose:
        logging.info(f"Collected {len(aValues)} values")

    # convert to numpy with inferred dtype
    if len(aValues) == 0:
        arr = np.array([])
    else:
        try:
            arr = np.array(aValues)
        except Exception:
            # fallback to object array
            arr = np.array(aValues, dtype=object)

    if bUnique:
        arr = np.unique(arr)

    if bVerbose:
        logging.info(f"Returning array with shape: {arr.shape}, dtype: {arr.dtype}")

    return arr


def get_field_and_value(
    sFilename_vector_in: str,
    iLayer: Optional[Union[int, str]] = None,
    bVerbose: bool = False,
) -> Tuple[List[str], List[Dict[str, Any]]]:
    """Return all field names and values from a specified layer of a vector.

    This helper opens a vector datasource and extracts all field names and their
    corresponding values for each feature in the layer. It supports optional layer
    selection and verbose logging.

    Parameters
    ----------
    sFilename_vector_in : str
        Path to the vector datasource.
    iLayer : int or str, optional
        Layer index (int) or name (str) to read from. Defaults to first layer (0).
    bVerbose : bool, optional
        If True, log progress information. Defaults to False.

    Returns
    -------
    Tuple[List[str], List[Dict[str, Any]]]
        A tuple containing:
        - List of field names (strings)
        - List of dictionaries, one per feature, with field names as keys and values as values
    """
    pDataset = None
    try:
        drv = get_vector_format_from_filename(sFilename_vector_in)
        if drv is not None:
            pDataset = drv.Open(sFilename_vector_in, 0)
    except Exception:
        pDataset = None

    if pDataset is None:
        try:
            pDataset = ogr.Open(sFilename_vector_in)
        except Exception as e:
            raise ValueError(f"Could not open {sFilename_vector_in}: {e}")

    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_vector_in}")

    # Get the first pLayer
    pLayer = pDataset.GetLayer(0)

    # Collect all field names and values
    aField_names = [field.name for field in pLayer.schema]
    aValues = list()

    for pFeature in pLayer:
        feature_dict = {}
        # Collect all field values for this feature
        for sFieldname in aField_names:
            feature_dict[sFieldname] = pFeature.GetField(sFieldname)

        # Add this feature's dictionary to the list
        aValues.append(feature_dict)

    return aField_names, aValues


def add_field_to_vector_file(
    sFilename_vector_in: str,
    aField: Union[str, List[str]],
    aValue: Union[Any, List[Any], Dict[str, Any], List[Dict[str, Any]]],
    iLayer: Optional[Union[int, str]] = None,
    bVerbose: bool = False,
) -> None:
    """Add fields to a vector file with specified values.

    This function opens a vector datasource for updating, creates new fields if they
    don't exist, and populates them with the provided values. It supports single or
    multiple fields, various value formats, optional layer selection, and verbose logging.

    Parameters
    ----------
    sFilename_vector_in : str
        Path to the vector datasource to update.
    aField : str or List[str]
        Field name(s) to add. If str, treated as single field.
    aValue : Any or List[Any] or Dict[str, Any] or List[Dict[str, Any]]
        Values to assign. Format depends on number of features and fields.
        - Single value: applied to all features for single field
        - List: one value per feature
        - Dict: field names as keys, values as lists per feature
        - List of dicts: one dict per feature with field values
    iLayer : int or str, optional
        Layer index (int) or name (str) to update. Defaults to first layer (0).
    bVerbose : bool, optional
        If True, log progress information. Defaults to False.

    Returns
    -------
    None
        Modifies the file in place.

    Notes
    -----
    - For Parquet files, in-place updates may not be supported.
    - Field types are inferred from the first value encountered.
    """
    # Open the vector file for updating
    pDataset = None
    if bVerbose:
        logging.info(f"Opening vector file for updating: {sFilename_vector_in}")
    try:
        drv = get_vector_format_from_filename(sFilename_vector_in)
        if drv is not None:
            pDataset = drv.Open(sFilename_vector_in, 1)  # 1 for update
    except Exception:
        pDataset = None

    if pDataset is None:
        try:
            pDataset = ogr.Open(sFilename_vector_in, 1)  # 1 for update
        except Exception as e:
            # Handle parquet files specially
            if sFilename_vector_in.lower().endswith(".parquet"):
                if bVerbose:
                    logging.warning(
                        f"Parquet files may not support in-place updates: {e}"
                    )
                return
            raise ValueError(f"Could not open {sFilename_vector_in} for updating: {e}")

    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_vector_in} for updating")

    # Get the specified layer
    if iLayer is not None:
        if isinstance(iLayer, int):
            pLayer = pDataset.GetLayer(iLayer)
        elif isinstance(iLayer, str):
            pLayer = pDataset.GetLayerByName(iLayer)
        else:
            raise ValueError("iLayer must be int or str")
        if pLayer is None:
            raise ValueError(f"Layer {iLayer} not found in {sFilename_vector_in}")
    else:
        pLayer = pDataset.GetLayer(0)

    if bVerbose:
        logging.info(f"Updating layer: {pLayer.GetName()}")

    # Ensure fields and values are both lists or arrays
    if not isinstance(aField, (list, np.ndarray)):
        aField = [aField]

    if not isinstance(aValue, (list, np.ndarray, dict)):
        aValue = [aValue]

    # Get field names that already exist in the layer
    aField_existing = [field.GetName() for field in pLayer.schema]

    if bVerbose:
        logging.info(f"Existing fields: {aField_existing}")
        logging.info(f"Adding fields: {aField}")

    # Process each field
    aValue_row = aValue[0] if len(aValue) > 0 else None
    for i, sFieldname in enumerate(aField):
        # Check if the field already exists
        if sFieldname not in aField_existing:
            # Determine the field type based on the value
            if isinstance(aValue_row, dict):
                sample_value = (
                    aValue_row.get(sFieldname, [None])[0]
                    if isinstance(aValue_row.get(sFieldname), list)
                    and len(aValue_row.get(sFieldname, [])) > 0
                    else aValue_row.get(sFieldname)
                )
            else:
                sample_value = aValue_row

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
            if bVerbose:
                logging.info(f"Created field '{sFieldname}' of type {sFieldtype}")

    # Reset reading
    pLayer.ResetReading()

    # Update each feature with the values
    pFeature = pLayer.GetNextFeature()
    nrow = len(aValue)
    if bVerbose:
        logging.info(f"Updating {nrow} features")

    if nrow == 1:
        aValue_row = aValue[0]
        while pFeature:
            for j, sFieldname in enumerate(aField):
                if isinstance(aValue_row, dict):
                    value = aValue_row.get(sFieldname)
                    if isinstance(value, list) and len(value) > 0:
                        pFeature.SetField(sFieldname, value[0])
                    else:
                        pFeature.SetField(sFieldname, value)
                else:
                    pFeature.SetField(sFieldname, aValue_row)
            pLayer.SetFeature(pFeature)
            pFeature = pLayer.GetNextFeature()
    else:
        i = 0
        while pFeature and i < nrow:
            aValue_row = aValue[i]
            for j, sFieldname in enumerate(aField):
                if isinstance(aValue_row, dict):
                    value = aValue_row.get(sFieldname)
                    if isinstance(value, list) and len(value) > 0:
                        pFeature.SetField(sFieldname, value[0])
                    else:
                        pFeature.SetField(sFieldname, value)
                else:
                    pFeature.SetField(sFieldname, aValue_row)
            pLayer.SetFeature(pFeature)
            pFeature = pLayer.GetNextFeature()
            i += 1

    # Commit changes and close
    pDataset.FlushCache()
    pDataset = None

    if bVerbose:
        logging.info("Field addition completed")
