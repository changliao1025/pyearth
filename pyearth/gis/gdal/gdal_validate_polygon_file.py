import os
from osgeo import ogr


def gdal_validate_polygon_file(filename: str) -> bool:
    """
    Check if a given polygon file is valid.

    Parameters:
    filename (str): The filename of the polygon file.

    Returns:
    bool: True if the file is a valid polygon file, False otherwise.
    """
    if not os.path.exists(filename):
        return False

    try:
        dataset = ogr.Open(filename)
        if dataset is None:
            return False

        layer = dataset.GetLayer()
        if layer is None:
            return False

        for feature in layer:
            geometry = feature.GetGeometryRef()
            if geometry is None:
                return False

            if geometry.GetGeometryName() != "POLYGON":
                return False

            if not geometry.IsValid():
                return False

        return True

    except Exception as e:
        print(f"An error occurred: {e}")
        return False
