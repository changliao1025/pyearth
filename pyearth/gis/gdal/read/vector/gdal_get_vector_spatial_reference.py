import os
from osgeo import ogr


def gdal_get_vector_spatial_reference_wkt(sFilename_in: str) -> str:
    """Return the spatial reference of a vector file as a WKT string.

    Parameters
    ----------
    sFilename_in : str
        Path to the vector file.

    Returns
    -------
    str
        The spatial reference in Well Known Text (WKT) format.

    Raises
    ------
    FileNotFoundError
        If the vector file does not exist.
    RuntimeError
        If OGR cannot open the vector file.
    ValueError
        If the layer cannot be accessed or lacks spatial reference metadata.
    """

    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"File {sFilename_in} does not exist.")

    dataset = None
    try:
        dataset = ogr.Open(sFilename_in, 0)  # 0 = read-only
        if dataset is None:
            raise RuntimeError(f"Unable to open vector file {sFilename_in}.")

        layer = dataset.GetLayer(0)
        if layer is None:
            raise ValueError(f"Unable to access layer from {sFilename_in}.")

        spatial_ref = layer.GetSpatialRef()
        if spatial_ref is None:
            raise ValueError(
                f"Vector file {sFilename_in} does not define spatial reference metadata."
            )

        try:
            projection_wkt = spatial_ref.ExportToWkt()
        except Exception as e:
            raise ValueError(f"Unable to export spatial reference to WKT: {str(e)}")

        if not projection_wkt:
            raise ValueError(
                f"Vector file {sFilename_in} contains empty spatial reference metadata."
            )

        return projection_wkt

    finally:
        dataset = None
