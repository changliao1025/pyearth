import os
from osgeo import gdal


def gdal_get_raster_spatial_reference_wkt(sFilename_in: str) -> str:
    """Return the spatial reference of a raster as a WKT string.

    Parameters
    ----------
    sFilename_in : str
        Path to the raster file readable by GDAL.

    Returns
    -------
    str
        The spatial reference in Well Known Text (WKT) format.

    Raises
    ------
    FileNotFoundError
        If the raster path does not exist.
    RuntimeError
        If GDAL cannot open the raster.
    ValueError
        If the raster lacks spatial reference information.
    """

    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"File {sFilename_in} does not exist.")

    dataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError(f"Unable to open raster {sFilename_in}.")

    try:
        spatial_ref = getattr(dataset, "GetSpatialRef", lambda: None)()
        if spatial_ref is not None:
            try:
                wkt = spatial_ref.ExportToWkt()
            except AttributeError:
                wkt = None
            else:
                if wkt:
                    return wkt

        wkt_fallback = dataset.GetProjectionRef() or dataset.GetProjection()
        if wkt_fallback:
            return wkt_fallback

        raise ValueError(
            f"Raster {sFilename_in} does not define spatial reference metadata."
        )
    finally:
        dataset = None
