def find_index_by_longitude_latitude(dLon_in,
                         dLat_in,
                         dLon_min_in,
                         dLat_max_in,
                         dResolution_x_in,
                         dResolution_y_in,
                         iFlag_center_in=None):
    """
    Find the index of a lcoation in a large 2D array, we requre the location is at center, the longitude should be in the same range (-180, 180) or (0, 360)

    Args:
        dLon_in (float): _description_
        dLat_in (float): _description_
        dLon_min_in (float): _description_
        dLat_max_in (float): _description_
        dResolution_x_in (float): The mesh cell resolution in x direction
        dResolution_y_in (float): The mesh cell resolution in y directioniption_
        iFlag_center_in (int, optional): Whether the max/min is at center. Defaults to None.
    """
    if iFlag_center_in is None:
        iFlag_center = 1
    else:
        iFlag_center = iFlag_center_in

    if iFlag_center == 1:  # when the max and min are center
        iRow_index = round((dLat_max_in - dLat_in) / dResolution_y_in)
        iColumn_index = round((dLon_in - (dLon_min_in)) / dResolution_x_in)
        pass
    else:  # when the max and min are upper corner
        iRow_index = round(
            ((dLat_max_in - 0.5*dResolution_y_in) - dLat_in) / dResolution_y_in)
        iColumn_index = round(
            (dLon_in - (dLon_min_in + 0.5*dResolution_x_in)) / dResolution_x_in)
        pass

    return iRow_index, iColumn_index
