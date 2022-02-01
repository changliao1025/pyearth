import numpy as np
def calculate_hexagon_area(dLength_edge_in):
    """
    Calculate the area of a hexagon grid

    Args:
        dLength_edge_in (float): The length of hexagon edge

    Returns:
        float: The hexagon area
    """
    area = 1.5 * np.sqrt(3.0) * dLength_edge_in
    return area