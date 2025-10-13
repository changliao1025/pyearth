import os, sys
import numpy as np
from collections import defaultdict

def extract_unique_vertices_and_connectivity(aCell_vertex_x, aCell_vertex_y, dTolerance=1e-10):
    """
    Extract unique vertices and connectivity from mesh cell coordinates.

    Parameters:
    -----------
    aCell_vertex_longitude : list of arrays
        List where each element is an array of x-coordinates for one mesh cell
    aCell_vertex_latitude : list of arrays
        List where each element is an array of y-coordinates for one mesh cell
    tolerance : float, optional
        Tolerance for considering vertices as identical (default: 1e-10)

    Returns:
    --------
    xv : numpy.ndarray
        Array of unique x-coordinates (vertices)
    yv : numpy.ndarray
        Array of unique y-coordinates (vertices)
    connectivity : numpy.ndarray
        2D array where connectivity[i] contains vertex indices for cell i,
        padded with -1 for variable-sized polygons
    vertex_to_index : dict
        Dictionary mapping (x,y) coordinate pairs to vertex indices

    Raises:
    -------
    ValueError
        If input lists have different lengths or contain invalid data
    """
    # Input validation
    if len(aCell_vertex_x) != len(aCell_vertex_y):
        raise ValueError("aCell_vertex_longitude and aCell_vertex_latitude must have the same length")

    if not aCell_vertex_x:
        return np.array([]), np.array([]), np.array([]).reshape(0, 0), {}

    # Estimate total number of vertices for pre-allocation
    nVertex_total_estimate = sum(len(coords) for coords in aCell_vertex_x)

    # Use numpy arrays for better performance
    all_x = np.empty(nVertex_total_estimate, dtype=np.float64)
    all_y = np.empty(nVertex_total_estimate, dtype=np.float64)

    # Collect all vertices efficiently
    idx = 0
    cell_sizes = []
    for x_coords, y_coords in zip(aCell_vertex_x, aCell_vertex_y):
        x_arr = np.asarray(x_coords, dtype=np.float64)
        y_arr = np.asarray(y_coords, dtype=np.float64)

        if len(x_arr) != len(y_arr):
            raise ValueError(f"Mismatched coordinate array lengths in cell: {len(x_arr)} vs {len(y_arr)}")

        cell_size = len(x_arr)
        cell_sizes.append(cell_size)

        all_x[idx:idx + cell_size] = x_arr
        all_y[idx:idx + cell_size] = y_arr
        idx += cell_size

    # Trim arrays to actual size
    all_x = all_x[:idx]
    all_y = all_y[:idx]

    # Create coordinate pairs as structured array for efficient comparison
    vertex_dtype = np.dtype([('x', np.float64), ('y', np.float64)])
    all_vertices = np.empty(len(all_x), dtype=vertex_dtype)
    all_vertices['x'] = all_x
    all_vertices['y'] = all_y

    # Find unique vertices using numpy (much faster than set operations)
    unique_vertices, inverse_indices = np.unique(all_vertices, return_inverse=True)

    # Extract coordinates
    xv = unique_vertices['x']
    yv = unique_vertices['y']

    # Create vertex lookup dictionary
    vertex_to_index = {(float(v['x']), float(v['y'])): i for i, v in enumerate(unique_vertices)}

    # Build connectivity efficiently
    max_vertices = max(cell_sizes) if cell_sizes else 0
    num_cells = len(aCell_vertex_x)
    connectivity = np.full((num_cells, max_vertices), -1, dtype=np.int32)

    # Use inverse indices to map vertices to indices
    start_idx = 0
    for i, cell_size in enumerate(cell_sizes):
        connectivity[i, :cell_size] = inverse_indices[start_idx:start_idx + cell_size]
        start_idx += cell_size

    return xv, yv, connectivity, vertex_to_index