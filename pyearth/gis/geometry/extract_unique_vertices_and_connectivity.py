import os, sys
import numpy as np

def extract_unique_vertices_and_connectivity(cell_x_coords_list, cell_y_coords_list):
    """
    Extract unique vertices and connectivity from mesh cell coordinates.

    Parameters:
    -----------
    cell_x_coords_list : list of arrays
        List where each element is an array of x-coordinates for one mesh cell
    cell_y_coords_list : list of arrays
        List where each element is an array of y-coordinates for one mesh cell

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
    """
    # Collect all vertices from all cells
    all_vertices_x = []
    all_vertices_y = []

    for x_coords, y_coords in zip(cell_x_coords_list, cell_y_coords_list):
        all_vertices_x.extend(x_coords)
        all_vertices_y.extend(y_coords)

    # Create coordinate pairs and find unique vertices
    vertex_pairs = list(zip(all_vertices_x, all_vertices_y))
    unique_vertex_pairs = list(set(vertex_pairs))

    # Extract unique x and y coordinates
    xv = np.array([pair[0] for pair in unique_vertex_pairs])
    yv = np.array([pair[1] for pair in unique_vertex_pairs])

    # Create vertex lookup dictionary for fast index mapping
    vertex_to_index = {vertex: idx for idx, vertex in enumerate(unique_vertex_pairs)}

    # Build connectivity array - maps each cell to its vertex indices
    connectivity_list = []
    max_vertices = 0

    for x_coords, y_coords in zip(cell_x_coords_list, cell_y_coords_list):
        cell_vertex_indices = []
        for x, y in zip(x_coords, y_coords):
            vertex_pair = (x, y)
            vertex_index = vertex_to_index[vertex_pair]
            cell_vertex_indices.append(vertex_index)
        connectivity_list.append(cell_vertex_indices)
        max_vertices = max(max_vertices, len(cell_vertex_indices))

    # Create 2D connectivity array padded with -1 for variable-sized polygons
    connectivity = np.full((len(connectivity_list), max_vertices), -1, dtype=np.int32)
    for i, cell_indices in enumerate(connectivity_list):
        connectivity[i, :len(cell_indices)] = cell_indices

    return xv, yv, connectivity, vertex_to_index