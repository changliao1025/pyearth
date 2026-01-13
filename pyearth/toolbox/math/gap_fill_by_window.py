import statistics
from typing import Optional, Union, Tuple
import numpy as np


def search_neighbors(
    iIndex_in: int, aArray_in: np.ndarray, iWindow_size_in: Optional[int] = None
) -> np.ndarray:
    """
    Extract neighboring values around a given index for gap filling operations.

    This function extracts a window of neighboring values centered on the specified
    index. For indices near the array boundaries, it returns all available neighbors
    within the specified window size, ensuring robust gap filling at edges.

    Parameters
    ----------
    iIndex_in : int
        The center index around which to extract neighbors. Must be a valid
        index within the input array bounds.
    aArray_in : np.ndarray
        Input 1D array from which to extract neighboring values. Should contain
        numeric data that may include NaN values for gap filling.
    iWindow_size_in : int, optional
        The number of neighbors to extract on each side of the center index.
        If None, defaults to 1. Must be non-negative.

    Returns
    -------
    np.ndarray
        Array of neighboring values. The length will be 2*window_size + 1 for
        interior points, or fewer for boundary points where neighbors are limited.
        Values are returned in their original order (left to right around center).

    Raises
    ------
    TypeError
        If input array is not a numpy array or index/window_size are not integers.
    ValueError
        If index is out of bounds or window_size is negative.
    IndexError
        If the input array is empty.

    Notes
    -----
    - For boundary handling: when the requested window extends beyond array bounds,
      only available neighbors are returned. This ensures gap filling works
      consistently across the entire array.
    - The returned array maintains the original data types and NaN values.
    - This function is designed for use in gap filling algorithms where spatial
      locality (neighboring values) is used to estimate missing data.

    Examples
    --------
    Extract 1 neighbor on each side (3 total values):

    >>> arr = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    >>> neighbors = search_neighbors(2, arr, 1)
    >>> print(neighbors)  # [2.0, nan, 4.0]

    Handle boundary case with limited neighbors:

    >>> neighbors = search_neighbors(0, arr, 2)
    >>> print(neighbors)  # [1.0, 2.0] (only 2 values available)
    """
    # Input validation
    if not isinstance(aArray_in, np.ndarray):
        raise TypeError(f"Input array must be a numpy array, got {type(aArray_in)}")

    if aArray_in.ndim != 1:
        raise ValueError(
            f"Input array must be 1-dimensional, got {aArray_in.ndim}D array"
        )

    if len(aArray_in) == 0:
        raise IndexError("Input array cannot be empty")

    if not isinstance(iIndex_in, int):
        raise TypeError(f"Index must be an integer, got {type(iIndex_in)}")

    if not (0 <= iIndex_in < len(aArray_in)):
        raise ValueError(
            f"Index {iIndex_in} is out of bounds for array of length {len(aArray_in)}"
        )

    # Handle window size
    if iWindow_size_in is None:
        iWindow_size = 1
    elif not isinstance(iWindow_size_in, int):
        raise TypeError(
            f"Window size must be an integer or None, got {type(iWindow_size_in)}"
        )
    elif iWindow_size_in < 0:
        raise ValueError(f"Window size must be non-negative, got {iWindow_size_in}")
    else:
        iWindow_size = iWindow_size_in

    # Calculate the actual window bounds, respecting array boundaries
    array_length = len(aArray_in)
    start_idx = max(0, iIndex_in - iWindow_size)
    end_idx = min(array_length, iIndex_in + iWindow_size + 1)

    # Extract the neighboring values
    neighbors = aArray_in[start_idx:end_idx]

    return neighbors


def gap_fill_by_window(
    aArray_in: np.ndarray,
    iWindow_size_in: Optional[int] = None,
    min_good_neighbors: Optional[int] = None,
) -> np.ndarray:
    """
    Fill gaps (NaN values) in a 1D array using neighboring values within a sliding window.

    This function performs gap filling by replacing NaN values with the median of
    neighboring non-NaN values within a specified window. Only gaps with sufficient
    neighboring data are filled to ensure reliable interpolation.

    Parameters
    ----------
    aArray_in : np.ndarray
        Input 1D array containing numeric data with potential NaN gaps to be filled.
    iWindow_size_in : int, optional
        Size of the search window (number of neighbors on each side of the gap).
        If None, defaults to 1. Must be non-negative.
    min_good_neighbors : int, optional
        Minimum number of non-NaN neighboring values required to perform gap filling.
        If None, defaults to window_size + 1. Must be positive.

    Returns
    -------
    np.ndarray
        Array with gaps filled where possible. The original array is not modified;
        a copy is returned with filled values.

    Raises
    ------
    TypeError
        If input array is not a numpy array or parameters are of wrong type.
    ValueError
        If array is not 1D, window size is negative, or min_good_neighbors is invalid.

    Notes
    -----
    - Gap filling only occurs when there are enough non-NaN neighbors within the window.
    - The median of neighboring values is used for robust gap filling.
    - Edge gaps may have fewer available neighbors, which is handled automatically.
    - If insufficient neighbors are available, the gap remains unfilled.
    - The function processes all gaps in the array, not just interior points.

    Examples
    --------
    Fill gaps in a simple array:

    >>> arr = np.array([1.0, np.nan, 3.0, np.nan, 5.0])
    >>> filled = gap_fill_by_window(arr, iWindow_size_in=1, min_good_neighbors=1)
    >>> print(filled)  # [1.0, 2.0, 3.0, 4.0, 5.0]

    Handle case with insufficient neighbors:

    >>> arr = np.array([np.nan, 2.0, 3.0])
    >>> filled = gap_fill_by_window(arr, iWindow_size_in=1, min_good_neighbors=2)
    >>> print(filled)  # [nan, 2.0, 3.0] (first gap unfilled due to insufficient neighbors)
    """
    # Input validation
    if not isinstance(aArray_in, np.ndarray):
        raise TypeError(f"Input array must be a numpy array, got {type(aArray_in)}")

    if aArray_in.ndim != 1:
        raise ValueError(
            f"Input array must be 1-dimensional, got {aArray_in.ndim}D array"
        )

    if len(aArray_in) == 0:
        raise ValueError("Input array cannot be empty")

    # Handle window size
    if iWindow_size_in is None:
        iWindow_size = 1
    elif not isinstance(iWindow_size_in, int):
        raise TypeError(
            f"Window size must be an integer or None, got {type(iWindow_size_in)}"
        )
    elif iWindow_size_in < 0:
        raise ValueError(f"Window size must be non-negative, got {iWindow_size_in}")
    else:
        iWindow_size = iWindow_size_in

    # Handle minimum good neighbors
    if min_good_neighbors is None:
        min_good_neighbors = iWindow_size + 1
    elif not isinstance(min_good_neighbors, int):
        raise TypeError(
            f"min_good_neighbors must be an integer or None, got {type(min_good_neighbors)}"
        )
    elif min_good_neighbors <= 0:
        raise ValueError(
            f"min_good_neighbors must be positive, got {min_good_neighbors}"
        )

    # Create a copy to avoid modifying the original
    aArray_out = aArray_in.copy()

    # Find all NaN positions
    nan_mask = np.isnan(aArray_out)
    nan_indices = np.where(nan_mask)[0]

    # Process each gap
    for gap_idx in nan_indices:
        # Get neighboring values
        neighbors = search_neighbors(gap_idx, aArray_out, iWindow_size)

        # Count good (non-NaN) neighbors
        good_neighbors = neighbors[~np.isnan(neighbors)]
        n_good = len(good_neighbors)

        # Only fill if we have enough good neighbors
        if n_good >= min_good_neighbors:
            try:
                # Use median for robust gap filling
                fill_value = np.median(good_neighbors)
                aArray_out[gap_idx] = fill_value
            except (ValueError, RuntimeWarning):
                # Handle edge cases where median calculation might fail
                continue

    return aArray_out
