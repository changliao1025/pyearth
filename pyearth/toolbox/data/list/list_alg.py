"""List min and max algorithms ignoring NaNs."""

from typing import Any, Sequence
import numpy as np


def list_max(pList: Sequence[Any]) -> float:
    """
    Find the maximum value across a nested sequence of array-like elements, ignoring NaNs.

    Parameters
    ----------
    pList : sequence
        Sequence of array-like sublists.

    Returns
    -------
    float
        Maximum value found, or np.nan if no valid numbers are present.
    """
    dmax_max = np.nan
    for sublist in pList:
        arr = np.asarray(sublist)
        if arr.size == 0:
            continue
        valid = arr[np.isfinite(arr)]
        if valid.size > 0:
            dmax = float(np.max(valid))
            if np.isnan(dmax_max) or dmax > dmax_max:
                dmax_max = dmax
    return dmax_max


def list_min(pList: Sequence[Any]) -> float:
    """
    Find the minimum value across a nested sequence of array-like elements, ignoring NaNs.

    Parameters
    ----------
    pList : sequence
        Sequence of array-like sublists.

    Returns
    -------
    float
        Minimum value found, or np.nan if no valid numbers are present.
    """
    dmin_min = np.nan
    for sublist in pList:
        arr = np.asarray(sublist)
        if arr.size == 0:
            continue
        valid = arr[np.isfinite(arr)]
        if valid.size > 0:
            dmin = float(np.min(valid))
            if np.isnan(dmin_min) or dmin < dmin_min:
                dmin_min = dmin
    return dmin_min
