import numpy as np
from typing import Union, Optional


def remap(
    x: Union[float, int],
    oMin: Union[float, int],
    oMax: Union[float, int],
    nMin: Union[float, int],
    nMax: Union[float, int],
) -> Optional[float]:
    """
    Remap a value from one range to another range.

    Linearly maps a value x from the input range [oMin, oMax] to the output range [nMin, nMax].
    Handles reversed ranges (where min > max) automatically.

    Parameters
    ----------
    x : Union[float, int]
        The value to remap.
    oMin : Union[float, int]
        Minimum value of the input range.
    oMax : Union[float, int]
        Maximum value of the input range.
    nMin : Union[float, int]
        Minimum value of the output range.
    nMax : Union[float, int]
        Maximum value of the output range.

    Returns
    -------
    Optional[float]
        The remapped value, or None if input/output ranges are invalid (zero range).

    Notes
    -----
    The function automatically detects and handles reversed ranges. For example:
    - Input range [0, 10] with output range [100, 0] will reverse the mapping
    - Input range [10, 0] with output range [0, 100] will also reverse appropriately

    If the input or output range has zero span (min == max), the function returns None
    to avoid division by zero.

    Examples
    --------
    Basic linear remapping:

    >>> remap(5, 0, 10, 0, 100)  # Maps 5 from [0,10] to [0,100]
    50.0

    Reverse output range:

    >>> remap(5, 0, 10, 100, 0)  # Maps 5 from [0,10] to [100,0] (reversed)
    50.0

    Reverse input range:

    >>> remap(5, 10, 0, 0, 100)  # Maps 5 from [10,0] to [0,100]
    50.0
    """
    # Input validation
    if oMin == oMax:
        print("Warning: Zero input range")
        return None
    if nMin == nMax:
        print("Warning: Zero output range")
        return None

    # Determine if ranges are reversed
    input_reversed = oMin > oMax
    output_reversed = nMin > nMax

    # Normalize ranges for calculation
    old_min = min(oMin, oMax)
    old_max = max(oMin, oMax)
    new_min = min(nMin, nMax)
    new_max = max(nMin, nMax)

    # Calculate the normalized position in the input range
    if input_reversed:
        # For reversed input range, invert the position
        normalized_pos = (old_max - x) / (old_max - old_min)
    else:
        normalized_pos = (x - old_min) / (old_max - old_min)

    # Map to output range
    if output_reversed:
        # For reversed output range, invert the result
        result = new_max - normalized_pos * (new_max - new_min)
    else:
        result = normalized_pos * (new_max - new_min) + new_min

    return float(result)
