import logging
from typing import Optional, Union

logger = logging.getLogger(__name__)


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
        The remapped value, or None if the input range is zero (oMin == oMax).
    """
    if oMin == oMax:
        logger.warning("Zero input range in remap: oMin == oMax (%s)", oMin)
        return None

    if nMin == nMax:
        return float(nMin)

    # Standard affine linear interpolation naturally handles inverted ranges
    result = nMin + (float(x) - float(oMin)) * (float(nMax) - float(nMin)) / (float(oMax) - float(oMin))
    return float(result)
