import numpy as np
from numpy.typing import NDArray

__all__ = ["remove_duplicate_closure"]


def remove_duplicate_closure(coords: NDArray[np.floating]) -> NDArray[np.floating]:
    """Remove duplicated closing vertex when first and last are equal."""
    if len(coords) > 1 and np.allclose(coords[0], coords[-1], atol=1.0e-12):
        return coords[:-1]
    return coords
