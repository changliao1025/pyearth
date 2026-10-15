"""Calculate angle between two vectors in degrees."""

import numpy as np

__all__ = ["calculate_angle_between_vectors_degrees"]


def calculate_angle_between_vectors_degrees(
    vector1: np.ndarray, vector2: np.ndarray
) -> float:
    """Return the angle between two vectors in any dimensional space, in degrees.

    Parameters
    ----------
    vector1 : np.ndarray
        First vector.
    vector2 : np.ndarray
        Second vector.

    Returns
    -------
    float
        Angle between the two vectors in degrees in range [0, 180].
        Returns 0.0 if either vector has zero length.
    """
    v1 = np.asarray(vector1, dtype=float)
    v2 = np.asarray(vector2, dtype=float)

    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    norm_product = norm1 * norm2

    if norm_product == 0.0 or np.isclose(norm_product, 0.0, atol=1e-15):
        return 0.0

    dot_product = np.dot(v1, v2)
    cosine_similarity = np.clip(dot_product / norm_product, -1.0, 1.0)
    angle_in_radians = np.arccos(cosine_similarity)
    return float(np.degrees(angle_in_radians))
