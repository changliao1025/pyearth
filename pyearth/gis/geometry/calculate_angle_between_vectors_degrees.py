import numpy as np


def calculate_angle_between_vectors_degrees(
    vector1: np.ndarray, vector2: np.ndarray
) -> float:
    """
    Return the angle between two vectors in any dimension space, in degrees.

    Parameters:
    vector1 (np.ndarray): The first vector.
    vector2 (np.ndarray): The second vector.

    Returns:
    float: The angle between the two vectors in degrees.
    """
    dot_product = np.dot(vector1, vector2)
    norm_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    cosine_similarity = np.clip(dot_product / norm_product, -1.0, 1.0)
    angle_in_radians = np.arccos(cosine_similarity)
    angle_in_degrees = np.degrees(angle_in_radians)

    return angle_in_degrees
