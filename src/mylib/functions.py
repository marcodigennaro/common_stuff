import numpy as np


def unit_vector(my_array: np.ndarray) -> np.ndarray:
    """
    Returns the unit vector of the vector.
    - Input: a N-dim vector
    - Output: scalar
    - Example:
        unit_vector([1,1,1]) = 1.7320508075688772
    """
    return my_array / np.linalg.norm(my_array)


def angle_between(v1: np.array , v2: np.array ) -> float:
    """
    Returns the angle in radians between vectors 'v1' and 'v2':
    - Inputs: two vectors
    - Output: scalar
    - Example:
        angle_between((1, 0, 0), (1, 0, 0)) = 0.0
        angle_between((1, 0, 0), (0, 1, 0)) = 1.5707963267948966
        angle_between((1, 0, 0), (-1, 0, 0)) = 3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))



