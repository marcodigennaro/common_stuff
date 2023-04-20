import numpy as np
from typing import List, Any, Union

Vector = List[float]


def unit_vector(my_array: Union[Vector, Any]) -> Vector:
    """
    Returns the unit vector of the vector.
    - Input: my_array is a N-dim vector
    - Output: normalized vector from my_array
    - Example:
        unit_vector([1,1,1]) = [ 0.5773, 0.5773, 0.5773 ]
    """
    norm = float(np.linalg.norm(my_array))
    norm_array = [elem/norm for elem in my_array]
    return norm_array


def angle_between(v1: Vector, v2: Vector) -> Any:
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
