
#from math import acos
import numpy as np
def calculate_angle_between_vectors_degrees(u, v):
    """Return the angle between two vectors in any dimension space,
    in degrees.
    """
    a = np.dot(u, v)
    b = np.linalg.norm(u)
    c = np.linalg.norm(v)
    d = a / (b * c)
    if d > 1:
        d = 1
    if d < -1:
        d = -1

    #e = acos(d) 
    e = np.arccos(d)
    
    f = np.degrees(e)

    return f