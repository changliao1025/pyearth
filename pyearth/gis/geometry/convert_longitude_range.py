import numpy as np
def convert_360_to_180(dLongitude_degree_in):
    """

    This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html


    Args:
        dLongitude_degree_in (float): The input longitude range from 0 to 360

    Returns:
        float: Longitude from -180 to 180
    """
    #a = int(dLongitude_degree_in /180)
    #dLongitude_out = dLongitude_degree_in - a*360.0

    dLongitude_out = ((dLongitude_degree_in + 180) % 360) - 180

    return dLongitude_out

#write a numpy version of the function

def convert_360_to_180_np(aLongitude_degree_in):

    """

    This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html


    Args:
        aLongitude_degree_in (numpy array): The input longitude range from 0 to 360

    Returns:
        numpy array: Longitude from -180 to 180
    """
    #aLongitude_degree_in = np.array(aLongitude_degree_in)
    #a = np.floor(aLongitude_degree_in /180)
    #aLongitude_out = aLongitude_degree_in - a*360.0

    aLongitude_degree_in = np.array(aLongitude_degree_in)
    aLongitude_out = ((aLongitude_degree_in + 180) % 360) - 180

    return aLongitude_out



def convert_180_to_360(dLongitude_degree_in):
    """

    This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html


    Args:
        dLongitude_degree_in (float): The input longitude range from -180 to 180

    Returns:
        float: Longitude from 0 to 360
    """
    dLongitude_out = (dLongitude_degree_in + 360.0) % 360.0

    return dLongitude_out