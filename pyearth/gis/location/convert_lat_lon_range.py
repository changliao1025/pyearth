def convert_360_to_180(dLongitude_in):
    """
   
    This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html
    

    Args:
        dLongitude_in (float): The input longitude range from 0 to 360

    Returns:
        float: Longitude from -180 to 180
    """
    a = int(dLongitude_in /180)
    dLongitude_out = dLongitude_in - a*360.0

    return dLongitude_out



def convert_180_to_360(dLongitude_in):
    """
   
    This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html
    

    Args:
        dLongitude_in (float): The input longitude range from -180 to 180

    Returns:
        float: Longitude from 0 to 360
    """
    dLongitude_out = (dLongitude_in + 360.0) % 360.0

    return dLongitude_out