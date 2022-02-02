def line_count(sFilename_in):
    """
    Count the line count in a text file

    Args:
        sFilename_in (string): The filename

    Returns:
        int: The total line count
    """
    ifs=open(sFilename_in, 'rb') 
    i=0 
    sLine0=(ifs.readline())#.rstrip()
    sLine= sLine0.decode("utf-8", 'ignore')
    while len(sLine) > 0:
        i = i+1
        sLine0=(ifs.readline())
        sLine= sLine0.decode("utf-8", 'ignore')
        
        
    return i