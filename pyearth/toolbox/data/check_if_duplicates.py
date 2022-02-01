def check_if_duplicates(aList_in):
    """
    Check if given list contains any duplicates
    Args:
        aList_in (list): List of native Python data type

    Returns:
        int: 1-no duplicate, 0-has duplicate
    """
    
    iFlag_unique = 1
    for elem in aList_in:
        if aList_in.count(elem) > 1:
            iFlag_unique = 0
            break
        else:
            pass
    
    return iFlag_unique