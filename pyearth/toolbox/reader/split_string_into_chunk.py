def split_string_into_chunk(sString_in,iChunk_size_in=None): 
    """
    Split a string into pieces

    Args:
        sString_in (string): The input strin
        iChunk_size_in (int, optional): the number of chuck. Defaults to None.

    Yields:
        list: Subset of the string
    """
    iTotel_size = len(sString_in)
    if iChunk_size_in is not None:
        iChunk_size = iChunk_size_in
    else:
        iChunk_size = 10
        
    iChunk_count = iTotel_size//iChunk_size
    for pos in range(0, iTotel_size, iChunk_count):
        yield sString_in[pos:pos+iChunk_count]
    