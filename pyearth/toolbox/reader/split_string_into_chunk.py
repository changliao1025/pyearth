"""Split string into chunks utility."""

from typing import Generator, Optional


def split_string_into_chunk(
    sString_in: str, iChunk_size_in: Optional[int] = None
) -> Generator[str, None, None]:
    """
    Split a string into pieces of given chunk length.

    Parameters
    ----------
    sString_in : str
        The input string.
    iChunk_size_in : int, optional
        Length of each chunk (default is 10). Must be positive.

    Yields
    ------
    str
        Chunks of the string.
    """
    if not isinstance(sString_in, str):
        raise TypeError(f"Input must be str, got {type(sString_in).__name__}")

    chunk_size = 10 if iChunk_size_in is None else int(iChunk_size_in)
    if chunk_size <= 0:
        raise ValueError(f"Chunk size must be positive, got {chunk_size}")

    for pos in range(0, len(sString_in), chunk_size):
        yield sString_in[pos : pos + chunk_size]
