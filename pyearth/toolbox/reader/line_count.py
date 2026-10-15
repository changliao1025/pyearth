"""Line count utility."""

import os


def line_count(sFilename_in: str) -> int:
    """
    Count the total number of lines in a text file.

    Parameters
    ----------
    sFilename_in : str
        Path to the text file.

    Returns
    -------
    int
        The total line count.
    """
    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"File not found: {sFilename_in}")

    count = 0
    with open(sFilename_in, "r", encoding="utf-8", errors="ignore") as ifs:
        for _ in ifs:
            count += 1
    return count
