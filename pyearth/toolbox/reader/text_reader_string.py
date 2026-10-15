"""Text reader to 2D numpy string array utility."""

import os
from typing import Optional, Union
import numpy as np


def text_reader_string(
    sFilename_in: str,
    ncolumn_in: Optional[int] = None,
    nrow_in: Optional[int] = None,
    cDelimiter_in: Optional[str] = None,
    iFlag_remove_quota: Optional[Union[int, bool]] = None,
    iSkipline_in: Optional[int] = None,
) -> np.ndarray:
    """
    Read a text-based file into a 2D NumPy array of strings.

    Parameters
    ----------
    sFilename_in : str
        Path to the text file.
    ncolumn_in : int, optional
        Expected number of columns. If None, inferred from the first parsed line.
    nrow_in : int, optional
        Maximum number of data rows to read. If None, reads all available lines.
    cDelimiter_in : str, optional
        Delimiter used to split tokens. If None, splits by arbitrary whitespace.
    iFlag_remove_quota : int or bool, optional
        If True (or 1), double-quotes are stripped from lines.
    iSkipline_in : int, optional
        Number of initial header lines to skip.

    Returns
    -------
    np.ndarray
        2D NumPy array containing the parsed string tokens.
    """
    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"File not found: {sFilename_in}")

    if ncolumn_in is not None and (not isinstance(ncolumn_in, (int, np.integer)) or int(ncolumn_in) <= 0):
        raise ValueError(f"Number of columns must be a positive integer, got {ncolumn_in}")

    target_cols = int(ncolumn_in) if ncolumn_in is not None else None
    target_rows = int(nrow_in) if nrow_in is not None else None
    strip_quotes = bool(iFlag_remove_quota)
    skip = int(iSkipline_in) if iSkipline_in is not None else 0

    rows = []
    with open(sFilename_in, "r", encoding="utf-8", errors="ignore") as f:
        for _ in range(skip):
            f.readline()

        for line in f:
            if target_rows is not None and len(rows) >= target_rows:
                break
            line_str = line.rstrip("\r\n")
            if strip_quotes:
                line_str = line_str.replace('"', "")

            if cDelimiter_in is not None:
                tokens = line_str.split(cDelimiter_in)
            else:
                tokens = line_str.split()

            if not tokens:
                continue

            if target_cols is not None:
                if len(tokens) < target_cols:
                    tokens = tokens + [""] * (target_cols - len(tokens))
                elif len(tokens) > target_cols:
                    tokens = tokens[:target_cols]

            rows.append(tokens)

    if not rows:
        col_count = target_cols if target_cols is not None else 0
        return np.empty((0, col_count), dtype=object)

    num_cols = target_cols if target_cols is not None else max(len(r) for r in rows)
    padded = [r + [""] * (num_cols - len(r)) for r in rows]
    return np.array(padded, dtype=object)
