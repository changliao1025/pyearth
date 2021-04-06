# Locals



def check_crs(crs):
    """Checks if the crs represents a valid grid, projection or ESPG string.
    Examples
    --------
    >>> p = check_crs('+units=m +init=epsg:26915')
    >>> p.srs
    '+units=m +init=epsg:26915 '
    >>> p = check_crs('wrong')
    >>> p is None
    True
    Returns
    -------
    A valid crs if possible, otherwise None
    """

    try:
        crs = crs.salem.grid  # try xarray
    except:
        pass

    if isinstance(crs, string_types):
        # necessary for python 2
        crs = str(crs)

    if isinstance(crs, pyproj.Proj) or isinstance(crs, Grid):
        out = crs
    elif isinstance(crs, dict) or isinstance(crs, string_types):
        if isinstance(crs, string_types):
            # quick fix for https://github.com/pyproj4/pyproj/issues/345
            crs = crs.replace(' ', '').replace('+', ' +')
        try:
            out = pyproj.Proj(crs, preserve_units=True)
        except RuntimeError:
            try:
                out = pyproj.Proj(init=crs, preserve_units=True)
            except RuntimeError:
                out = None
    else:
        out = None
    return out