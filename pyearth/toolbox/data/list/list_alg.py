import numpy as np


def list_max(pList):
    nlist = len(pList)
    dmax_max = np.nan  # Initialize with nan
    for i in range(nlist):
        sublist = pList[i]
        sublist = np.array(sublist)
        dmax = np.nanmax(sublist)
        if np.isnan(dmax_max):
            dmax_max = dmax  # First valid value (could still be nan)
        elif not np.isnan(dmax) and dmax > dmax_max:
            dmax_max = dmax

    return dmax_max


def list_min(pList):
    nlist = len(pList)
    dmin_min = np.nan  # Initialize with nan
    for i in range(nlist):
        sublist = pList[i]
        sublist = np.array(sublist)
        dmin = np.nanmin(sublist)
        if np.isnan(dmin_min):
            dmin_min = dmin  # First valid value (could still be nan)
        elif not np.isnan(dmin) and dmin < dmin_min:
            dmin_min = dmin

    return dmin_min
