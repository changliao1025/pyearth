import numpy as np
def list_max(pList):
    nlist = len(pList)
    for i in range(nlist):
        sublist = pList[i]
        sublist = np.array(sublist)
        dmax = np.nanmax(sublist)
        if i == 0:
            dmax_max = dmax
        else:
            if dmax > dmax_max:
                dmax_max = dmax

    return dmax_max

def list_min(pList):
    nlist = len(pList)
    for i in range(nlist):
        sublist = pList[i]
        sublist = np.array(sublist)
        dmin = np.nanmin(sublist)
        if i == 0:
            dmin_min = dmin
        else:
            if dmin < dmin_min:
                dmin_min = dmin

    return dmin_min