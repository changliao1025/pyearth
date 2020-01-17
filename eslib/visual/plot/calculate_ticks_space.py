

import math as Math
import numpy as np

def calculate_ticks_space (aData, nstep_in = 5):
    """
    Computes domain with given nstep_in encompassing series aData
    @ params
    aData    - Required - A list-like object of integers or floats
    nstep_in - Optional - Tick frequency
    """
    #xMax, xMin = Math.ceil(max(aData)), Math.floor(min(aData))
    xMax = Math.ceil( np.nanmax(aData) )
    xMin = Math.floor( np.nanmin(aData) )
    #dummy1 = (xMax % nstep_in) 
    #dummy2 = (xMin % nstep_in)
    #dMax =  xMax + abs( dummy1- nstep_in) + (nstep_in if (xMax % nstep_in != 0) else 0)
    #dMin = xMin - abs(dummy2)
    #dSpace = (dMax - dMin)/nstep_in
    dSpace = (xMax - xMin)/ (nstep_in-2)

    dMax = xMax + Math.ceil( dSpace )
    dMin =  Math.floor( xMin - dSpace) 
    return (dSpace, dMin, dMax)