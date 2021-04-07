import math as Math
import numpy as np

def calculate_ticks_space (aData, nstep_in = 5, iFlag_small_in = None):
    """
    Computes domain with given nstep_in encompassing series aData
    @ params
    aData    - Required - A list-like object of integers or floats
    nstep_in - Optional - Tick frequency
    """

    if iFlag_small_in is not None:
        iFlag_small = 1
    else:
        iFlag_small = 0
    if iFlag_small == 1:
        xMax = np.nanmax(aData) 
        xMin =  np.nanmin(aData) 
        dSpace = (xMax - xMin)/ (nstep_in-2)
        dMax = xMax * 1.2
        dMin =  xMin * 0.8

    else: 
        
        xMax = Math.ceil( np.nanmax(aData) )
        xMin = Math.floor( np.nanmin(aData) )
      
        dSpace = (xMax - xMin)/ (nstep_in-2)

        dMax = xMax + Math.ceil( dSpace )
        dMin =  Math.floor( xMin - dSpace) 
    return (dSpace, dMin, dMax)