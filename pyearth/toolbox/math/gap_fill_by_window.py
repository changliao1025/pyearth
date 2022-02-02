
import numpy as np
import statistics 
from numpy  import array
from pyearth.system.define_global_variables import *

def search_neighbors(iIndex_in, aArray_in, iWindow_size_in = None):
    """
    Gap filling a 2D array

    Args:
        iIndex_in (int): [description]
        aArray_in (numpy.array): Input data
        iWindow_size_in (int, optional): Search window. Defaults to None.

    Returns:
        numpy.array: Gap filled array
    """
    missing_value = np.nan
    if iWindow_size_in is not None:
        iWindow_size = iWindow_size_in           
    else:
        iWindow_size = 1

    iLength =  len(aArray_in)
    iStart = iWindow_size
    iEnd = iLength - iWindow_size
    iLength_out = iWindow_size * 2 + 1
    aArray_out = np.full( iLength_out , np.nan , dtype = float )    
    if iIndex_in >= iStart and iIndex_in < iEnd  :
        aArray_out=  aArray_in[ (iIndex_in - iWindow_size) : (iIndex_in + iWindow_size+1)   ]     
 
    return aArray_out


def gap_fill_by_window (aArray_in, iWindow_size_in = None):
    """
    Gap filling a 2D array

    Args:
        aArray_in (numpy.array): 2D array
        iWindow_size_in (int, optional): search window. Defaults to None.

    Returns:
        numpy.array: Gap filled 2D array
    """
    iLength = len(aArray_in)
    iWindow_size = iWindow_size_in
    aArray_out = aArray_in

    iStart = iWindow_size
    iEnd = iLength - iWindow_size
   
    for iIndex in range( iStart , iEnd  ):
        dummy = aArray_in[iIndex]
        if ( np.isnan(dummy) )  :
            aNeighbors = search_neighbors(iIndex, aArray_out, iWindow_size_in = iWindow_size)
            
            aNeighbors = aNeighbors.reshape( len(aNeighbors) )
            nan_index = np.where( np.isnan(aNeighbors) == True)

            nan_count  = len(nan_index[0])
             #only if half of the aNeighbors are qualified till we use this
             #gap fill
            if nan_count <= (iWindow_size +1):
               #aNeighbors=gap_fill(aNeighbors,/spline)

                good_index = np.where( np.isnan(aNeighbors) == False)

                aNeighbors_good = aNeighbors[ good_index ]
                good_count = len(aNeighbors_good)
                
                if good_count <  10:
                    aArray_out[iIndex] = statistics.median(aNeighbors_good)
                    pass
                   
                else:
                    aArray_out[iIndex] = statistics.median(aNeighbors_good)
                    pass
                
                

            else:
                
                pass
             #good elements count is not enough

            
        else:
   
            pass
     
        

    return aArray_out