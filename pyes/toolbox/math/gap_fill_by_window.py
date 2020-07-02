
import statistics 
import math
import sys
import numpy as np
from numpy  import array

#use global method
sSystem_paths = os.environ['PATH'].split(os.pathsep)
sys.path.extend(sSystem_paths)
from eslib.system import define_global_variables
from eslib.system.define_global_variables import *

#import the eslib library
sPath_library_python = sWorkspace_code +  slash + 'python' + slash + 'library' + slash + 'eslib_python'
sys.path.append(sPath_library_python)
from eslib.toolbox.math.search_neighbors import search_neighbors

def gap_fill_by_window (aArray_in, iWindow_size_in = None):
    """aArray_in,
     
    iWindow_size_in =  None
    """
    iLength = len(aArray_in)

    iWindow_size = iWindow_size_in

    aArray_out = aArray_in

    iStart = iWindow_size
    iEnd = iLength - iWindow_size
    #print(  np.nansum(aArray_in) )
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
                    #print(aArray_out[iIndex])
                else:
                    aArray_out[iIndex] = statistics.median(aNeighbors_good)
                    #print(aArray_out[iIndex])
                #print(  np.nansum(aArray_out) )

            else:
                #print(nan_count)
                pass
             #good elements count is not enough

            #print(  np.nansum(aArray_out) )
        else:
            #print(dummy)
            pass
     
        

    #print(  np.nansum(aArray_out) )
    return aArray_out