import numpy as np

def search_neighbors(iIndex_in, aArray_in, iWindow_size_in = None):
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