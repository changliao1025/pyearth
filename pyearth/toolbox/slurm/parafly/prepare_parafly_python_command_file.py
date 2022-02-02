import os,sys

from pyearth.system.define_global_variables import *

def prepare_parafly_python_command_file(iIndex_start,     iIndex_end,    nThread,     sFilename_parafly,     sFilename_python):
    """
    Prepare a parafly file which runs python script

    Args:
        iIndex_start (int): The start index
        iIndex_end (int): The end index
        nThread (int): Number of thread
        sFilename_parafly (string): The python binary
        sFilename_python (string): The parafly filename
    """
    
    ofs =  open(sFilename_parafly,"w")  #write mode 
    nTask = iIndex_end - iIndex_start + 1
    nChunkPerThread = nTask // nThread 
    
    for iRank in range(nThread):
        if iRank == 0:
            iStart = (nThread-1) * nChunkPerThread + iIndex_start
            iEnd = iIndex_end      
        else:
            iStart = (iRank-1) * nChunkPerThread + iIndex_start
            iEnd = (iRank) * nChunkPerThread + iIndex_start - 1
            
        sStart = "{:0d}".format( iStart ) 
        sEnd = "{:0d}".format( iEnd ) 
        sLine = 'python ' + sFilename_python \
            + ' --iIndex_start ' + sStart \
            + ' --iIndex_end '  + sEnd + '\n'
       
        ofs.write(sLine) 
    
    ofs.close()
    
    return
