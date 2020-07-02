import os,sys

sSystem_paths = os.environ['PATH'].split(os.pathsep)
sys.path.extend(sSystem_paths)
from eslib.system.define_global_variables import *

def prepare_parafly_python_command_file(iIndex_start, iIndex_end,\
    nThread, \
    sFilename_parafly, \
    sFilename_python):
    
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
if __name__ == '__main__':
    slurm_prepare_parafly_python_command()