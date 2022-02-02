import os, sys
from pyearth.system.define_global_variables import *

def slurm_prepare_job_script_python(iIndex_start, iIndex_end,      sBasename_job,     sDirectory_job,     sFilename_python,       iWalltime_in = None,     nNode_in = None,     nThread_in = None,   sAccount =None,    sEmail=None,      sQueue_in = None):
    """
    Prepare a slurm job script for python without checkpoint

    Args:
        iIndex_start ([type]): The start index
        iIndex_end ([type]): The end index
        sBasename_job ([type]): The basename of job
        sDirectory_job ([type]): [description]
        sFilename_python ([type]): [description]
        iWalltime_in ([type], optional): [description]. Defaults to None.
        nNode_in ([type], optional): [description]. Defaults to None.
        nThread_in ([type], optional): [description]. Defaults to None.
        sAccount ([type], optional): [description]. Defaults to None.
        sEmail ([type], optional): [description]. Defaults to None.
        sQueue_in ([type], optional): [description]. Defaults to None.
    """
    if nNode_in is not None:
        iNode = nNode_in
    else:
        iNode = 1
    if nThread_in is not None:
        nThread = nThread_in
    else:
        nThread = 24

    if sQueue_in is not None:
        sQueue = sQueue_in
    else:
        sQueue = 'short'
    if iWalltime_in is not None:
        iWalltime = iWalltime_in
    else:
        iWalltime = 2


    nTask = iIndex_end - iIndex_start + 1
    nChunkPerThread = nTask // nThread

    sNode =  "{:0d}".format(iNode  )
    sTask =   "{:0d}".format(nTask  )
    sWalltime ="{:0d}".format(iWalltime  )
    os.chdir(sDirectory_job)

    for iRank in range(1, nThread+1):
        if iRank == 1:
            iStart = (nThread-1) * nChunkPerThread + iIndex_start
            iEnd = iIndex_end
        else:
            iStart = (iRank-2) * nChunkPerThread + iIndex_start
            iEnd = (iRank-1) * nChunkPerThread + iIndex_start - 1
            sStart = "{:0d}".format(iStart  )
            sEnd =  "{:0d}".format(iEnd  )

        sFilename_job = sDirectory_job + slash + sBasename_job + "{:03d}".format(iRank) + '.job'
        pFile =  open(sFilename_job,"w")  #write mode
        sLine = '#!/bin/bash' + '\n'
        pFile.write( sLine )

        if sAccount is not None:
            sLine = '#SBATCH --account=' + sAccount  + '\n'
            pFile.write( sLine )

        sLine = '#SBATCH --begin=now+1minutes' + '\n'
        #pFile.write( sLine )

        sLine = '#SBATCH --cpus-per-task=1 ' + '\n'
        pFile.write( sLine )

        sLine = '#SBATCH --dependency=singleton ' + '\n'
        #pFile.write( sLine )
        sLine = '#SBATCH --error=stderr_%j.err' + '\n'
        pFile.write( sLine )
        sLine = '#SBATCH --job-name=' \
            + sBasename_job + "{:03d}".format(iRank) + '  # create a name for your job' + '\n'
        pFile.write( sLine )
        sLine = '#SBATCH --mail-type=ALL' + '\n'
        pFile.write( sLine )

        if sEmail is not None:
            sLine = '#SBATCH --mail-user=' + sEmail + '\n'
            pFile.write( sLine )

        sLine = '#SBATCH --nodes=' + sNode + ' # node count' + '\n'
        pFile.write( sLine )
        sLine = '#SBATCH --ntasks=1'  + ' # total number of tasks' + '\n'
        pFile.write( sLine )
        sLine = '#SBATCH --output=stdout_%j.out' + '\n'
        pFile.write( sLine )

        sLine = '#SBATCH --partition=' + sQueue + '\n'  #can be improved here
        pFile.write( sLine )

        sLine = '#SBATCH --time=' + sWalltime +':00:00      # total run time limit (HH:MM:SS)' + '\n'
        pFile.write( sLine )

        sLine = 'module purge' + '\n'
        pFile.write( sLine )
        sLine = 'module load anaconda3/2019.03' + '\n'
        pFile.write( sLine )
        sLine = 'source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh' + '\n'
        pFile.write( sLine )
        sLine = 'unset PYTHONHOME' + '\n'
        pFile.write( sLine )
        sLine = 'conda activate gdalenv' + '\n'
        pFile.write( sLine )

        sLine = 'python ' + sFilename_python \
            + ' --iIndex_start ' + sStart \
            + ' --iIndex_end '  + sEnd + '\n'
        pFile.write( sLine )

        sLine = 'echo " Job " ' + '${SLURM_JOBID}' + ' is launched' + '\n'
        pFile.write( sLine )

        sLine = 'conda deactivate' + '\n'
        pFile.write( sLine )
        sLine = 'echo "Finished"'      + '\n'
        pFile.write( sLine )
        pFile.close()

    return
