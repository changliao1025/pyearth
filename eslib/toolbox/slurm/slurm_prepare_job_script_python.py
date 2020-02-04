import os
import subprocess

def slurm_prepare_job_script_python(iStart, iEnd, \
        sDirectory_job, \
        sFilename_checkpoint, \
        sFilename_job, \
        sFilename_python,\
        nNode_in = None, \
        nTask_in=None, \
        sPath_in = None):
    if nNode_in is not None:
        iNode = nNode_in            
    else:
        iNode = 1
    if nTask_in is not None:
        nTask = nTask_in            
    else:
        nTask = 40
    if sPath_in is not None:
        sPath = sPath_in            
    else:
        sPath = '.'

    sStart = "{:0d}".format(iStart  )
    sEnd =  "{:0d}".format(iEnd  )
    sNode =  "{:0d}".format(iNode  )
    sTask =   "{:0d}".format(nTask  )
    os.chdir(sPath)
    
    pFile =  open(sFilename_job,"w")  #write mode 
    sLine = '#!/bin/bash' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --account=e3sm' + '\n'
    pFile.write( sLine ) 

    sLine = '#SBATCH --begin=now+1minutes' + '\n'
    pFile.write( sLine ) 

    sLine = '#SBATCH --cpus-per-task=1 ' + '\n'
    pFile.write( sLine ) 

    sLine = '#SBATCH --dependency=singleton ' + '\n'
    pFile.write( sLine )
    sLine = '#SBATCH --error=stderr_%j' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --job-name=resubmit' + sStart + '_' + sEnd +   '# create a name for your job' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --mail-type=ALL' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --mail-user=chang.liao@pnnl.gov' + '\n'
    pFile.write( sLine ) 

    sLine = '#SBATCH --nodes=' + sNode + ' # node count' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --ntasks=' + sTask + ' # total number of tasks' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --output=stdout_%j.out' + '\n'
    pFile.write( sLine ) 

    sLine = '#SBATCH --partition=short' + '\n'  #can be improved here
    pFile.write( sLine ) 

    sLine = '#SBATCH --time=2:00:00      # total run time limit (HH:MM:SS)' + '\n'
    pFile.write( sLine ) 

    sLine = 'module purge' + '\n'
    pFile.write( sLine ) 
    sLine = 'module load anaconda3/2019.03' + '\n'
    pFile.write( sLine ) 
    sLine = 'source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh' + '\n'
    pFile.write( sLine ) 
    sLine = 'unset PYTHONHOME' + '\n'
    pFile.write( sLine ) 
    sLine = 'conda activate mpienv' + '\n'
    pFile.write( sLine ) 
    sLine = 'module load gcc/8.1.0' + '\n'
    pFile.write( sLine ) 
    sLine = 'module load openmpi/3.1.3' + '\n'
    pFile.write( sLine ) 
    sLine = 'export MPICC="$(which mpicc)"' + '\n'
    pFile.write( sLine ) 
    
    sLine = 'SLURM_SUBMIT_DIR=' + sDirectory_job + '\n'
    pFile.write( sLine ) 
    
    sLine = 'echo " Job " ' + '${SLURM_JOBID}' + ' is launched' + '\n'
    pFile.write( sLine ) 

    sLine = 'mpiexec --n ' +sTask + ' python ' + sFilename_python  + ' --sPath_job=' + sDirectory_job + '\n'
    pFile.write( sLine ) 

    sLine = 'conda deactivate' + '\n'
    pFile.write( sLine ) 

    #now prepare the resubmit part
    sLine = 'iIndex="$(sed -n' + '1p' + sFilename_checkpoint + ' | xargs)"' + '\n'
    pFile.write( sLine ) 
    sLine = 'iIndex="$(sed -n' + '2p' + sFilename_checkpoint + ' | xargs)"' + '\n'
    pFile.write( sLine ) 
    sLine = 'iIndex="$(sed -n' + '3p' + sFilename_checkpoint + ' | xargs)"' + '\n'
    pFile.write( sLine ) 

    sLine = 'if (($iIndex == 0));then' + '\n'
    pFile.write( sLine )
    sLine = '    echo "All tasks are finished"'      + '\n'
    pFile.write( sLine )
    sLine = 'else'      + '\n'
    pFile.write( sLine )
    sLine = '    echo "Job $SLURM_JOBID at $(date) will be re-submitted with new indices:, $sStart, $sEnd"'      + '\n'
    pFile.write( sLine )
    sLine = '    SLURM_SCRIPT_NAME="resubmit_${iIndex}.job"'      + '\n'
    pFile.write( sLine )

    sLine = '    sbatch $SLURM_SUBMIT_DIR/$SLURM_SCRIPT_NAME'      + '\n'
    pFile.write( sLine ) 
    sLine = 'fi'      + '\n'
    pFile.write( sLine ) 
    
    sLine = 'echo "Finished"'      + '\n'
    pFile.write( sLine ) 
    pFile.close() 
    
    return