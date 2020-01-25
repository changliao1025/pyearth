import os
import subprocess
def submit_python_job(iStart, iEnd, sFilename_python, sFilename_job, iNode_in = None, nTask_in=None, sPath_in = None):
    if iNode_in is not None:
        iNode = iNode_in            
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
    sFilename_job = sPath + slash + 'submit.job'
    #print(sFilename_job, time.time())
    pFile =  open(sFilename_job,"w")  #write mode 
    sLine = '#!/bin/bash' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --account=e3sm' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --output=stdout_mpi' + sStart + '_' + sEnd  + '.out' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --error=stderr_mpi' + sStart + '_' + sEnd  + 'err' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --mail-type=ALL' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --mail-user=chang.liao@pnnl.gov' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --partition=short' + '\n'  #can be improved here
    pFile.write( sLine ) 
    sLine = '#SBATCH --job-name=resubmit' + sStart + '_' + sEnd +   '# create a name for your job' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --nodes=' + sNode + ' # node count' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --ntasks=' + sTask + ' # total number of tasks' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --cpus-per-task=1        # cpu-cores per task' + '\n'
    pFile.write( sLine ) 
    sLine = '#SBATCH --time=2:00:00          # total run time limit (HH:MM:SS)' + '\n'
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
    #sLine = 'JOB_DIRECTORY=/people/liao313/workspace/python/lab/lab_python/slurm/' + '\n'
    #pFile.write( sLine ) 
    #sLine = 'cd $JOB_DIRECTORY' + '\n'
    #pFile.write( sLine ) 
    #we make an assumption here that the python script must accept at least three parameters, start index, end index and the job path
    sLine = 'mpiexec --n ' +sTask + ' python ' + sFilename_python + ' --iStart_index=' + sStart +' --iEnd_index=' + sEnd + ' --sPath=' + sFilename_job + '\n'
    pFile.write( sLine ) 
    sLine = 'conda deactivate' + '\n'
    pFile.write( sLine ) 
    sLine = 'echo "Finished"'      + '\n'
    pFile.write( sLine ) 
    pFile.close() 

    sLine = 'sbatch ' + sFilename_job
    subprocess.call(sLine, shell=True)
    return