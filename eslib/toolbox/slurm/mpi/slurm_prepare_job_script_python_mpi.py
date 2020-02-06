import os

#with mpi, we assume slave are uniform distributed, and the checkpoint is optional
def slurm_prepare_job_script_python_mpi(iIndex_start, iIndex_end, \
        
        sBasename_job, \
        sBasename_python,\
        sDirectory_job, \
        sDirectory_python,\
        iFlag_checkpoint_in=None, \
        iWalltime_in=None,\
        nNode_in = None, \
        nTask_in=None, \
        sBasename_checkpoint_in=None, \
        sJob_name_in= None, \    
        sQueue_in=None):
    if iFlag_checkpoint_in is not None:
        iFlag_checkpoint = 1            
    else:
        iFlag_checkpoint = 0
    if iWalltime_in is not None:
        iWalltime = iWalltime_in            
    else:
        iWalltime = 2
    if nNode_in is not None:
        iNode = nNode_in            
    else:
        iNode = 1
    if nTask_in is not None:
        nTask = nTask_in            
    else:
        nTask = 40
    if sBasename_checkpoint_in is not None:
        sBasename_checkpoint = sBasename_checkpoint_in            
    else:
        sBasename_checkpoint = 'checkpoint.txt'
    if sJob_name_in is not None:
        sJob_name = sJob_name_in            
    else:
        sJob_name = 'auto_resubmit'
    if sQueue_in is not None:
        sQueue = sQueue_in            
    else:
        sQueue = 'short'
    
    sStart = "{:0d}".format(iIndex_start  )
    sEnd =  "{:0d}".format(iIndex_end  )
    sWalltime ="{:0d}".format(iWalltime  )
    sNode =  "{:0d}".format(iNode  )
    sTask =   "{:0d}".format(nTask  )
    
    os.chdir(sDirectory_job)
    
    ofs =  open(sBasename_job,"w")  #write mode 
    sLine = '#!/bin/bash' + '\n'
    ofs.write( sLine ) 
    sLine = '#SBATCH --account=e3sm' + '\n'
    ofs.write( sLine ) 

    sLine = '#SBATCH --begin=now+1minutes' + '\n'
    ofs.write( sLine ) 

    sLine = '#SBATCH --cpus-per-task=1 ' + '\n'
    ofs.write( sLine ) 

    sLine = '#SBATCH --dependency=singleton ' + '\n'
    ofs.write( sLine )
    sLine = '#SBATCH --error=stderr_%j.err' + '\n'
    ofs.write( sLine ) 
    sLine = '#SBATCH --job-name=' + sJob_name + '  # create a name for your job' + '\n'
    ofs.write( sLine ) 
    sLine = '#SBATCH --mail-type=ALL' + '\n'
    ofs.write( sLine ) 
    sLine = '#SBATCH --mail-user=chang.liao@pnnl.gov' + '\n'
    ofs.write( sLine ) 

    sLine = '#SBATCH --nodes=' + sNode + ' # node count' + '\n'
    ofs.write( sLine ) 
    sLine = '#SBATCH --ntasks=' + sTask + ' # total number of tasks' + '\n'
    ofs.write( sLine ) 
    sLine = '#SBATCH --output=stdout_%j.out' + '\n'
    ofs.write( sLine ) 

    sLine = '#SBATCH --partition=' + sQueue + '\n'  #can be improved here
    ofs.write( sLine ) 

    sLine = '#SBATCH --time=' + sWalltime +':00:00      # total run time limit (HH:MM:SS)' + '\n'
    ofs.write( sLine ) 

    sLine = 'module purge' + '\n'
    ofs.write( sLine ) 
    sLine = 'module load anaconda3/2019.03' + '\n'
    ofs.write( sLine ) 
    sLine = 'source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh' + '\n'
    ofs.write( sLine ) 
    sLine = 'unset PYTHONHOME' + '\n'
    ofs.write( sLine ) 
    sLine = 'conda activate mpienv' + '\n'
    ofs.write( sLine ) 
    sLine = 'module load gcc/8.1.0' + '\n'
    ofs.write( sLine ) 
    sLine = 'module load openmpi/3.1.3' + '\n'
    ofs.write( sLine ) 
    sLine = 'export MPICC="$(which mpicc)"' + '\n'
    ofs.write( sLine ) 
    
    sLine = 'sDirectory_job=' + sDirectory_job + '\n'
    ofs.write( sLine ) 
    sLine = 'sDirectory_python=' + sDirectory_python + '\n'
    ofs.write( sLine ) 
    sLine = 'sBasename_python=' + sBasename_python + '\n'
    ofs.write( sLine ) 
    sLine = 'sFilename_python=' + '$sDirectory_python' + "/" + '$sBasename_python' + '\n'
    ofs.write( sLine ) 

    sLine = 'echo " Job " ' + '${SLURM_JOBID}' + ' is launched' + '\n'
    ofs.write( sLine ) 

    sLine = 'mpiexec --n ' + sTask + ' python  $sFilename_python '\
         + ' --sDirectory_job=$sDirectory_job' + '\n'
    ofs.write( sLine ) 

    sLine = 'conda deactivate' + '\n'
    ofs.write( sLine ) 
    
    #now prepare the resubmit part
    sLine = 'iIndex="$(sed -n ' + ' 1p ' + sBasename_checkpoint + ' | xargs)"' + '\n'
    ofs.write( sLine ) 
    sLine = 'iIndex="$(sed -n ' + ' 2p ' + sBasename_checkpoint + ' | xargs)"' + '\n'
    ofs.write( sLine ) 
    sLine = 'iIndex="$(sed -n ' + ' 3p ' + sBasename_checkpoint + ' | xargs)"' + '\n'
    ofs.write( sLine ) 

    sLine = 'if (($iIndex == 0));then' + '\n'
    ofs.write( sLine )
    sLine = '    echo "All tasks are finished"'      + '\n'
    ofs.write( sLine )
    sLine = 'else'      + '\n'
    ofs.write( sLine )
    sLine = '    echo "Job $SLURM_JOBID at $(date) will be re-submitted with new indices:, $sStart, $sEnd"'      + '\n'
    ofs.write( sLine )
    sLine = '    sBasename_job="resubmit_${iIndex}.job"'      + '\n'
    ofs.write( sLine )

    sLine = '    sbatch $sDirectory_job/$sBasename_job'      + '\n'
    ofs.write( sLine ) 
    sLine = 'fi'      + '\n'
    ofs.write( sLine ) 
    
    sLine = 'echo "Finished"'      + '\n'
    ofs.write( sLine ) 
    ofs.close() 
    
    return