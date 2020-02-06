import os
import subprocess

def slurm_prepare_job_script_parafly(sDirectory_job, \
        sBasename_job, \
        sBasename_parafly, \
        sJob_name, \
        iWalltime_in = None, \
        nNode_in = None, \
        nTask_in=None, \
        sQueue_in=None):
    if nNode_in is not None:
        iNode = nNode_in            
    else:
        iNode = 1
    if nTask_in is not None:
        nTask = nTask_in            
    else:
        nTask = 40
    if nNode_in is not None:
        iNode = nNode_in            
    else:
        iNode = 1
    if sQueue_in is not None:
        sQueue = sQueue_in            
    else:
        sQueue = 'short'
    if iWalltime_in is not None:
        iWalltime = iWalltime_in            
    else:
        iWalltime = 2
   
    sNode =  "{:0d}".format(iNode  )
    sTask =   "{:0d}".format(nTask  )
    sWalltime ="{:0d}".format(iWalltime  )
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
    sLine = 'module load parafly/2013' + '\n'
    ofs.write( sLine ) 
    sLine = 'module load anaconda3/2019.03' + '\n'
    ofs.write( sLine ) 
    sLine = 'source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh' + '\n'
    ofs.write( sLine ) 
    sLine = 'unset PYTHONHOME' + '\n'
    ofs.write( sLine ) 
    sLine = 'conda activate mpienv' + '\n'
    ofs.write( sLine )   
    sFilename_parafly = sBasename_parafly

    sLine = 'ParaFly -c ' + sFilename_parafly + ' -CPU -failed_cmds rerun.txt' + '\n'
    ofs.write( sLine ) 
    
    
    

    sLine = 'echo " Job " ' + '${SLURM_JOBID}' + ' is launched' + '\n'
    ofs.write( sLine ) 

    sLine = 'conda deactivate' + '\n'
    ofs.write( sLine ) 
    
    
    sLine = 'echo "Finished"'      + '\n'
    ofs.write( sLine ) 
    ofs.close() 
    
    return