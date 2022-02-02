import os

def prepare_parafly_slurm_job_script(sBasename_job,  sBasename_parafly,   sDirectory_job, sEmail, iWalltime_in = None, nNode_in = None, nThread_in=None, sJob_name_in =None,  sPython_env_in =None, sQueue_in=None):
    """
    prepare a job script for parafly

    Args:
        sBasename_job (string): The base job name
        sBasename_parafly (string): The parafly script filename
        sDirectory_job (string): The path where parafly runs the script
        sEmail (string): Email address for job notification
        iWalltime_in (int, optional): Walltime in hour. Defaults to None.
        nNode_in (int, optional): NUmber of node. Defaults to None.
        nThread_in (int, optional): Number of node. Defaults to None.
        sJob_name_in (string, optional): Job name. Defaults to None.
        sPython_env_in (string, optional): Conda python environment. Defaults to None.
        sQueue_in (string, optional): HPC queue or partition. Defaults to None.
    """
    if iWalltime_in is not None:
        iWalltime = iWalltime_in            
    else:
        iWalltime = 2
    if nNode_in is not None:
        iNode = nNode_in            
    else:
        iNode = 1
    if nThread_in is not None:
        nThread = nThread_in            
    else:
        nThread = 40
   
    if sJob_name_in is not None:
        sJob_name = sJob_name_in            
    else:
        sJob_name = 'parafly'
    if sPython_env_in is not None:
        sPython_env = sPython_env_in            
    else:
        sPython_env = 'base'
        
    if sQueue_in is not None:
        sQueue = sQueue_in            
    else:
        sQueue = 'short'
        
    sWalltime ="{:0d}".format(iWalltime  )
    sNode =  "{:0d}".format(iNode  )
    sThread =   "{:0d}".format(nThread  )
    
    os.chdir(sDirectory_job)
    
    ofs =  open(sBasename_job,"w")  #write mode 
    sLine = '#!/bin/bash' + '\n'
    ofs.write( sLine ) 
    sLine = '#SBATCH --account=e3sm' + '\n'
    ofs.write( sLine ) 

    #sLine = '#SBATCH --begin=now+1minutes' + '\n'
    #ofs.write( sLine ) 

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
    sLine = '#SBATCH --mail-user=' + sEmail + '\n'
    ofs.write( sLine ) 

    sLine = '#SBATCH --nodes=' + sNode + ' # node count' + '\n'
    ofs.write( sLine ) 
    sLine = '#SBATCH --ntasks=' + sThread + ' # total number of tasks' + '\n'
    ofs.write( sLine ) 
    sLine = '#SBATCH --output=stdout_%j.out' + '\n'
    ofs.write( sLine ) 

    sLine = '#SBATCH --partition=' + sQueue + '\n'  #can be improved here
    ofs.write( sLine ) 
    sLine = '#SBATCH --time=' + sWalltime +':00:00   # total run time limit (HH:MM:SS)' + '\n'
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
    sLine = 'conda activate ' + sPython_env + '\n'
    ofs.write( sLine )   

    sLine = 'ParaFly -c ' + sBasename_parafly + ' -CPU ' + sThread + ' -failed_cmds rerun.txt' + '\n'
    ofs.write( sLine ) 
    
    sLine = 'echo " Job " ' + '${SLURM_JOBID}' + ' is launched' + '\n'
    ofs.write( sLine ) 

    sLine = 'conda deactivate' + '\n'
    ofs.write( sLine ) 
    
    sLine = 'echo "Finished"'      + '\n'
    ofs.write( sLine ) 
    ofs.close() 
    
    return