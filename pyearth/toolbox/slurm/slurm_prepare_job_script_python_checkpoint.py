import os

#normal job, no checkpoint
def slurm_prepare_job_script_python(iStart,  iEnd,  sBasename_checkpoint, sBasename_job, sBasename_python, sDirectory_job, sDirectory_python, sJob_name, iFlag_resubmit =None,  iWalltime_in = None, nNode_in = None, nTask_in=None, sAccount = None, sEmail=None, sQueue_in=None):
    """
    Prepare a slurm job script for python with checkpoint

    Args:
        iStart (int): Checkpoint start index
        iEnd (int): Checkpoint end index
        sBasename_checkpoint (string): Basename 
        sBasename_job (string): Job name
        sBasename_python (string): Python path
        sDirectory_job (string): The path where contains the bash script
        sDirectory_python (string): The folder of Python
        sJob_name (string): Job name
        iFlag_resubmit (int, optional): Flag to indicate this is a resubmit. Defaults to None.
        iWalltime_in (int, optional): The walltime in hour. Defaults to None.
        nNode_in (int, optional): Number of node. Defaults to None.
        nTask_in (int, optional): Number of task. Defaults to None.
        sAccount (string, optional): Cluster account. Defaults to None.
        sEmail (string, optional): Email address. Defaults to None.
        sQueue_in (string, optional): Cluster partition. Defaults to None.
    """

    if nNode_in is not None:
        iNode = nNode_in
    else:
        iNode = 1

    if nTask_in is not None:
        nTask = nTask_in
    else:
        nTask = 10



    if sQueue_in is not None:
        sQueue = sQueue_in
    else:
        sQueue = 'short'

    if iWalltime_in is not None:
        iWalltime = iWalltime_in
    else:
        iWalltime = 2

    sStart = "{:0d}".format(iStart  )
    sEnd =  "{:0d}".format(iEnd  )
    sNode =  "{:0d}".format(iNode  )
    sTask =   "{:0d}".format(nTask  )
    sWalltime ="{:0d}".format(iWalltime  )
    os.chdir(sDirectory_job)

    pFile =  open(sBasename_job,"w")  #write mode
    sLine = '#!/bin/bash' + '\n'
    pFile.write( sLine )

    if sAccount is not None:
        sLine = '#SBATCH --account=' + sAccount  + '\n'
        pFile.write( sLine )

    sLine = '#SBATCH --begin=now+1minutes' + '\n'
    pFile.write( sLine )

    sLine = '#SBATCH --cpus-per-task=1 ' + '\n'
    pFile.write( sLine )

    sLine = '#SBATCH --dependency=singleton ' + '\n'
    pFile.write( sLine )
    sLine = '#SBATCH --error=stderr_%j.err' + '\n'
    pFile.write( sLine )

    if sJob_name is not None:
        sLine = '#SBATCH --job-name=' + sJob_name + '  # create a name for your job' + '\n'
        pFile.write( sLine )

    sLine = '#SBATCH --mail-type=ALL' + '\n'
    pFile.write( sLine )

    if sEmail is not None:
        sLine = '#SBATCH --mail-user=' + sEmail + '\n'
        pFile.write( sLine )

    sLine = '#SBATCH --nodes=' + sNode + ' # node count' + '\n'
    pFile.write( sLine )
    sLine = '#SBATCH --ntasks=' + sTask + ' # total number of tasks' + '\n'
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
    #only need to manuually source if the cluster requires
    #sLine = 'source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh' + '\n'
    #pFile.write( sLine )
    sLine = 'unset PYTHONHOME' + '\n'
    pFile.write( sLine )
    sLine = 'conda activate mpienv' + '\n'
    pFile.write( sLine )
    sLine = 'sDirectory_job=' + sDirectory_job + '\n'
    pFile.write( sLine )
    sLine = 'sDirectory_python=' + sDirectory_python + '\n'
    pFile.write( sLine )
    sLine = 'sBasename_python=' + sBasename_python + '\n'
    pFile.write( sLine )
    sLine = 'sFilename_python=' + '$sDirectory_python' + "/" + '$sBasename_python' + '\n'
    pFile.write( sLine )

    sLine = 'echo " Job " ' + '${SLURM_JOBID}' + ' is launched' + '\n'
    pFile.write( sLine )

    sLine = 'conda deactivate' + '\n'
    pFile.write( sLine )

    if iFlag_resubmit ==1: #advance job for checkpoint
        #now prepare the resubmit part
        sLine = 'iIndex="$(sed -n ' + ' 1p ' + sBasename_checkpoint + ' | xargs)"' + '\n'
        pFile.write( sLine )
        sLine = 'iIndex="$(sed -n ' + ' 2p ' + sBasename_checkpoint + ' | xargs)"' + '\n'
        pFile.write( sLine )
        sLine = 'iIndex="$(sed -n ' + ' 3p ' + sBasename_checkpoint + ' | xargs)"' + '\n'
        pFile.write( sLine )

        sLine = 'if (($iIndex == 0));then' + '\n'
        pFile.write( sLine )
        sLine = '    echo "All tasks are finished"'      + '\n'
        pFile.write( sLine )
        sLine = 'else'      + '\n'
        pFile.write( sLine )
        sLine = '    echo "Job $SLURM_JOBID at $(date) will be re-submitted with new indices:, $sStart, $sEnd"'        + '\n'
        pFile.write( sLine )
        sLine = '    sBasename_job="resubmit_${iIndex}.job"'      + '\n'
        pFile.write( sLine )

        sLine = '    sbatch $sDirectory_job/$sBasename_job'      + '\n'
        pFile.write( sLine )
        sLine = 'fi'      + '\n'
        pFile.write( sLine )

    else:
        pass

    sLine = 'echo "Finished"'      + '\n'
    pFile.write( sLine )
    pFile.close()

    return
