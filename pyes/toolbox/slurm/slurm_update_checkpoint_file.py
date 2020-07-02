#save the checkpoint information
def slurm_update_checkpoint_file(iIndex_restart, \
    iIndex_start, \
    iIndex_end,\
    sFilename_checkpoint):
    #write checkpoint file
    ofs = open(sFilename_checkpoint, "w")  #write mode 
    sLine = "{:0d}".format( iIndex_restart ) + '\n'
    ofs.write( sLine ) 
    sLine = "{:0d}".format( iIndex_start )  + '\n'
    ofs.write( sLine ) 
    sLine = "{:0d}".format( iIndex_end )  + '\n'
    ofs.write( sLine ) 
    ofs.close() 
    #if(iIndex_restart == 0):
    #    print('Checkpoint file is resetted to: ', iIndex_restart,  iIndex_start, iIndex_end)
    #else:
    #    print('Checkpoint file is updated as: ', iIndex_restart,  iIndex_start, iIndex_end)
        
    return