
def slurm_update_checkpoint_file(iIndex_restart, \
    iIndex_start, \
    iIndex_end,\
    sFilename_checkpoint):
    """
    #write checkpoint file
    # save the checkpoint information
    """
    ofs = open(sFilename_checkpoint, "w")  #write mode 
    sLine = "{:0d}".format( iIndex_restart ) + '\n'
    ofs.write( sLine ) 
    sLine = "{:0d}".format( iIndex_start )  + '\n'
    ofs.write( sLine ) 
    sLine = "{:0d}".format( iIndex_end )  + '\n'
    ofs.write( sLine ) 
    ofs.close() 
    
        
    return