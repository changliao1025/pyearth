
def slurm_update_checkpoint_file(iIndex_restart,     iIndex_start,     iIndex_end,    sFilename_checkpoint):
    """
    Write and save the checkpoint information

    Args:
        iIndex_restart (int): [description]
        iIndex_start (int): [description]
        iIndex_end (int): [description]
        sFilename_checkpoint (string): [description]
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