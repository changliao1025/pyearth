import os
import numpy as np
def text_reader_string( sFilename_in,\
     ncolumn_in = None, \
     nrow_in = None, \
     delimiter_in = None, \
     remove_quota = None, \
     skipline_in =  None ):
    """
    sFilename_in,\
    ncolumn_in = None, \
    nrow_in = None, \
    delimiter_in = None, \
    skipline_in =  None
    """
    #print(sFilename_in)
    #print(os.path.isfile(sFilename_in))
    aData_out = -1
    if os.path.isfile(sFilename_in):
        #print('file exist')
        #check the parameter
        if ncolumn_in is not None:
            iFlag_column = 1
            ncolumn_out = ncolumn_in
        else:
            iFlag_column = 0
        if nrow_in is not None :
            #there is nothing
            nrow_out = nrow_in
            #print(' ')
        else :
            #get the total line number
            ifs = open(sFilename_in, "r")
            nrow_out = len(ifs.readlines())
            ifs.close()
        sLine=' '
        ifs = open(sFilename_in, "r")

        if skipline_in is not None:
            nrow_out =  nrow_out - skipline_in
            for i in range(skipline_in):
                ifs.readline()
        else:
            pass
        if remove_quota is not None:
            iFlag_remove_quota = 1
        else:
            iFlag_remove_quota = 0
        #get delimiter
        if delimiter_in is not None:
            iFlag_delimiter = 1
        else:
            iFlag_delimiter = 0
            #delimiter_in = ' '
        sLine = (ifs.readline()).rstrip()
        
        if iFlag_delimiter == 1:
            if iFlag_column == 1:
                pass
            else :
                dummy = sLine.split(delimiter_in)
                
                ncolumn_out = len(dummy)
            #check ncolumn_in count
            if ncolumn_out < 1 :
                print( 'The file has no data!' )
                        
            #aData_out = [[0 for x in range(ncolumn_out)] for y in range(nrow_out)]   
            aData_out = np.full( (nrow_out, ncolumn_out), '  ' , dtype=object )          
            dummy = sLine.split(delimiter_in)
           
            aData_out[0] =  dummy

            for iRow in range(1, nrow_out):
                sLine=(ifs.readline()).rstrip()
                if iFlag_remove_quota ==1:
                    sLine.replace('"','')
                else:
                    pass    
                aData_out[iRow] = sLine.split(delimiter_in)
        else :
            if iFlag_column == 1:
                pass
            else :
                dummy = sLine.split()
                ncolumn_out = len(dummy)
            #check ncolumn_in count
            if ncolumn_out < 1 :                
                return aData_out                      
            #aData_out = [[0 for x in range(ncolumn_out)] for y in range(nrow_out)] 
            aData_out = np.full( (nrow_out, ncolumn_out), '  ', dtype=object )
            if iFlag_remove_quota ==1:
                sLine=sLine.replace('"','')
            else:
                pass 

            dummy = sLine.split()
           
            aData_out[0] =  dummy

            for iRow in range(1,nrow_out):
                dummy1 = ifs.readline()
                #print('New line is: ' + dummy1)
                dummy2 = dummy1.rstrip()
                if iFlag_remove_quota ==1:
                    dummy2=dummy2.replace('"','')
                else:
                    pass
                dummy3 = dummy2.split()
                #print(dummy3)
                aData_out[iRow] = dummy3
        #print(aData_out)
        ifs.close()

    else :
        print('file does not exist')
    return aData_out
