import os
def text_reader_string( sFilename_in,\
     ncolumn_in = None, \
     nrow_in = None, \
     delimiter_in = None, \
     skipline_in =  None ):
    #print(sFilename_in)
    #print(os.path.isfile(sFilename_in))
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
        #get delimiter
        if delimiter_in is not None:
            iFlag_delimiter = 1
        else:
            iFlag_delimiter = 0
            delimiter_in = ' '
        sLine = (ifs.readline()).rstrip()

        if iFlag_delimiter == 0:
            if iFlag_column == 1:
                pass
            else :
                dummy = sLine.split(delimiter_in)
                ncolumn_out = len(dummy)
            #check ncolumn_in count
            if ncolumn_out < 1 :
                print( 'The file has no data!' )
                        
            aData_out = [[0 for x in range(ncolumn_out)] for y in range(nrow_out)]             
            dummy = sLine.split(delimiter_in)
           
            aData_out[0] =  dummy

            for iRow in range(1, nrow_out):
                sLine=(ifs.readline()).rstrip()
                aData_out[iRow] = sLine.split(delimiter_in)
        else :
            if iFlag_column == 1:
                pass
            else :
                dummy = sLine.split(delimiter_in)
                ncolumn_out = len(dummy)
            #check ncolumn_in count
            if ncolumn_out < 1 :
                
                pass                      
            aData_out = [[0 for x in range(ncolumn_out)] for y in range(nrow_out)] 

            dummy = sLine.split(delimiter_in)
           
            aData_out[0] =  dummy

            for iRow in range(1,nrow_out):
                dummy1 = ifs.readline()
                #print('New line is: ' + dummy1)
                dummy2 = dummy1.rstrip()
                dummy3 = dummy2.split(delimiter_in)
                #print(dummy3)
                aData_out[iRow] = dummy3
        #print(aData_out)
        ifs.close()

    else :
        print('file does not exist')
    return aData_out
