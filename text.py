import os


def text_reader_string( sFilename_in, ncolumn_in = None,  nrow_in = None, delimiter_in = None, skipline_in =  None ):
    print(sFilename_in)
    print(os.path.isfile(sFilename_in))
    if os.path.isfile(sFilename_in):
        print('file exist')


        #check the parameter
        if ncolumn_in is not None:
            iFlag_column = 1
        else:
            iFlag_column = 0
        if nrow_in is not None :
            #there is nothing
            print(' ')
        else :
            #get the total line number
            ifs = open(sFilename_in, "r")
            nrow_in = len(ifs.readlines())
            ifs.close()
       
        line=' '
        ifs = open(sFilename_in, "r")

        if skipline_in is not None:
            nrow_in =  nrow_in - skipline_in
            for i in range(skipline_in):
                ifs.readline()
            
        else:
            print('stupid')
        #get delimiter
        if delimiter_in is not None:
            iFlag_delimiter = 1
            
        else:
            iFlag_delimiter = 0
            delimiter_in = ' '
        line = (ifs.readline()).rstrip()

        if iFlag_delimiter == 0:
            if iFlag_column == 1:
                print('stupid')
            else :
                temp = line.split(delimiter_in)
                
                ncolumn_in = len(temp)
      
            #check ncolumn_in count
            if ncolumn_in < 1 :
                print( 'The file has no data!' )
                        
            data_out = [[0 for x in range(ncolumn_in)] for y in range(nrow_in)] 

            print(type(data_out))
            print(data_out)
            temp = line.split(delimiter_in)
            print(temp)
            print(type(temp))
            data_out[0] =  temp

            for row in range(1, nrow_in):
                line=(ifs.readline()).rstrip()
                data_out[row] = line.split(delimiter_in)
            
      

        else :
            if iFlag_column == 1:
                print('stupid')
            else :
                temp = line.split(delimiter_in)
                
                ncolumn_in = len(temp)
      
            #check ncolumn_in count
            if ncolumn_in < 1 :
                print( 'The file has no data!' )                        
            data_out = [[0 for x in range(nrow_in)] for y in range(ncolumn_in)] 

            temp = line.split(delimiter_in)
            data_out[0, :] =  temp

            for row in range(1,nrow_in):
                line = (ifs.readline()).rstrip()
                data_out[row, :] = line.split(delimiter_in)



        print(data_out)
        ifs.close()

    else :
        print('file does not exist')

    return data_out