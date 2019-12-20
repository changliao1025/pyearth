import os
def read_configuration_file(sFilename_configuration_in):
    if os.path.isfile(sFilename_configuration_in):
        aConfig_out = {}
        ifs = open(sFilename_configuration_in, 'r')
        for sLine in ifs:
            sDummy = sLine.split(',')
            if (len(sDummy) == 2):
                print(sDummy)
                sKey = (sDummy[0]).strip()
                sValue = (sDummy[1]).strip()
                aConfig_out[sKey] = sValue
            else:
                pass
        ifs.close()
        return aConfig_out
    else:
        print('File does not exist: ' + sFilename_configuration_in)
        return -1
    
    
    