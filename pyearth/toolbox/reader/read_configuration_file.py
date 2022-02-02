import os
def read_configuration_file(sFilename_configuration_in):
    """
    Read a text based configuration file

    Args:
        sFilename_configuration_in (string): The filename of the configuration file

    Returns:
        dict: The content in the configuration file
    """

    if os.path.exists(sFilename_configuration_in):
        pass
    else:
        print('The xml file does not exist!')
        return
 
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
 
    
    
    
    