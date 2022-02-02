import os
import numpy as np
import xml.etree.ElementTree as ET

def parse_xml_file(sFilename_xml_in):
    """
    Parse an XML file

    Args:
        sFilename_xml_in (string): The XML filename

    Returns:
        dict: The content in the XML file
    """

    if os.path.exists(sFilename_xml_in):
        pass
    else:
        print('The xml file does not exist!')
        return
    tree = ET.parse(sFilename_xml_in)
    root = tree.getroot()
    namelist = {}
    keys_str = []

    for entry in root.findall('./group/entry'):
        key = entry.get('id')
        dtype = entry.find('type').text
        if dtype=='char':
            keys_str.append(key)
        if 'value' in entry.keys():
            # single value setting
            if dtype=='integer':
                value = int(entry.get('value'))
            elif dtype=='real':
                value = float(entry.get('value'))
            elif dtype=='logical':
                if entry.get('value')=='TRUE':
                    value = True
                else:
                    value = False
            else:
                value = entry.get('value')
            namelist[key] = value
        else:
            # value is an array
            values = []
            for value in entry.findall('./values'):
                key = entry.get('id')
                dtype = entry.find('type').text
                if dtype=='integer':
                    values.append(int(value.text))
                elif dtype=='real':
                    values.append(float(value.text))
                elif dtype=='logical':
                    if value.text=='TRUE':
                        values.append(True)
                    else:
                        values.append(False)
                else:
                    values.append(value.text)
            namelist[key] = np.array(values, dtype=np.float64, order='F')

    # fill environment variable values
    for key in keys_str:
        keys_str_iter = list(set(keys_str)-set([key]))
        for okey in keys_str_iter:
            string = namelist[okey]
            if string.find('$'+key)>=0:
                new_string = string.replace('$'+key, namelist[key])
                namelist[okey] = new_string
                
    return namelist