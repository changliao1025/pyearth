import os, sys
from netCDF4 import Dataset
from pyearth.system.define_global_variables import *

def add_variable_to_netcdf(sFilename_old, sFilename_new, aData_in, sVariable_in, sUnit_in, aDimension_in):
    if os.path.exists(sFilename_old):
        print("Yep, I can read that file!")
    else:
        print("Nope, the path doesn't reach your file. Go research filepath in python")
        exit

    pDatasets_in = Dataset(sFilename_old)
    netcdf_format = pDatasets_in.file_format
    #output file
    pDatasets_out = Dataset(sFilename_new, "w", format=netcdf_format)
    aDimension_key=list()
    aDimension_value=list()
    for sKey, iValue in pDatasets_in.dimensions.items():
        dummy = len(iValue)
        if not iValue.isunlimited():            
            aDimension_key.append(sKey)
            aDimension_value.append(sKey)
            pDatasets_out.createDimension(sKey, dummy)            
        else:
            pDatasets_out.createDimension(sKey, dummy )

    #
    aDimension_list=list()
    for d in aDimension_in:
        index = aDimension_value.index(d)
        aDimension_list.append(aDimension_key[index])

    aDimension_tuple = tuple(aDimension_list)

    #Copy variables
    for sKey, aValue in pDatasets_in.variables.items():        
        
        # we need to take care of rec dimension
        dummy = aValue.dimensions
        outVar = pDatasets_out.createVariable(sKey, aValue.datatype, dummy )        
        for sAttribute in aValue.ncattrs():            
            outVar.setncatts( { sAttribute: aValue.getncattr(sAttribute) } )

        outVar[:] = aValue[:]
        # close the output file
    #add new variable 


    pVar3 = pDatasets_out.createVariable(sVariable_in, 'f4', aDimension_tuple  ) 
    pVar3[:] = aData_in
    pVar3.description = sVariable_in
    pVar3.unit = sUnit_in
    pDatasets_out.close()

    return


def add_multiple_variable_to_netcdf(sFilename_old, sFilename_new, aData_in, aVariable_in, aUnit_in, aDimension_in):
    if os.path.exists(sFilename_old):
        print("Yep, I can read that file!")
    else:
        print("Nope, the path doesn't reach your file. Go research filepath in python")
        exit

    pDatasets_in = Dataset(sFilename_old)
    netcdf_format = pDatasets_in.file_format
    #output file
    pDatasets_out = Dataset(sFilename_new, "w", format=netcdf_format)
    aDimension_key=list()
    aDimension_value=list()
    for sKey, iValue in pDatasets_in.dimensions.items():
        dummy = len(iValue)
        if not iValue.isunlimited():            
            aDimension_key.append(sKey)
            aDimension_value.append(dummy)
            pDatasets_out.createDimension(sKey, dummy)            
        else:
            pDatasets_out.createDimension(sKey, dummy )

    

    #Copy variables
    for sKey, aValue in pDatasets_in.variables.items():      
        # we need to take care of rec dimension
        dummy = aValue.dimensions
        outVar = pDatasets_out.createVariable(sKey, aValue.datatype, dummy )        
        for sAttribute in aValue.ncattrs():            
            outVar.setncatts( { sAttribute: aValue.getncattr(sAttribute) } )

        outVar[:] = aValue[:]
        # close the output file
    #add new variable 

    nVariable = len(aData_in)
    for iVariable in range(nVariable):
        aDimension_variable = aDimension_in[iVariable]
        sVariable = aVariable_in[iVariable]
        sUnit=aUnit_in[iVariable]
        aData = aData_in[iVariable]
        aDimension_list=list()
        for d in aDimension_variable:
            index = aDimension_value.index(d)
            aDimension_list.append(aDimension_key[index])

        aDimension_tuple = tuple(aDimension_list)

        pVar3 = pDatasets_out.createVariable(sVariable, 'f4', aDimension_tuple  ) 
        pVar3[:] = aData
        pVar3.description = sVariable
        pVar3.unit = sUnit
        
    pDatasets_out.close()

    return