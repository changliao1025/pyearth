import os
from pyearth.system.define_global_variables import *

def add_variable_to_netcdf(sFilename_old, sFilename_new, aData_in, aVariable_in, aUnit_in, aaDimension_in):

    try:
        import netCDF4 as nc
    except ImportError as e:
        raise ImportError("The package 'netCDF4' is required for this function to run.") from e

    if os.path.exists(sFilename_old):
        print("Yep, I can read that file!")
    else:
        print("Nope, the path doesn't reach your file. Go research filepath in python")
        exit

    if os.path.exists(sFilename_new):
        os.remove(sFilename_new)

    pDatasets_in = nc.Dataset(sFilename_old)
    netcdf_format = pDatasets_in.file_format
    #output file
    pDatasets_out = nc.Dataset(sFilename_new, "w", format=netcdf_format)
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
    aaDimension_list=list()
    for aDimension_in in aaDimension_in:
        aDimension_list=list()
        for d in aDimension_in:
            index = aDimension_value.index(d)
            aDimension_list.append(aDimension_key[index])
        aaDimension_list.append(aDimension_list)

    aaDimension_tuple = tuple(aaDimension_list)

    #Copy variables
    for sKey, aValue in pDatasets_in.variables.items():
        # we need to take care of rec dimension
        dummy = aValue.dimensions
        outVar = pDatasets_out.createVariable(sKey, aValue.datatype, dummy, fill_value=-9999  )
        for sAttribute in aValue.ncattrs():
            if sAttribute != '_FillValue' and sAttribute != 'missing_value':
                outVar.setncatts( { sAttribute: aValue.getncattr(sAttribute) } )


        outVar[:] = aValue[:]
        # close the output file
    #add new variable

    nData = len(aData_in)
    for iData in range(nData):
        sVariable_in = aVariable_in[iData]
        sUnit_in = aUnit_in[iData]
        aData_in = aData_in[iData]
        aDimension_tuple = aaDimension_tuple[iData]

        pVar3 = pDatasets_out.createVariable(sVariable_in, 'f4', aDimension_tuple  )
        pVar3[:] = aData_in
        pVar3.description = sVariable_in
        pVar3.unit = sUnit_in
    pDatasets_out.close()

    return


def add_multiple_variable_to_netcdf(sFilename_old, sFilename_new, aData_in, aVariable_in, aUnit_in, aDimension_in):
    try:
        import netCDF4 as nc
    except ImportError as e:
        raise ImportError("The package 'netCDF4' is required for this function to run.") from e

    if os.path.exists(sFilename_old):
        print("Yep, I can read that file!")
    else:
        print("Nope, the path doesn't reach your file. Go research filepath in python")
        exit
    missing_value=-9999
    pDatasets_in = nc.Dataset(sFilename_old)
    netcdf_format = pDatasets_in.file_format
    #output file
    pDatasets_out = nc.Dataset(sFilename_new, "w", format=netcdf_format)
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

        pVar3 = pDatasets_out.createVariable(sVariable, 'f4', aDimension_tuple,fill_value=missing_value  )
        pVar3[:] = aData
        pVar3.description = sVariable
        pVar3.unit = sUnit

    pDatasets_out.close()

    return