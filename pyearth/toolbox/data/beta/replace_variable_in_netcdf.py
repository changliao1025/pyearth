import os, sys
from netCDF4 import Dataset
def replace_variable_in_netcdf(sFilename_old, sFilename_new, aData_in, sVariable_in):

    if os.path.exists(sFilename_old):
        print("Yep, I can read that file!")
    else:
        print("Nope, the path doesn't reach your file. Go research filepath in python")
        exit
    pDatasets_in = Dataset(sFilename_old)
    netcdf_format = pDatasets_in.file_format
    #output file
    pDatasets_out = Dataset(sFilename_new, "w", format=netcdf_format)
    for sKey, iValue in pDatasets_in.dimensions.items():
        dummy = len(iValue)
        if not iValue.isunlimited():            
            pDatasets_out.createDimension(sKey, dummy)            
        else:
            pDatasets_out.createDimension(sKey, dummy )

    # Copy variables
    for sKey, aValue in pDatasets_in.variables.items():        
        print(aValue.datatype)
        print( aValue.dimensions)
        # we need to take care of rec dimension
        dummy = aValue.dimensions
        if(sKey != sVariable_in): #only copy other dataset first
            outVar = pDatasets_out.createVariable(sKey, aValue.datatype, dummy )        
            for sAttribute in aValue.ncattrs():            
                outVar.setncatts( { sAttribute: aValue.getncattr(sAttribute) } )

            outVar[:] = aValue[:]
        else:
            pDataType = aValue.datatype
            pDimension = aValue.dimensions
            pUnit = aValue.units

            pass
        # close the output file
    #replace variable 
    pVar3 = pDatasets_out.createVariable(sVariable_in, pDataType, pDimension  ) 
    pVar3[:] = aData_in
    pVar3.description = sVariable_in
    pVar3.units = pUnit
    pDatasets_out.close()
    return