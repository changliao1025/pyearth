import os
import numpy
import sys
from netCDF4 import Dataset
#import library
sSystem_paths = os.environ['PATH'].split(os.pathsep)
sys.path.extend(sSystem_paths)
#import global variable
from eslib.system import define_global_variables
from eslib.system.define_global_variables import *
def add_variable_to_netcdf(sFilename_old, sFilename_new, aData_in, sVariable_in, sUnit_in, iDimension_in):
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
         
        if not iValue.isunlimited():
            pDatasets_out.createDimension(sKey, len(iValue) )
            if len(iValue) == iDimension_in:
                sDimension = sKey
                print('found')
        else:
            pDatasets_out.createDimension(sKey, len(iValue) )
        
    # Copy variables
    for sKey, aValue in pDatasets_in.variables.items():        
        print(aValue.datatype)
        print( aValue.dimensions)
        # we need to take care of rec dimension
        dummy = aValue.dimensions
        outVar = pDatasets_out.createVariable(sKey, aValue.datatype,   dummy  )
        
        for sAttribute in aValue.ncattrs():
            
            outVar.setncatts( { sAttribute: aValue.getncattr(sAttribute) } )

          
        outVar[:] = aValue[:]
        # close the output file
    #add new variable 
    pVar3 = pDatasets_out.createVariable(sVariable_in, 'f4', sDimension) 
    pVar3[:] = aData_in
    pVar3.description = sVariable_in
    pVar3.unit = sUnit_in

    pDatasets_out.close()

