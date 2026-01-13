import os


def replace_variable_in_netcdf(sFilename_old, sFilename_new, aData_in, aVariable_in):
    try:
        import netCDF4 as nc
    except ImportError as e:
        raise ImportError(
            "The package 'netCDF4' is required for this function to run."
        ) from e

    if os.path.exists(sFilename_old):
        print("Yep, I can read that file: " + sFilename_old)
    else:
        print("Nope, the path doesn't reach your file. Go research filepath in python")
        exit
    pDatasets_in = nc.Dataset(sFilename_old)
    netcdf_format = pDatasets_in.file_format
    # output file
    pDatasets_out = nc.Dataset(sFilename_new, "w", format=netcdf_format)
    for sKey, iValue in pDatasets_in.dimensions.items():
        dummy = len(iValue)
        if not iValue.isunlimited():
            pDatasets_out.createDimension(sKey, dummy)
        else:
            pDatasets_out.createDimension(sKey, dummy)

    # copy global attribute

    for attr in pDatasets_in.ncattrs():
        # print(':::GlobalAtt:', attr,' Val:', getattr(pDatasets_in,attr))
        pDatasets_out.setncattr(attr, pDatasets_in.getncattr(attr))

    # Copy variables
    aAttribute = list()
    aDimension = list()
    aDataType = list()
    for sKey, aValue in pDatasets_in.variables.items():

        # we need to take care of rec dimension

        if sKey in aVariable_in:  # only copy other dataset first
            aDataType.append(aValue.datatype)
            aDimension.append(aValue.dimensions)
            # pUnit = aValue.units
            for sAttribute in aValue.ncattrs():
                if sAttribute == "_FillValue" or sAttribute == "missing_value":
                    pass
                else:
                    aAttribute.append(sAttribute)
            pass
        else:

            outVar = pDatasets_out.createVariable(
                sKey, aValue.datatype, aValue.dimensions
            )
            for sAttribute in aValue.ncattrs():
                outVar.setncatts({sAttribute: aValue.getncattr(sAttribute)})

            outVar[:] = aValue[:]

    # close the output file
    # replace variable
    for i in range(len(aVariable_in)):
        sVariable_in = aVariable_in[i]
        aData_in0 = aData_in[i]
        pVar3 = pDatasets_out.createVariable(
            sVariable_in, aDataType[i], aDimension[i], fill_value=-9999
        )
        pVar3[:] = aData_in0
        pValue = pDatasets_in.variables[sVariable_in]
        for sAttribute in pValue.ncattrs():
            if sAttribute == "_FillValue" or sAttribute == "missing_value":
                pass
            else:
                pVar3.setncatts({sAttribute: pValue.getncattr(sAttribute)})
        pVar3.missing_value = -9999

    # pVar3.description = sVariable_in
    # pVar3.units = pUnit
    pDatasets_out.close()
    return
