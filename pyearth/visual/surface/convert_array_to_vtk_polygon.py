import os, sys
import numpy as np
from pyevtk.hl import  unstructuredGridToVTK
from pyevtk.vtk import  VtkQuad

#convert a point based data to unstructured polygon vtk 
def convert_array_to_vtk_polygon(pArray_in, pLabel_in, aX_in, aY_in, pZ_in, sFilename_vtk):
    

    dMin_x = np.min(aX_in)
    dMax_x = np.max(aX_in)

    dMin_y = np.min(aY_in)
    dMax_y = np.max(aY_in)
    pShape_x = aX_in.shape
    #get resolution in x
    dResolution_x = ( dMax_x - dMin_x) / (pShape_x[1] -1 )
    pShape_y = aY_in.shape
    #get resolution in y
    dResolution_y = ( dMax_y - dMin_y ) / (pShape_y[0] -1 )
    pShape_z = pZ_in.shape
    #center
    nrow_center = pShape_y[0]
    ncolumn_center = pShape_x[1]
    nlayer_center = pShape_z[0] - 1
    #vertex
    nrow_vertex = nrow_center + 1
    ncolumn_vertex = ncolumn_center + 1
    nlayer_vertex = pShape_z[0] #including the base/reference 
    x = np.zeros(( nlayer_vertex , nrow_vertex, ncolumn_vertex))
    y = np.zeros(( nlayer_vertex , nrow_vertex, ncolumn_vertex))
    z = np.full( ( nlayer_vertex , nrow_vertex, ncolumn_vertex), np.nan, dtype=float )


    X = np.arange( dMin_x - 0.5 * dResolution_x, dMax_x + dResolution_x, dResolution_x, dtype='float64')
    Y = np.arange( dMax_y + 0.5 * dResolution_y, dMin_y - dResolution_y, -dResolution_y, dtype='float64')
    aGrid_x, aGrid_y = np.meshgrid(X, Y)


    iFlag =1
    for iLayer in np.arange(1, nlayer_vertex+1, 1):
        #x and y are ready 
        x[iLayer-1, :,: ] = aGrid_x
        y[iLayer-1, :,: ] = aGrid_y

        aZ = np.full( ( nrow_vertex, ncolumn_vertex), np.nan, dtype=float )

        aZ_subset = pZ_in[iLayer-1]
        #now deal with edges, there are two possible arrangments
        if iFlag ==1 : #this is a easy approach
            aZ[0:pShape_y[0],0:pShape_x[1] ] = aZ_subset
            #right edge
            for iRow in np.arange(0, nrow_center, 1 ):
                aZ[iRow, ncolumn_vertex-1] = aZ_subset[iRow, ncolumn_vertex-2]
                pass
            #bottm edge
            for iColumn in np.arange(0, ncolumn_center, 1 ):
                aZ[nrow_vertex-1, iColumn] = aZ_subset[nrow_vertex-2, iColumn]
            #lower right 
            aZ[nrow_center ,ncolumn_center  ] = aZ_subset[nrow_vertex-2, ncolumn_vertex-2]           
            pass
        else:
            #in this case, we set the origin at lower left
            pass

        #set z
        z[iLayer-1,:,:] = aZ     

        pass
       
    # Define connectivity or vertices that belongs to each element
    conn = np.full(( nrow_center * ncolumn_center * nlayer_center * 4), 0, dtype=int)

    #without bottom layer
    for iLayer in np.arange(1, nlayer_center + 1, 1):
        start_id = (iLayer-1) * nrow_vertex * ncolumn_vertex
        start_index = (iLayer-1) * nrow_center * ncolumn_center * 4

        for iRow in np.arange(1, nrow_center+1, 1):
            for iColumn in np.arange(1, ncolumn_center+1, 1 ):
                
                new_index = (iRow-1) * ncolumn_center * 4 + (iColumn-1) * 4

                conn[start_index + new_index] = (iRow-1) * ncolumn_vertex + iColumn-1  + start_id
                conn[start_index + new_index+1] = (iRow-1) * ncolumn_vertex + iColumn + start_id
                conn[start_index + new_index+2] = iRow * ncolumn_vertex + iColumn + start_id
                conn[start_index + new_index+3] = iRow * ncolumn_vertex + iColumn-1 + start_id

                pass  
        pass

    #we don not consider bottom base layer connectivity

    # Define offset of last vertex of each element
    offset = np.arange( 4, (4*ncolumn_center*nrow_center*nlayer_center+1), 4 )    
    # Define cell types
    ctype = np.full( (nrow_center, ncolumn_center, nlayer_center), VtkQuad.tid )


    x.shape =  (nrow_vertex) * (ncolumn_vertex) * (nlayer_vertex)
    y.shape =  (nrow_vertex) * (ncolumn_vertex) * (nlayer_vertex)
    z.shape =  (nrow_vertex) * (ncolumn_vertex) * (nlayer_vertex)

    celldata = {}
    nCelldata = len(pLabel_in)
    ncell = nrow_center * ncolumn_center * nlayer_center

    for i in np.arange(nCelldata):
        sLabel = pLabel_in[i]
        #reshape 
        data=pArray_in[i]
        data.shape = ncell
        celldata[sLabel] = data
        pass
    

    conn= np.reshape(conn, (4*nrow_center*ncolumn_center*nlayer_center))
    #save as vtk file


    unstructuredGridToVTK(sFilename_vtk, x, y, z, connectivity = conn, \
                              offsets = offset, cell_types = ctype, cellData = celldata, pointData = None)



    return