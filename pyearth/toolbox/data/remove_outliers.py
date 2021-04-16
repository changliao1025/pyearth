import numpy as np
def remove_outliers(aData_in, outlierConstant):
    """
    Remove the outliers within the datesets
    """
    #remove nan first
    b = np.array(aData_in)
    c = b.ravel() #flatted array
    good_index =np.where( np.isfinite(c) == True  )
    a = c[good_index]
    upper_quartile = np.percentile(a, 95)
    lower_quartile = np.percentile(a, 5)
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    quartileSet = (lower_quartile + IQR, upper_quartile - IQR)    

    dummy_index = np.where( (aData_in>=quartileSet[0]) & (aData_in<=quartileSet[1]) )

    return aData_in[dummy_index]