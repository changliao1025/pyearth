import numpy as np
def remove_outliers(x, outlierConstant):
    #remove nan first
    b = np.array(x)
    c = b.ravel() 
    good_index =np.where( np.isfinite(c) == True  )
    a = c[good_index]
    upper_quartile = np.percentile(a, 75)
    lower_quartile = np.percentile(a, 25)
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    quartileSet = (lower_quartile + IQR, upper_quartile - IQR)    

    dummy_index = np.where( (x>=quartileSet[0]) & (x<=quartileSet[1]) )

    return x[dummy_index]