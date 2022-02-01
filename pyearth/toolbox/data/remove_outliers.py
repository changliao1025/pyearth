import copy
import numpy as np
def remove_outliers(aData_in, dOutlier_percentage):
    """
    Remove the outliers within the datesets
    Args:
        aData_in (numpy.array): The array data
        dOutlier_percentage (float): The percentage to be marked as outlier

    Returns:
        numpy.array: The data without outlier
    """
    
    #remove nan first    
    a = copy.deepcopy(aData_in) 
    b = np.array(a)
    c = b.ravel() #flatted array
    good_index =np.where( np.isfinite(c) == True  )
    d = c[good_index]
    upper_quartile = np.percentile(d, 95)
    lower_quartile = np.percentile(d, 5)
    dQR = (upper_quartile - lower_quartile) * dOutlier_percentage
    quartileSet = (lower_quartile + dQR, upper_quartile - dQR)    

    dummy_index = np.where( (b>=quartileSet[0]) & (b<=quartileSet[1]) )

    return aData_in[dummy_index]