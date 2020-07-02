import numpy as np
def remove_outliers(x, outlierConstant):
    a = np.array(x)
    upper_quartile = np.percentile(a, 75)
    lower_quartile = np.percentile(a, 25)
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    quartileSet = (lower_quartile + IQR, upper_quartile - IQR)
    #resultList = []
    #for y in a.tolist():
    #    if y >= quartileSet[0] and y <= quartileSet[1]:
    #        resultList.append(y)

    dummy_index = np.where( (x>=quartileSet[0]) & (x<=quartileSet[1]) )

    return x[dummy_index]