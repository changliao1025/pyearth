import numpy as np

from scipy.stats import gaussian_kde

def kde_support(aData_in, bw_in, iGridsize_in, cut_in, aClip_in):
    """
    Establish support for a kernel density estimate.

    Args:
        aData_in (numpy.array): Input data
        bw_in (string): bandwidth
        iGridsize_in (int): Gridsize
        cut_in (float): Factor
        aClip_in (numpy.array): Bound

    Returns:
        numpy.array: kde
    """
    support_min = max(aData_in.min() - bw_in * cut_in, aClip_in[0])
    support_max = min(aData_in.max() + bw_in * cut_in, aClip_in[1])
    return np.linspace(support_min, support_max, iGridsize_in)

def scipy_bivariate_kde(aX_in, aY_in, bw_in, iGridsize_in, cut_in, aClip_in):
    """
    Compute a bivariate kde using scipy.

    Args:
        aX_in (numpy.array): X array    
        aY_in (numpy.array): Y array
        bw_in (string): bandwidth
        iGridsize_in (int): Gridsize
        cut_in (float): Factor
        aClip_in (numpy.array): Bound

    Raises:
        

    Returns:
        numpy.array: kde
    """    
    

    data = np.c_[aX_in, aY_in]
    kde = gaussian_kde(data.T, bw_method=bw_in)
    data_std = data.std(axis=0, ddof=1)
    if isinstance(bw_in, str):
        bw_in = "scotts" if bw_in == "scott" else bw_in
        bw_x = getattr(kde, "%s_factor" % bw_in)() * data_std[0]
        bw_y = getattr(kde, "%s_factor" % bw_in)() * data_std[1]
    elif np.isscalar(bw_in):
        bw_x, bw_y = bw_in, bw_in
    else:
        msg = ("Cannot specify a different bandwidth for each dimension "
               "with the scipy backend. You should install statsmodels.")
        raise ValueError(msg)
        
    x_support = kde_support(data[:, 0], bw_x, iGridsize_in, cut_in, aClip_in[0])
    y_support = kde_support(data[:, 1], bw_y, iGridsize_in, cut_in, aClip_in[1])
    xx, yy = np.meshgrid(x_support, y_support)
    z = kde([xx.ravel(), yy.ravel()]).reshape(xx.shape)
    return xx, yy, z