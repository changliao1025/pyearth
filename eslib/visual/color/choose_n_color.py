import numpy as np
def choose_n_color(nColor_in, pColorMap_in):
    kn = np.arange(nColor_in) 

    aColor = [pColorMap_in(float(k)/kn.max()) for k in kn]
    return aColor