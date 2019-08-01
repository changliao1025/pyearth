import numpy as np
def choose_n_color(nColor_in, pColorMap_in):
    kn = np.arange(nColor_in)
    aColor=[]
    for k in kn:
        h,binEdges=np.histogram(data1[k])
        color=cmap(float(k)/kn.max())
        aColor.append[color]
    return aColor
