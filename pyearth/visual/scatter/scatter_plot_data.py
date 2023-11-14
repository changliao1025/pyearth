
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyearth.visual.scatter.scatter_lowess import scatter_lowess

def scatter_plot_data(aData_x,  
                      aData_y, 
                      sFilename_out,  
                      iFlag_scientific_notation_x_in=None, 
                      iFlag_scientific_notation_y_in=None, 
                      iSize_x_in = None,  
                      iSize_y_in = None,   
                      iDPI_in = None , 
                      iFlag_log_x_in = None, 
                      iFlag_log_y_in = None, 
                      iFlag_lowess_in = None, 
                      dMin_x_in = None,  
                      dMax_x_in = None,  
                      dMin_y_in = None,  
                      dMax_y_in = None,  
                      dSpace_x_in = None,  
                      dSpace_y_in = None,  
                      sFormat_x_in =None, 
                      sFormat_y_in =None, 
                      sLabel_x_in =None, 
                      sLabel_y_in = None ,  
                      sLabel_legend_in = None, 
                      sTitle_in = None):
    try:
        import scipy
    except ImportError as e:
        raise ImportError("The package 'scipy' is required for this function to run.") from e

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 12

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 12

    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_scientific_notation_x_in is not None:
        iFlag_scientific_notation_x = iFlag_scientific_notation_x_in
    else:
        iFlag_scientific_notation_x = 0

    if iFlag_scientific_notation_y_in is not None:
        iFlag_scientific_notation_y = iFlag_scientific_notation_y_in
    else:
        iFlag_scientific_notation_y = 0

    if iFlag_log_x_in is not None:
        iFlag_log_x = iFlag_log_x_in
    else:
        iFlag_log_x = 0

    if iFlag_log_y_in is not None:
        iFlag_log_y = iFlag_log_y_in
    else:
        iFlag_log_y = 0

    
    if iFlag_log_y_in is not None:
        iFlag_lowess = iFlag_lowess_in
    else:
        iFlag_lowess = 0 

    if sLabel_x_in is not None:
        sLabel_X = sLabel_x_in
    else:
        sLabel_X = ''

    if sLabel_y_in is not None:
        sLabel_Y = sLabel_y_in
    else:
        sLabel_Y = ''

    if sLabel_legend_in is not None:
        sLabel_legend = sLabel_legend_in
    else:
        sLabel_legend = ''

    if sTitle_in is not None:
        sTitle = sTitle_in
    else:
        sTitle = ''

    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )

    left, width = 0.1, 0.75
    bottom, height = 0.1, 0.75
    spacing = 0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.15]
    rect_histy = [left + width + spacing, bottom, 0.15, height]

    
    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    iFlag_histgram=1
    if iFlag_histgram ==1:
        ax_histx = plt.axes(rect_histx)
        ax_histx.tick_params(direction='in', labelbottom=False)
        ax_histy = plt.axes(rect_histy)
        ax_histy.tick_params(direction='in', labelleft=False)


    nPoint = len(aData_x)
    x_min = np.nanmin(aData_x)
    x_max = np.nanmax(aData_x)
    y_min = np.nanmin(aData_y)
    y_max =np.ceil(np.nanmax(aData_y))

    x = aData_x
    y = aData_y

    aLegend_artist = []
    aLegend_label=[]

    
    cmap = plt.get_cmap('BuPu')

    sc= ax_scatter.scatter(x, y,  alpha=0.5,cmap=cmap)
    #ax_scatter.set_facecolor('silver')
    aLegend_artist.append(sc)
    sLabel = sLabel_legend
    aLegend_label.append(sLabel)

    ax_scatter.axis('on')
    ax_scatter.grid(which='major', color='grey', linestyle='--', axis='y')

    ax_scatter.tick_params(axis="x", labelsize=13)
    ax_scatter.tick_params(axis="y", labelsize=13)

    ax_scatter.set_xmargin(0.05)
    ax_scatter.set_ymargin(0.05)

    ax_scatter.set_xlabel(sLabel_X,fontsize=12)
    ax_scatter.set_ylabel(sLabel_Y,fontsize=12)
    
    # round to nearest years...

    if sFormat_x_in is not None:
        sFormat_x=sFormat_x_in
        ax_scatter.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(sFormat_x))

        #ax_scatter.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1e'))'%.1f'

    if sFormat_y_in is not None:
        #sFormat_y = sFormat_y_in
        #ax_scatter.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(sFormat_y))
        sFormat_y_dummy =  sFormat_y_in.replace("{", "{x")
        ax_scatter.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter(  sFormat_y_dummy ) ) 

        #ax_scatter.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))'%.1f'

    ax_scatter.tick_params(axis='y', pad=8)
    if dMin_x_in is not None:
        dMin_x = dMin_x_in
    else:
        dMin_x = x_min

    if dMax_x_in is not None:
        dMax_x = dMax_x_in
    else:
        dMax_x = x_max

    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        dMin_y = y_min

    if dMax_y_in is not None:
        dMax_y = dMax_y_in
    else:
        dMax_y = y_max

    if dSpace_x_in is not None:
        dSpace_x = dSpace_x_in
    else:
        dSpace_x = (dMax_x - dMin_x) / 4.0

    if dSpace_y_in is not None:
        dSpace_y = dSpace_y_in
    else:
        dSpace_y = (dMax_y - dMin_y) / 4.0

    ax_scatter.set_xlim( dMin_x, dMax_x )
    ax_scatter.set_ylim( dMin_y, dMax_y)


    #ax_scatter.xaxis.set_major_locator(mpl.ticker.MaxNLocator(prune='upper', nbins=4))
    ax_scatter.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base = dSpace_x))


    if iFlag_log_x ==1:
        aLabel_x = []
        for i in np.arange( dMin_x, dMax_x +1, 1 ):
            sTicklabel = r'$10^{{{}}}$'.format(int(i))
            aLabel_x.append(sTicklabel)
            pass

        ax_scatter.set_xticks(np.arange( dMin_x, dMax_x +1, 1 ))
        ax_scatter.set_xticklabels(aLabel_x)
    else:
        if iFlag_scientific_notation_x ==1:
            formatter = mpl.ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            #formatter.set_powerlimits((-1,6)) # you might need to change here
            ax_scatter.xaxis.set_major_formatter(formatter)
        else:
            pass


    ax_scatter.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base = dSpace_y))


    if iFlag_log_y ==1:
        aLabel_y = []
        for i in np.arange( dMin_y, dMax_y +1, 1 ):
            sTicklabel = r'$10^{{{}}}$'.format(int(i))
            aLabel_y.append(sTicklabel)
            pass

        ax_scatter.set_yticks(np.arange( dMin_y, dMax_y +1, 1 ))
        ax_scatter.set_yticklabels(aLabel_y)
        pass
    else:
        if iFlag_scientific_notation_y ==1:
            formatter = mpl.ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            #formatter.set_powerlimits((-1,6)) # you might need to change here
            ax_scatter.yaxis.set_major_formatter(formatter)

        pass

    dRatio = (float(iSize_y)/iSize_x) / ( (dMax_y-dMin_y )/ ( dMax_x-dMin_x ) )
    
    
    if(iFlag_lowess==1):      

        y_sm, y_std, order = scatter_lowess(aData_x, aData_y, f=1./3.)
        ax_scatter.plot(x[order], y_sm[order], color='tomato')
        sc_lowess =ax_scatter.fill_between(x[order],  
                                y_sm[order] - 1.96*y_std[order],  
                                y_sm[order] + 1.96*y_std[order],  
                                alpha=0.3)

        aLegend_artist.append(sc_lowess)

        sLabel_legend_lowess = 'LOWESS uncertainty'
      
        aLegend_label.append(sLabel_legend_lowess)

    ax_scatter.set_aspect(dRatio, 'box')  #this one set the y / x ratio
    ax_scatter.legend(aLegend_artist, aLegend_label, 
                      loc="upper right", fontsize=12)

    
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
    
    
    print(slope, intercept, r_value, p_value, std_err)
    sText = r'R: ' + "{:.2f}".format( r_value )
    ax_scatter.text(0.05, 0.95, sText,  
    verticalalignment='top', horizontalalignment='left', 
            transform=ax_scatter.transAxes,  
            color='black', fontsize=12)
    
    sText = 'P-value: ' + "{:.2E}".format( p_value )
    ax_scatter.text(0.05, 0.9, sText,  
    verticalalignment='top', horizontalalignment='left', 
            transform=ax_scatter.transAxes,  
            color='black', fontsize=12)

    ax_scatter.tick_params(which='both', # Options for both major and minor ticks
                           top='off', # turn off top ticks
                           left='off', # turn off left ticks
                           right='off',  # turn off right ticks
                           bottom='off') # turn off bottom ticks
    
    if iFlag_histgram ==1:
        density = scipy.stats.gaussian_kde(x)
        xx = np.linspace(dMin_x, dMax_x,1000)
        yy = density(xx)
        ax_histx.plot(xx,yy, color='navy')
        ax_histx.fill_between(xx, yy, 0, linewidth=3,  color = 'lightblue')

        ax_histx.set_xlim( dMin_x, dMax_x )
        ax_histx.set_ylim( 0, auto=None )

        ax_histx.axis('on')
        ax_histx.grid(which='major', color='white', linestyle='-', axis='x')
        ax_histx.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base = dSpace_x/2.0))
        ax_histx.spines['right'].set_visible(False)
        ax_histx.spines['top'].set_visible(False)
        ax_histx.spines['bottom'].set_visible(False)
        ax_histx.spines['left'].set_visible(False)
        ax_histx.tick_params(which='both', # Options for both major and minor ticks
                             top='off', # turn off top ticks
                             left='off', # turn off left ticks
                             right='off',  # turn off right ticks
                             bottom='off') # turn off bottom ticks

        ax_histx.axes.get_yaxis().set_visible(False)
        ax_histx.tick_params(axis='x', colors='white')

        #y margin
    
        density = scipy.stats.gaussian_kde(y)
        xx = np.linspace(dMin_y, dMax_y,1000)
        yy = density(xx)
        xx, yy = yy, xx
        ax_histy.plot(xx,yy, color='navy')
        ax_histy.fill_betweenx(yy, 0, xx, linewidth=3,  color = 'lightblue')

        ax_histy.set_xlim(0, auto=None)
        ax_histy.set_ylim(dMin_y, dMax_y)

        ax_histy.axis('on')
        ax_histy.grid(which='major', color='white', linestyle='-', axis='y')
        ax_histy.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base = dSpace_y/2.0))
        ax_histy.spines['right'].set_visible(False)
        ax_histy.spines['top'].set_visible(False)
        ax_histy.spines['bottom'].set_visible(False)
        ax_histy.spines['left'].set_visible(False)
        ax_histy.axes.get_xaxis().set_visible(False)


        ax_histy.tick_params(axis='y', colors='white')

        ax_histy.tick_params(which='both', # Options for both major and minor ticks
                             top='off', # turn off top ticks
                             left='off', # turn off left ticks
                             right='off',  # turn off right ticks
                             bottom='off') # turn off bottom ticks

    ax_scatter.set_title( sTitle, loc='center', fontsize=15)
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
   
