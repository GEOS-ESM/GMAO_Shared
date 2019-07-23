'''
Different utils used by plotters.
'''

import mpl_toolkits.basemap as bm
import matplotlib.pyplot as pl
from matplotlib import colors, cm
import scipy as sp

# Make spectral plot
def my_psd(x):
    P,f=pl.psd(x); pl.clf()
    T=2./f; ind= T>=12.
    pl.plot(T[ind]/12,P[ind]);
    ax=pl.gca(); 
    ax.set_xscale('log');
    #ax.set_xlim(1,20); ax.set_ylim(0,20); 
    #tick=[1,2,3,5,10,20]; tickl=[str(i) for i in tick]
    #ax.set_xticks(tick); ax.set_xticklabels(tickl);
    ax.set_xlabel('Period, years'); ax.set_ylabel('Power');
    ax.set_title('Power spectral density');
    
# Contour plot
def contour(x,y,z,func=pl.contourf,black=None,**opts):
    '''
    Adds a "black" functionality to default contour function
    '''
    if black!=None:
        clevs=opts.get('levels',None)
        if clevs != None:
            min=clevs[0]; max=clevs[-1]
        else:
            min=sp.ma.minimum(z); max=sp.ma.maximum(z)
        norm=opts.get('norm',colors.normalize(min,max));
        cmap=opts.get('cmap',MyCmap(cm.get_cmap(),black=norm(black)))
        opts['norm']=norm; opts['cmap']=cmap
    cs=func(x,y,z,**opts)
    
    return cs

