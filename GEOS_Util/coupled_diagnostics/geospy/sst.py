import importlib
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean
import geosdset

def plot_clim(exp, ds):
    var='TS'
    FREEZE=273.16
    sst=ds[var]; sst-=FREEZE
    sst_clim=sst.groupby('time.season').mean('time')

    subplot_kws={'projection': ccrs.PlateCarree(),
                 'facecolor':'grey'}

    cbar_kwargs={'orientation': 'horizontal',
              'shrink': 0.8,
              'extend': 'both'}

    plotopts={'cmap': cmocean.cm.thermal, 
              'levels': np.arange(0.0,31.0,2.0),
              'transform': ccrs.PlateCarree(),
              'cbar_kwargs': cbar_kwargs}    

    ff=pl.figure(1) 
    ax=pl.axes(**subplot_kws)
    pp=sst_clim.sel(season='DJF').plot.contourf(ax=ax,**plotopts)
    ax.coastlines()
    ax.stock_img()
    gl=ax.gridlines()
    gl.ylabbels_left=True
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
#    pl.figure(2)
#    sst_clim.sel(season='JJA').plot(**plotopts)
#    pl.figure(3)
#    sst.mean('time').plot(**plotopts)
    pl.show()

def plot_diff(exp, ds1, ds2):
    pass

def plots(exps, dsets):
    obs=[]
    plot_clim(exps[0] ,dsets[0])

    for ds in dsets:
        plot_diff(exps[0], dsets[0], ds)
        
    for ds in obs:
        plot_diff(exps[0], dsets[0], ds)

if __name__=='__main__':
    import sys
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    plots(exps,dsets)
    
