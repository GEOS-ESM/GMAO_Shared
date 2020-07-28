import importlib
import matplotlib.pyplot as pl
import cartopy.crs as ccrs
import cmocean
import geosdset

def plot_clim(exp, ds):
    var='TS'
    FREEZE=273.16
    sst=ds[var]; sst-=FREEZE
    sst_clim=sst.groupby('time.season').mean('time')

    plotopts={'cmap': 'jet', 
              'vmin': -2, 
              'vmax': 32,
              'transform': ccrs.PlateCarree()}
    
    subplot_kws=dict(projection=ccrs.PlateCarree(),
                 facecolor='grey')

    cbar_kwargs={'orientation': 'horizontal',
              'shrink': 0.8,
              'extend': 'both'}

    pl.figure(1)
    sst_clim.sel(season='DJF').plot(**plotopts,subplot_kws=subplot_kws,cbar_kwargs=cbar_kwargs)
    pl.figure(2)
    sst_clim.sel(season='JJA').plot(**plotopts,subplot_kws=subplot_kws,cbar_kwargs=cbar_kwargs)
    pl.figure(3)
    sst.mean('time').plot(**plotopts,subplot_kws=subplot_kws,cbar_kwargs=cbar_kwargs)
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
    
