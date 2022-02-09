#!/usr/bin/env python3

'''
Calculates AMOC stream function.
'''
import sys
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import geosdset, plots

plotname='AMOC'
defaults={'name': 'ty_trans', 
          'colname': 'ocean_month', 
          'coltype': 'MOM'}

def mkamoc(exp,dset):
    vardata=exp['plots'].get(plotname,defaults)
    varname=vardata['name']

    amoc=dset[varname]*dset['atl_mask'].values
    amoc=amoc.sum('xt_ocean',skipna=True)
    # We need to integrate transport from the bottom using reverse cumsum
    amoc=-amoc[:,::-1,:].cumsum('st_ocean',skipna=True)[:,::-1,:]

    # Add GM transport
    try:
        gm=dset[varname]*dset['atl_mask'].values
        gm=gm.sum('xt_ocean',skipna=True)
    except KeyError:
        gm=0.0

    return amoc+gm

def plot_amoc(exp,da):

    cbar_kwargs={'label': 'SV',
                 'orientation': 'horizontal',
                 'shrink': 0.8}
    
    fill_opts={'yincrease': False,
               'levels': np.arange(-20.0,20.1,2.0),
               'cmap': cmocean.cm.balance,
               'cbar_kwargs': cbar_kwargs}

    contour_opts={'levels': np.arange(-20.0,20.1,2.0),
                  'colors': 'black',
                  'yincrease': False}

    clab_opts={}
    plotter=plots.Plot2d(fill_opts=fill_opts, contour_opts=contour_opts, clab_opts=clab_opts)

    pl.figure(1); pl.clf()
    ax=plotter.contour(da.sel(time=slice(*exp['dates'])).mean('time'))
    ax.set_title(f'{exp["expid"]} AMOC stream function')
    ax.set_xlim(-30.0,65.0)
    ax.set_ylabel('depth')
    ax.set_xlabel('')
    ax.xaxis.set_major_formatter(plots.LATITUDE_FORMATTER)
    ax.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/atl_moc.png')
    pl.show()

def plot_amoc_ind(exp,da):
    ind40=da.sel(yu_ocean=slice(40.,50.), st_ocean=slice(1000.,1500)).mean(('yu_ocean','st_ocean'))
    ind26=da.sel(yu_ocean=26., st_ocean=1000., method='nearest')
    plotter=plots.Plot1d()
    pl.figure(2); pl.clf();
    plotter.line_opts['label']='40N-50N'
    ax=plotter.line(ind40)
    plotter.line_opts['label']='26N'
    plotter.line(ind26,ax=ax)
    ax.legend()
    ax.set_title(f'{exp["expid"]} AMOC index')
    ax.set_ylabel('SV')
    ax.set_xlabel('')
    ax.set_ylim(-5,40)
    pl.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/amoc_ind.png')    
    pl.show()


def mkplots(exps,dsets):
    amoc=mkamoc(exps[0], dsets[0])

    plot_amoc(exps[0],amoc)
    plot_amoc_ind(exps[0],amoc)

def main(exps):
    dsets=geosdset.load_data(exps, plotname, defaults)
    mkplots(exps,dsets)
    geosdset.close(dsets)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    main(exps)
