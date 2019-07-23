#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
from importlib import import_module
import netCDF4 as nc
import scipy as sp
import matplotlib.pyplot as pl
from matplotlib import colors
from mpl_toolkits.basemap.cm import s3pcpn_l, sstanom
import g5lib.plotters as ptrs
from  g5lib import cmaps as g5cmaps
from g5lib import g5dset

def read_precip(ctl, varname, dates):
    

    tprec=ctl(varname,dates=dates, levs=(1000,))
    tprec.data*=3600*24
            
    return tprec

# Read data
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl(exp,'geosgcm_surf')

exp.clim=read_precip(exp.ctl, 'TPREC', dates=exp.dates).clim(12)

exp.clim.shiftgrid(30.);

ind=[0,1,11]; exp.djf=exp.clim.subset(tind=ind).mean(0); exp.djf.name+=', DJF'
ind=[5,6,7]; exp.jja=exp.clim.subset(tind=ind).mean(0); exp.jja.name+=', JJA'
exp.am=exp.clim.mean(0); exp.am.name+=', Annual Mean'

# Read experiment to compare
exp1=g5dset.read_exp(exp.cmpexp)
exp1.ctl=g5dset.Ctl(exp1,'geosgcm_surf')

exp1.clim=read_precip(exp1.ctl, 'TPREC', dates=exp1.dates).clim(12)

exp1.clim.shiftgrid(30.);
if exp1.clim.dims[2:]!=exp.clim.dims[2:]:
    exp1.clim.regrid(exp.clim.grid, newmask=exp.clim.data.mask)

ind=[0,1,11]; exp1.djf=exp1.clim.subset(tind=ind).mean(0); exp1.djf.name+=' DJF'
ind=[5,6,7]; exp1.jja=exp1.clim.subset(tind=ind).mean(0); exp1.jja.name+=' JJA'
exp1.am=exp1.clim.mean(0); exp1.am.name+=' Annual Mean'

# Read validation data set
obs=import_module('gpcp')

obs.clim=obs.ctl('precip').clim(12)
obs.clim.shiftgrid(30.);

# Regrid to model grid
obs.clim.regrid(exp.clim.grid, newmask=exp.clim.data.mask)

ind=[0,1,11]; obs.djf=obs.clim.subset(tind=ind).mean(0)
obs.djf.name=obs.ctl.name+' '+obs.djf.name+' DJF'
ind=[5,6,7]; obs.jja=obs.clim.subset(tind=ind).mean(0)
obs.jja.name=obs.ctl.name+' '+obs.jja.name+' JJA'
obs.am=obs.clim.mean(0); obs.am.name=obs.ctl.name+' '+obs.am.name+' Annual Mean'

# Plots

def plot_field(field, fig, clevs, cmap, fill_range=None):
    pl.figure(fig); pl.clf()
    n=colors.Normalize()
    n.autoscale(clevs)
    if fill_range is not None:
        m=g5cmaps.FilledCmap(cmap, fill_range=n(fill_range))
    else: 
        m=cmap
    p=ptrs.GeoPlotter(copts=dict(levels=clevs, cmap=m, norm=n))
    p(field, stat=True)
    pl.show()
    
clevs1=(.2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.3, 2.6, 3, \
                 3.3, 3.6, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 11, 13, 15)
clevs2=(-8, -7, -6, -5, -4, -3, -2, -1, -.5, 0., .5, 1, 2, 3, 4, 5, 6, 7, 8)
frange=(-0.5,0.5)

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass
    
#DJF
plot_field(exp.djf,1,clevs1,s3pcpn_l)
pl.savefig(path+'/precip_djf.png')

dif=exp.djf.subset(); dif.data-=obs.djf.data
dif.name=exp.ctl.name+'-'+obs.ctl.name+' DJF'
plot_field(dif,2,clevs2,sstanom, frange)
pl.savefig(path+'/precip-obs_djf.png')

dif=exp.djf.subset(); dif.data-=exp1.djf.data
dif.name=exp.ctl.name+'-'+exp1.ctl.name+' DJF'
plot_field(dif,3,clevs2,sstanom, frange)
pl.savefig(path+'/precip-dif_djf.png')

#JJA
plot_field(exp.jja,1,clevs1,s3pcpn_l)
pl.savefig(path+'/precip_jja.png')

dif=exp.jja.subset(); dif.data-=obs.jja.data
dif.name=exp.ctl.name+'-'+obs.ctl.name+' JJA'
plot_field(dif,2,clevs2,sstanom, frange)
pl.savefig(path+'/precip-obs_jja.png')

dif=exp.jja.subset(); dif.data-=exp1.jja.data
dif.name=exp.ctl.name+'-'+exp1.ctl.name+' JJA'
plot_field(dif,3,clevs2,sstanom, frange)
pl.savefig(path+'/precip-dif_jja.png')

#AM
plot_field(exp.am,1,clevs1,s3pcpn_l)
pl.savefig(path+'/precip_am.png')

dif=exp.am.subset(); dif.data-=obs.am.data
dif.name=exp.ctl.name+'-'+obs.ctl.name+' Annual Mean'
plot_field(dif,2,clevs2,sstanom, frange)
pl.savefig(path+'/precip-obs_am.png')

dif=exp.am.subset(); dif.data-=exp1.am.data
dif.name=exp.ctl.name+'-'+exp1.ctl.name+' Annual Mean'
plot_field(dif,3,clevs2,sstanom, frange)
pl.savefig(path+'/precip-dif_am.png')



