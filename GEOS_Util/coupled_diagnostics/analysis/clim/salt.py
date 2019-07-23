#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
from importlib import import_module
import scipy as sp
import matplotlib.pyplot as pl
from matplotlib import ticker, mlab, colors
from matplotlib.cm import jet
from mpl_toolkits.basemap.cm import sstanom
import g5lib.plotters as ptrs
import g5lib.domain as domain
from  g5lib import cmaps as g5cmaps
from g5lib import g5dset

varname='S'
# Read variable
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn3d')

if exp.ctl.grid['lev'][-1] < 0.0:
   exp.ctl.grid['lev'][:]*=-1 

dd=exp.ctl.domain
dd['dates']=exp.dates
ind=domain.Domain(**dd)(exp.ctl.grid, exp.ctl.time)
exp.am=exp.ctl.tave(varname,ind['tind']) 

exp.am.shiftgrid(30.)

exp.lat_depth=exp.am.ave(3)
exp.eq_depth=exp.am(lats=(-2.1,2.1)).ave(2)

exp.lat_depth.name=exp.am.name+', Annual Mean'
exp.eq_depth.name=exp.am.name+', Eq. Annual Mean'


# Read experiment to compare
exp1=g5dset.read_exp(exp.cmpexp)
exp1.ctl=g5dset.Ctl(exp1,'geosgcm_ocn3d')

if exp1.ctl.grid['lev'][-1] < 0.0:
   exp1.ctl.grid['lev'][:]*=-1 

dd=exp1.ctl.domain
dd['dates']=exp1.dates
ind=domain.Domain(**dd)(exp1.ctl.grid, exp1.ctl.time)
exp1.am=exp1.ctl.tave(varname, ind['tind']) 

exp1.am.shiftgrid(30.)
# If dimensions do not match, regrid
if exp1.am.dims[2:] != exp.am.dims[2:]:
    exp1.am.regrid(exp.am.grid)

exp1.lat_depth=exp1.am.ave(3)
exp1.eq_depth=exp1.am(lats=(0,))

# If levels do not match, interpolate
if exp1.lat_depth.dims[1] !=  exp.lat_depth.dims[1]:
    exp1.lat_depth.vinterp(exp.lat_depth.grid,newmask=exp.lat_depth.data.mask)
    exp1.eq_depth.vinterp(exp.eq_depth.grid,newmask=exp.eq_depth.data.mask)

# Read vaidation data and interpolate to exp grid
obs=import_module('levitus')

obs.am=obs.ctl('salt').ave(0)
obs.am.shiftgrid(30.)
obs.am.regrid(exp.am.grid)

obs.lat_depth=obs.am.ave(3)
obs.eq_depth=obs.am(lats=(0,))

obs.lat_depth.vinterp(exp.lat_depth.grid,newmask=exp.lat_depth.data.mask)
obs.eq_depth.vinterp(exp.eq_depth.grid,newmask=exp.eq_depth.data.mask)
###################### Do plots #######################################################

clevs=sp.arange(33.,36.1,0.2)

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass
    

def plot_field(field, fig, clevs, cmap, fill_range=None):
    pl.figure(fig)
    pl.clf()
    n=colors.Normalize()
    n.autoscale(clevs)
    if fill_range is not None:
        m=g5cmaps.FilledCmap(cmap, fill_range=n(fill_range))
    else: 
        m=cmap
    p=ptrs.Plotter2d(copts=dict(levels=clevs, cmap=m, norm=n))
    p(field)
    p.method=pl.contour
    p.copts=dict(levels=clevs[0::2], colors='black')
    p(field)
    ax=p.axis
    ax.set_ylabel('depth, m'); ax.invert_yaxis()
    return p

p=plot_field(exp.lat_depth, 1, clevs, jet)
ax=p.axis; ax.set_ylim(3000., 0.) 
ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/'+varname+'_lat_depth.png')

clevs1=sp.arange(-2,2.1,0.2)
frange=(-0.2,0.2)
dif=exp.lat_depth.subset(); dif.data-=exp1.lat_depth.data
dif.name=exp.ctl.name+'-'+exp1.ctl.name+' '+varname+', Annual Mean'
p=plot_field(dif, 2, clevs1, sstanom, frange)
ax=p.axis; ax.set_ylim(3000., 0.)
ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/'+varname+'_dif_lat_depth.png')


dif=exp.lat_depth.subset(); dif.data-=obs.lat_depth.data
dif.name=exp.ctl.name+'-'+obs.ctl.name+' '+varname+', Annual Mean'
p=plot_field(dif, 3, clevs1, sstanom, frange)
ax=p.axis; ax.set_ylim(3000., 0.)
ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/'+varname+'-obs_lat_depth.png')

clevs=sp.arange(34.,36.1,0.2)

p=plot_field(exp.eq_depth, 4, clevs, jet)
ax=p.axis; ax.set_ylim(500., 0.)
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/'+varname+'_eq_depth.png')

dif=exp.eq_depth.subset(); dif.data-=exp1.eq_depth.data
dif.name=exp.ctl.name+'-'+exp1.ctl.name+' '+varname+', Eq. Annual Mean'
p=plot_field(dif, 5, clevs1, sstanom, frange)
ax=p.axis; ax.set_ylim(500., 0.)
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/'+varname+'_dif_eq_depth.png')

dif=exp.eq_depth.subset(); dif.data-=obs.eq_depth.data
dif.name=exp.ctl.name+'-'+obs.ctl.name+' '+varname+', Eq. Annual Mean'
p=plot_field(dif, 6, clevs1, sstanom, frange)
ax=p.axis; ax.set_ylim(500., 0.)
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/'+varname+'-obs_eq_depth.png')

