#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import os, sys
from importlib import import_module
import scipy as sp
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import ticker, mlab, colors
from matplotlib.cm import jet
import matplotlib as mpl
mpl.rcParams['image.cmap'] = 'jet'
from mpl_toolkits.basemap.cm import sstanom
import g5lib.plotters as ptrs
import g5lib.domain as domain
from  g5lib import cmaps as g5cmaps
from g5lib import g5dset
from g5lib import mappers, plotters



# Read vaidation data and interpolate to exp grid
obs=import_module('WOA13')

#print obs.swctl.time

# Read variable
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn2d')

dd=exp.ctl.domain
dd['dates']=exp.dates
ind=domain.Domain(**dd)(exp.ctl.grid, exp.ctl.time)

varname='TS'
exp.sst=exp.ctl(varname,dates=exp.dates).clim(12)
exp.sst.name='SST'

ind=[0,1,11]
exp.sst.djf=exp.sst.subset(tind=ind).mean(0); exp.sst.djf.name+=' DJF'
exp.sst.djf.shiftgrid(30.)

#swmam=obs.rad('SWDN_MOD').subset(tind=ind).mean(0); swmam.name+=' NH, MAM'
#lwmam=obs.rad('LWDN_MOD').subset(tind=ind).mean(0); lwmam.name+=' NH, MAM'

sst=obs.woa('t_an').subset(tind=ind,kind=[1]).mean(0); sst.name+=' DJF'
sst.data += 273.15

#obs.am=obs.ctl('temp').ave(0)
sst.shiftgrid(30.)
sst.regrid(exp.sst.djf.grid)
#print obs.jja.grid.dims


###################### Do plots #######################################################
# Plots

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

alevs=np.arange(271, 281, 1)
hlevs=np.arange(-3, 3.5, 0.5)

popts=dict(parallels=sp.arange(-90.,120.,15.),labels=[0,0,0,0])
mopts=dict(meridians=sp.arange(0.,420.,90.),labels=[1,1,0,1])
nhmap=mappers.BaseMapper(projection='npstere',lon_0=0,boundinglat=45, mopts=mopts, popts=popts)
shmap=mappers.BaseMapper(projection='spstere',lon_0=0,boundinglat=-50, mopts=mopts, popts=popts)


def plot_map(fig,F1,map,levs):
    pp=plotters.GeoPlotter(map=map)
    pp.copts.update(levels=levs)
    pp.cbar_opts.update(shrink=0.5)

    pl.figure(fig); pl.clf()
    pp(F1)
    pp.map.fillcontinents()
    pp.map.drawcountries()
    ax=pl.gca(); ax.set_title(F1.name)
    pl.tight_layout()
    pl.show()



plot_map(1,exp.sst.djf,shmap,alevs)
pl.savefig(path+'/sst_sha_djf.png')

plot_map(2,sst,shmap,alevs)
pl.savefig(path+'/sst_woa13_sha_djf.png')

dif=exp.sst.djf.subset(); dif.data-=sst.data
dif.name=exp.ctl.name+'-'+sst.name
plot_map(3,dif,shmap,hlevs)
pl.savefig(path+'/sst_diff_woa13_sha_djf.png')

