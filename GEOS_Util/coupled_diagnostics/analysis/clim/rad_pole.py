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
from mpl_toolkits.basemap.cm import sstanom
import g5lib.plotters as ptrs
import g5lib.domain as domain
from  g5lib import cmaps as g5cmaps
from g5lib import g5dset
from g5lib import mappers, plotters



# Read vaidation data and interpolate to exp grid
#obs=import_module('COREII')
obs=import_module('JRA55do')

#print obs.swctl.time

# Read variable
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl(exp,'geosgcm_rad')

dd=exp.ctl.domain
dd['dates']=exp.dates
ind=domain.Domain(**dd)(exp.ctl.grid, exp.ctl.time)

varname='SWGDWN'
exp.sw=exp.ctl(varname,dates=exp.dates).clim(12)
exp.lwnet=exp.ctl('FLNS',dates=exp.dates).clim(12)
exp.lwup=exp.ctl('SFCEM',dates=exp.dates).clim(12)
exp.lw=exp.lwnet.subset()
exp.lw.data+=exp.lwup.data
exp.lw.name='LWGDWN'

ind=[5,6,7]
exp.sw.jja=exp.sw.subset(tind=ind).mean(0); exp.sw.jja.name+=' JJA'
exp.sw.jja.shiftgrid(30.)
exp.lw.jja=exp.lw.subset(tind=ind).mean(0); exp.lw.jja.name+=' JJA'
exp.lw.jja.shiftgrid(30.)

#sw=obs.rad('SWDN_MOD').subset(tind=ind).mean(0); sw.name+=' NH, JJA'
#lw=obs.rad('LWDN_MOD').subset(tind=ind).mean(0); lw.name+=' NH, JJA'

sw=obs.rad('rsds').subset(tind=ind).mean(0); sw.name+=' JJA'
lw=obs.rad('rlds').subset(tind=ind).mean(0); lw.name+=' JJA'

ind=[2,3,4]
exp.sw.mam=exp.sw.subset(tind=ind).mean(0); exp.sw.mam.name+=' MAM'
exp.sw.mam.shiftgrid(30.)
exp.lw.mam=exp.lw.subset(tind=ind).mean(0); exp.lw.mam.name+=' MAM'
exp.lw.mam.shiftgrid(30.)

#swmam=obs.rad('SWDN_MOD').subset(tind=ind).mean(0); swmam.name+=' NH, MAM'
#lwmam=obs.rad('LWDN_MOD').subset(tind=ind).mean(0); lwmam.name+=' NH, MAM'

swmam=obs.rad('rsds').subset(tind=ind).mean(0); swmam.name+=' MAM'
lwmam=obs.rad('rlds').subset(tind=ind).mean(0); lwmam.name+=' MAM'

ind=[0,1,11]
exp.sw.djf=exp.sw.subset(tind=ind).mean(0); exp.sw.djf.name+=' DJF'
exp.sw.djf.shiftgrid(30.)
exp.lw.djf=exp.lw.subset(tind=ind).mean(0); exp.lw.djf.name+=' DJF'
exp.lw.djf.shiftgrid(30.)

#swmam=obs.rad('SWDN_MOD').subset(tind=ind).mean(0); swmam.name+=' NH, MAM'
#lwmam=obs.rad('LWDN_MOD').subset(tind=ind).mean(0); lwmam.name+=' NH, MAM'

swdjf=obs.rad('rsds').subset(tind=ind).mean(0); swdjf.name+=' DJF'
lwdjf=obs.rad('rlds').subset(tind=ind).mean(0); lwdjf.name+=' DJF'

#obs.am=obs.ctl('temp').ave(0)
sw.shiftgrid(30.)
sw.regrid(exp.sw.jja.grid)
lw.shiftgrid(30.)
lw.regrid(exp.lw.jja.grid)

swmam.shiftgrid(30.)
swmam.regrid(exp.sw.mam.grid)
lwmam.shiftgrid(30.)
lwmam.regrid(exp.lw.mam.grid)

swdjf.shiftgrid(30.)
swdjf.regrid(exp.sw.djf.grid)
lwdjf.shiftgrid(30.)
lwdjf.regrid(exp.lw.djf.grid)
#print obs.jja.grid.dims


###################### Do plots #######################################################
# Plots

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

alevs=np.arange(0, 375, 25)
alevs1=np.arange(250, 360, 10)
hlevs=np.arange(-50,60, 10)
hlevs1=np.arange(-30,35, 5)

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

plot_map(1,exp.sw.jja,nhmap,alevs)
pl.savefig(path+'/swgdwn_nha_jja.png')

plot_map(2,sw,nhmap,alevs)
#pl.savefig(path+'/swgdwn_core2_nha_jja.png')
pl.savefig(path+'/swgdwn_jra55do_nha_jja.png')

dif=exp.sw.jja.subset(); dif.data-=sw.data
dif.name=exp.ctl.name+'-'+sw.name
plot_map(3,dif,nhmap,hlevs)
#pl.savefig(path+'/swgdwn_diff_core2_nha_jja.png')
pl.savefig(path+'/swgdwn_diff_jra55do_nha_jja.png')

plot_map(4,exp.lw.jja,nhmap,alevs1)
pl.savefig(path+'/lwgdwn_nha_jja.png')

plot_map(5,lw,nhmap,alevs1)
#pl.savefig(path+'/lwgdwn_core2_nha_jja.png')
pl.savefig(path+'/lwgdwn_jra55do_nha_jja.png')

dif=exp.lw.jja.subset(); dif.data-=lw.data
dif.name=exp.ctl.name+'-'+lw.name
plot_map(6,dif,nhmap,hlevs1)
#pl.savefig(path+'/lwgdwn_diff_core2_nha_jja.png')
pl.savefig(path+'/lwgdwn_diff_jra55do_nha_jja.png')

plot_map(1,exp.sw.mam,nhmap,alevs)
pl.savefig(path+'/swgdwn_nha_mam.png')

plot_map(2,swmam,nhmap,alevs)
#pl.savefig(path+'/swgdwn_core2_nha_mam.png')
pl.savefig(path+'/swgdwn_jra55do_nha_mam.png')

dif=exp.sw.mam.subset(); dif.data-=swmam.data
dif.name=exp.ctl.name+'-'+swmam.name
plot_map(3,dif,nhmap,hlevs)
#pl.savefig(path+'/swgdwn_diff_core2_nha_mam.png')
pl.savefig(path+'/swgdwn_diff_jra55do_nha_mam.png')

plot_map(4,exp.lw.mam,nhmap,alevs1)
pl.savefig(path+'/lwgdwn_nha_mam.png')

plot_map(5,lwmam,nhmap,alevs1)
#pl.savefig(path+'/lwgdwn_core2_nha_mam.png')
pl.savefig(path+'/lwgdwn_jra55do_nha_mam.png')

dif=exp.lw.mam.subset(); dif.data-=lwmam.data
dif.name=exp.ctl.name+'-'+lwmam.name
plot_map(6,dif,nhmap,hlevs1)
#pl.savefig(path+'/lwgdwn_diff_core2_nha_mam.png')
pl.savefig(path+'/lwgdwn_diff_jra55do_nha_mam.png')

plot_map(1,exp.sw.djf,shmap,alevs)
pl.savefig(path+'/swgdwn_sha_djf.png')

plot_map(2,swdjf,shmap,alevs)
#pl.savefig(path+'/swgdwn_core2_nha_mam.png')
pl.savefig(path+'/swgdwn_jra55do_sha_djf.png')

dif=exp.sw.djf.subset(); dif.data-=swdjf.data
dif.name=exp.ctl.name+'-'+swdjf.name
plot_map(3,dif,shmap,hlevs)
#pl.savefig(path+'/swgdwn_diff_core2_nha_mam.png')
pl.savefig(path+'/swgdwn_diff_jra55do_sha_djf.png')

plot_map(4,exp.lw.djf,shmap,alevs1)
pl.savefig(path+'/lwgdwn_sha_djf.png')

plot_map(5,lwdjf,shmap,alevs1)
#pl.savefig(path+'/lwgdwn_core2_nha_mam.png')
pl.savefig(path+'/lwgdwn_jra55do_sha_djf.png')

dif=exp.lw.djf.subset(); dif.data-=lwdjf.data
dif.name=exp.ctl.name+'-'+lwdjf.name
plot_map(6,dif,shmap,hlevs1)
#pl.savefig(path+'/lwgdwn_diff_core2_nha_mam.png')
pl.savefig(path+'/lwgdwn_diff_jra55do_sha_djf.png')

dif=exp.lw.djf.subset(); dif.data+=exp.sw.djf.subset().data; dif.data-=swdjf.data+lwdjf.data
dif.name=exp.ctl.name+'-JRA55do'
plot_map(7,dif,shmap,hlevs1)
#pl.savefig(path+'/lwgdwn_diff_core2_nha_mam.png')
pl.savefig(path+'/rad_diff_jra55do_sha_djf.png')
