#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys, cPickle
from importlib import import_module
import matplotlib.pyplot as pl
from matplotlib import dates
import scipy as sp
import g5lib.plotters as ptrs
from g5lib import g5dset

exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.CtlTripolar(exp,'geosgcm_seaice')

# NH

A=exp.ctl('AREA')
Anh=A(lats=(0,90))
Ash=A(lats=(-90,0))

tmask=exp.ctl('TMASK')
aice=exp.ctl('AICE'); aice.data*=tmask.data
eice=aice.subset(); eice.data=sp.ma.where(aice.data<0.15,0.,1.)
hice=exp.ctl('HICE'); hice.data*=tmask.data

nha=aice(lats=(0,90)).aint(area=Anh)
nha.data/=1e6; nha.name=exp.ctl.name + ' NH Ice Area'
nhv=hice(lats=(0,90)).aint(area=Anh)
nhv.name=exp.ctl.name + ' NH Ice Volume'
nhh=nhv.subset(); nhh.data/=(nha.data*1e6);
nhh.name=exp.ctl.name + ' NH Ice Thickness'
nhe=eice(lats=(0,90)).aint(area=Anh)
nhe.data/=1e6; nhe.name=exp.ctl.name + ' NH Ice Extent'

sha=aice(lats=(-90,0)).aint(area=Ash)
sha.data/=1e6; sha.name=exp.ctl.name + ' SH Ice Area'
shv=hice(lats=(-90,0)).aint(area=Ash)
shv.name=exp.ctl.name + ' SH Ice Volume'
shh=shv.subset(); shh.data/=(sha.data*1e6);
shh.name=exp.ctl.name + ' SH Ice Thickness'
she=eice(lats=(-90,0)).aint(area=Ash)
she.data/=1e6; she.name=exp.ctl.name + ' SH Ice Extent'

Nt=exp.ctl.time.size
xx=nha.subset()
nsidc=import_module('nsidc')
piomas=import_module('piomas')

# Draw plots
path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

pnh=ptrs.Plotter1d(style='b-',lopts=dict(linewidth=1))
psh=ptrs.Plotter1d(style='g-',lopts=dict(linewidth=1))
pnhm=ptrs.Plotter1d(style='b-')
pshm=ptrs.Plotter1d(style='g-')
pobsnh=ptrs.Plotter1d(style='bo--',lopts=dict(linewidth=1))
pobssh=ptrs.Plotter1d(style='go--',lopts=dict(linewidth=1))

if(nha.time.size > 2500):
    myloc=dates.YearLocator((nha.time.size/1000)*10)

pl.figure(1); pl.clf()
pnh.axis=psh.axis=pnhm.axis=pshm.axis=pl.gca()
pnh(nha); psh(sha)
xx.data[:,0,0,0]=sp.tile(nsidc.anh_clim, int(Nt/12)+1)[:Nt]*1e6
pobsnh(xx)
xx.data[:,0,0,0]=sp.tile(nsidc.ash_clim, int(Nt/12)+1)[:Nt]*1e6
pobssh(xx)
pnhm(nha.time_mean(12)); pshm(sha.time_mean(12))
xx.data[:,0,0,0]=sp.tile(nsidc.anh_clim, int(Nt/12)+1)[:Nt]*1e6
pobsnh(xx)
xx.data[:,0,0,0]=sp.tile(nsidc.ash_clim, int(Nt/12)+1)[:Nt]*1e6
pobssh(xx)
ax=pl.gca(); ax.set_title(exp.ctl.name+' Ice Area'); ax.set_ylabel('km$^2$')
if(nha.time.size > 2500):
    ax.xaxis.set_major_locator(myloc)
ax.legend(('NH','SH','NSIDC clim'),loc='upper left')
pl.grid(); pl.tight_layout(); pl.show(); pl.savefig(path+'/aice.png')

pl.figure(2); pl.clf()
pnh.axis=psh.axis=pnhm.axis=pshm.axis=pl.gca()
pnh(nhv); psh(shv)
xx.data[:,0,0,0]=sp.tile(piomas.vnh_clim, int(Nt/12)+1)[:Nt]*1e12
pobsnh(xx)
pnhm(nhv.time_mean(12)); pshm(shv.time_mean(12))
ax=pl.gca(); ax.set_title(exp.ctl.name+' Ice Volume'); ax.set_ylabel('m$^3$')
if(nha.time.size > 2500):
    ax.xaxis.set_major_locator(myloc)
ax.legend(('NH','SH','PIOMAS clim'),loc='upper left')
pl.grid(); pl.tight_layout(); pl.show(); pl.savefig(path+'/volice.png')

pl.figure(3); pl.clf()
pnh.axis=psh.axis=pnhm.axis=pshm.axis=pl.gca()
pnh(nhh); psh(shh)
pnhm(nhh.time_mean(12)); pshm(shh.time_mean(12))
ax=pl.gca(); ax.set_title(exp.ctl.name+' Ice Thickness'); ax.set_ylabel('m')
if(nha.time.size > 2500):
    ax.xaxis.set_major_locator(myloc)
ax.legend(('NH','SH'),loc='upper left')
pl.grid(); pl.tight_layout(); pl.show(); pl.savefig(path+'/hice.png')

pl.figure(4); pl.clf()
pnh.axis=psh.axis=pnhm.axis=pshm.axis=pl.gca()
pnh(nhe); psh(she)
xx.data[:,0,0,0]=sp.tile(nsidc.enh_clim, int(Nt/12)+1)[:Nt]*1e6
pobsnh(xx)
xx.data[:,0,0,0]=sp.tile(nsidc.esh_clim, int(Nt/12)+1)[:Nt]*1e6
pobssh(xx)
pnhm(nhe.time_mean(12)); pshm(she.time_mean(12))
ax=pl.gca(); ax.set_title(exp.ctl.name+' Ice Extent'); ax.set_ylabel('km$^2$')
if(nhe.time.size > 2500):
    ax.xaxis.set_major_locator(myloc)
ax.legend(('NH','SH','NSIDC clim'),loc='upper left')
pl.grid(); pl.tight_layout(); pl.show(); pl.savefig(path+'/eice.png')
