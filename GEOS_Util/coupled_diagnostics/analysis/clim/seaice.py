#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
import scipy as sp
import matplotlib.pyplot as pl
from g5lib import mappers, plotters, field, g5dset

# Read data
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.CtlTripolar(exp,'geosgcm_seaice')

exp.aice=exp.ctl('AICE',dates=exp.dates).clim(12)
exp.aice.data=sp.ma.masked_values(exp.aice.data,5e15); exp.aice.data=sp.ma.masked_values(exp.aice.data,0.0)
exp.hice=exp.ctl('HICE',dates=exp.dates).clim(12)
exp.hice.data=sp.ma.masked_values(exp.hice.data,5e15); exp.hice.data=sp.ma.masked_values(exp.hice.data,0.0)
exp.ui=exp.ctl('UI',dates=exp.dates).clim(12); exp.ui.data.mask=exp.aice.data.mask
exp.vi=exp.ctl('VI',dates=exp.dates).clim(12); exp.vi.data.mask=exp.aice.data.mask
exp.uice=field.cmplx(exp.ui,exp.vi)

exp.nha=exp.aice(lats=(0,90))
exp.nhh=exp.hice(lats=(0,90))
exp.nhu=exp.uice(lats=(0,90))
exp.sha=exp.aice(lats=(-90,0))
exp.shh=exp.hice(lats=(-90,0))
exp.shu=exp.uice(lats=(-90,0))

ind=[0,1,2]
exp.nha.jfm=exp.nha.subset(tind=ind).mean(0); exp.nha.jfm.name+=' NH, JFM'
exp.nhh.jfm=exp.nhh.subset(tind=ind).mean(0); exp.nhh.jfm.name+=' NH, JFM'
exp.nhu.jfm=exp.nhu.subset(tind=ind).mean(0); exp.nhu.jfm.name+=' NH, JFM'
exp.sha.jfm=exp.sha.subset(tind=ind).mean(0); exp.sha.jfm.name+=' SH, JFM'
exp.shh.jfm=exp.shh.subset(tind=ind).mean(0); exp.shh.jfm.name+=' SH, JFM'
exp.shu.jfm=exp.shu.subset(tind=ind).mean(0); exp.shu.jfm.name+=' SH, JFM'

ind=[6,7,8]
exp.nha.jas=exp.nha.subset(tind=ind).mean(0); exp.nha.jas.name+=' NH, JAS'
exp.nhh.jas=exp.nhh.subset(tind=ind).mean(0); exp.nhh.jas.name+=' NH, JAS'
exp.nhu.jas=exp.nhu.subset(tind=ind).mean(0); exp.nhu.jas.name+=' NH, JAS'
exp.sha.jas=exp.sha.subset(tind=ind).mean(0); exp.sha.jas.name+=' SH, JAS'
exp.shh.jas=exp.shh.subset(tind=ind).mean(0); exp.shh.jas.name+=' SH, JAS'
exp.shu.jas=exp.shu.subset(tind=ind).mean(0); exp.shu.jas.name+=' NH, JAS'

# Plots

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

alevs=sp.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 0.99])
hlevs=sp.array([0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])

popts=dict(parallels=sp.arange(-90.,120.,15.),labels=[0,0,0,0])
mopts=dict(meridians=sp.arange(0.,420.,90.),labels=[1,1,0,1])
nhmap=mappers.BaseMapper(projection='npstere',lon_0=0,boundinglat=50, mopts=mopts, popts=popts)
shmap=mappers.BaseMapper(projection='spstere',lon_0=0,boundinglat=-50, mopts=mopts, popts=popts)


def plot_map(fig,F1,F2,map,levs):
    pp=plotters.GeoPlotter(map=map)
    pp.copts.update(levels=levs)
    pp.cbar_opts.update(shrink=0.5)

    pl.figure(fig); pl.clf()
    pp(F1)
    pp.map.fillcontinents()
    pp.map.drawcountries()
    pp.method=pp.map.quiver
    pp.copts={}
    Nq=int(sp.sqrt(F2.dims[-1]*F2.dims[-2])/50)
    qq=pp(F2,skip=Nq)
    pl.quiverkey(qq,0,-0.15,0.5,r'$0.5 \frac{m}{s}$')
    ax=pl.gca(); ax.set_title(F1.name)
    pl.tight_layout()
    pl.show()

plot_map(1,exp.nha.jfm,exp.nhu.jfm,nhmap,alevs)
pl.savefig(path+'/seaice_nha_jfm.png')

plot_map(2,exp.nha.jas,exp.nhu.jas,nhmap,alevs)
pl.savefig(path+'/seaice_nha_jas.png')

plot_map(3,exp.sha.jfm,exp.shu.jfm,shmap,alevs)
pl.savefig(path+'/seaice_sha_jfm.png')

plot_map(4,exp.sha.jas,exp.shu.jas,shmap,alevs)
pl.savefig(path+'/seaice_sha_jas.png')

plot_map(5,exp.nhh.jfm,exp.nhu.jfm,nhmap,hlevs)
pl.savefig(path+'/seaice_nhh_jfm.png')

plot_map(6,exp.nhh.jas,exp.nhu.jas,nhmap,hlevs)
pl.savefig(path+'/seaice_nhh_jas.png')

plot_map(7,exp.shh.jfm,exp.shu.jfm,shmap,hlevs)
pl.savefig(path+'/seaice_shh_jfm.png')

plot_map(8,exp.shh.jas,exp.shu.jas,shmap,hlevs)
pl.savefig(path+'/seaice_shh_jas.png')
