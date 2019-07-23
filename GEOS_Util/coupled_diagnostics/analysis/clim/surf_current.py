#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
import scipy as sp
import matplotlib.pyplot as pl
from mpl_toolkits.basemap.cm import s3pcpn_l
from g5lib import plotters, field, g5dset

# Read data
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn2d')

varname='Usurf'

exp.var=exp.ctl('US',dates=exp.dates).clim(12)
exp.var.data=exp.var.data+1j*exp.ctl('VS',dates=exp.dates).clim(12).data
exp.var.shiftgrid(30.)

ind=[0,1,2]
exp.var.djf=exp.var.subset(tind=ind).mean(0); exp.var.djf.name+=', DJF'

ind=[6,7,8]
exp.var.jja=exp.var.subset(tind=ind).mean(0); exp.var.jja.name+=', JJA'

# Plots

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

def plot_map(figure,F):
    Nq=F.dims[-2]/30
    copts=dict(levels=(0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.8,1.5), cmap=s3pcpn_l)
    p=plotters.GeoPlotter(copts=copts)
    pl.figure(figure); pl.clf()
    p(field.absolute(F))
    p.method=p.map.streamplot
    p.copts=dict(density=(10,5), color='black')
    p(F)
    pl.tight_layout()
    pl.show()

plot_map(1,exp.var.djf)
pl.savefig(path+'/'+varname+'_djf.png')

plot_map(2,exp.var.jja)
pl.savefig(path+'/'+varname+'_jja.png')

