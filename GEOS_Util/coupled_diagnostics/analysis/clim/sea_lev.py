#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
import scipy as sp
import matplotlib.pyplot as pl
from g5lib import plotters, g5dset

# Read data
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn2d')

varname='SLV'

exp.var=exp.ctl(varname,dates=exp.dates, levs=(0,)).clim(12)
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

clevs=sp.arange(-2.,2.1,0.5)

pp=plotters.GeoPlotter()
pp.copts.update(levels=clevs)

pl.figure(1); pl.clf()
pp(exp.var.djf)
pl.tight_layout()
pl.show()
pl.savefig(path+'/'+varname+'_djf.png')

pl.figure(2); pl.clf()
pp(exp.var.jja)
pl.tight_layout()
pl.show()
pl.savefig(path+'/'+varname+'_jja.png')

