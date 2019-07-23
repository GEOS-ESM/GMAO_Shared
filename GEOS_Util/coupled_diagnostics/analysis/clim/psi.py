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

varname='PSI'

exp.var=exp.ctl(varname,dates=exp.dates, levs=(0,)).clim(12)
exp.var.shiftgrid(30.)

exp.var.am=exp.var.mean(0); exp.var.am.name+=', Annual Mean'

# Plots

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

#clevs=sp.arange(-100.,101.,10.0)

pp=plotters.GeoPlotter()
#pp.copts.update(levels=clevs)

pl.figure(1); pl.clf()
pp(exp.var.am)
pl.show()
pl.savefig(path+'/'+varname+'.png')

