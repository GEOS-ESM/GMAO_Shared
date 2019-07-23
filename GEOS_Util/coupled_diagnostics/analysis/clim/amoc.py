#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys, cPickle
import netCDF4 as nc
import matplotlib.pyplot as pl
import scipy as sp
from matplotlib import dates, ticker
from mpl_toolkits.basemap.cm import sstanom
import g5lib.bilinear_scale as bscale
import g5lib.plotters as ptrs
from g5lib import g5dset

exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.CtlMOM(exp)

f=nc.Dataset(exp.ctl.basin_mask)

mm={'tmask':2,'AMOC_MASK':1}
for key in mm:
    try:
        x=f.variables[key]; x.set_auto_maskandscale(True)
        bmask=sp.ma.masked_not_equal(x[:],mm[key]).mask
    except KeyError:
        pass
        
ty=exp.ctl.subset(iind=0)
gm=exp.ctl.subset(iind=0)

for tt,xx in enumerate(exp.ctl.time):
    var=exp.ctl.fromfile('ty_trans', tind=tt)
    mask=sp.logical_or(var.data.mask,bmask[sp.newaxis,sp.newaxis])
    var.data.mask=mask
    ty.data[tt,:,:,0]=var.data.sum(3)

    try:
        var=exp.ctl.fromfile('ty_trans_gm',tind=tt)
    except KeyError:
        print 'No ty_trans_gm is found.'
    var.data.mask=mask
    gm.data[tt,:,:,0]=var.data.sum(3)

km=ty.grid['lev'].size
for l in reversed(xrange(km-1)):
    ty.data[:,l,:,:]=ty.data[:,l:l+2,:,:].sum(1)

ty.data+=gm.data; ty.data*=-1

ty.name=exp.ctl.name+' AMOC index'

exp.amoc_ind=ty(lats=(40., 50), levs=(1000., 1500.)).ave(2).ave(1)
exp.amoc_ind26=ty(lats=(26.,), levs=(1000.,))
exp.amoc=ty(dates=exp.dates,lats=(-30., 60.)).clim(12).mean(0); exp.amoc.name=exp.ctl.name+' AMOC'

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

pl.figure(1); pl.clf()
p=exp.amoc_ind.d();
p.lopts.update(color='red'); p(exp.amoc_ind26)
p.style='--'
p(exp.amoc_ind26.time_mean(12))
p.lopts.update(color='blue')
p(exp.amoc_ind.time_mean(12))
ax=p.axis
ax.set_ylim(-5,40)
if(exp.amoc_ind.time.size > 2500):
    myloc=dates.YearLocator((exp.amoc_ind.time.size/1000)*10)
    ax.xaxis.set_major_locator(myloc)
ax.set_ylabel('Sv')
ax.legend(('40N-50N','26N'))
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/amoc_ind.png')

clevs=sp.arange(-20,21,2)

x=exp.amoc

matplotlib.scale.register_scale(bscale.BiLinearScale)
pl.figure(2); pl.clf();
p=ptrs.Plotter2d(copts=dict(levels=clevs, cmap=sstanom),
                 cbar_opts=dict(orientation='horizontal'))
p(x)
del p.copts['cmap']
p.method=pl.contour
p.copts.update(colors='black')
p(x)
ax=p.axis; ax.invert_yaxis();
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.yaxis.set_major_locator(ticker.FixedLocator((200,400,600,800,1000,2000,3000,4000,5000)))
ax.set_yscale('bilinear',threshold=1000.,multiplier=5.)
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/atl_moc.png')



