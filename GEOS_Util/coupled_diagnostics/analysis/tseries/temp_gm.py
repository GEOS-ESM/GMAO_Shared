#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
import scipy as sp
import matplotlib.pyplot as pl
from matplotlib import dates, ticker, colors
from mpl_toolkits.basemap.cm import sstanom
import g5lib.plotters as ptrs
from g5lib import cmaps as g5cmaps
import g5lib.bilinear_scale as bscale
from g5lib import g5dset

matplotlib.scale.register_scale(bscale.BiLinearScale)

exp=g5dset.read_exp(sys.argv[1])
varname='T'
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn3d')
exp.gm=exp.ctl.subset(iind=0,jind=0)

if exp.gm.grid['lev'][-1] < 0.0:
   exp.gm.grid['lev'][:]*=-1

thresh = 1000.0

if exp.ctl.time.size % 12 == 0:
   size = exp.ctl.time.size
else:
   size =  exp.ctl.time.size - exp.ctl.time.size % 12 

exp.gm = exp.gm.subset(tind=slice(0,size))

for tt in xrange(size):
    exp.gm.data[tt]=exp.ctl.fromfile(varname,tind=tt).aave().data[0]
exp.gam = exp.gm.time_mean(12)

exp.gam.data-=exp.gam.subset(tind=slice(0,1)).data

exp.upper = exp.gam(levs=(0,thresh))
exp.lower = exp.gam(levs=(thresh,6000))

# Plot

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

fig=pl.figure(1)
pl.clf()
clevs=sp.arange(-0.8,0.81,0.1)
n=colors.Normalize()
n.autoscale(clevs)
m=g5cmaps.FilledCmap(sstanom, fill_range=n((-0.1,0.1)))
exp.upper.name=exp.ctl.name +' '+varname+' Anomaly, Global Mean'
p=ptrs.Plotter2d(copts=dict(levels=clevs, cmap=m, norm=n), 
                # cbar_opts=dict(orientation='horizontal', pad=0.2),
                 cbar_opts=None,
                 axes=('time', 'lev'))
ax1 = pl.axes([0.1, 0.5, 0.8, 0.4])
#p(exp.gam)
p(exp.upper)
p.method=pl.contour
del p.copts['cmap']
p.copts.update(levels=clevs[0::2], colors='black')
#p(exp.gam)
p(exp.upper)
ax1.invert_yaxis(); 
ax1.set_xticks([])
ax1.set_xticklabels([])


exp.lower.name=''
p2=ptrs.Plotter2d(copts=dict(levels=clevs, cmap=m, norm=n), 
                #cbar_opts=dict(orientation='horizontal', pad=0.2),
                 cbar_opts=None,
                 axes=('time', 'lev'))
ax2 = pl.axes([0.1, 0.1, 0.8, 0.4])
#p(exp.gam)
p2(exp.lower)
#units=getattr(exp.lower,'units',"")
#coloraxis = [0.9, 0.1, 0.03, 0.8]
#cx = fig.add_axes(coloraxis, label=units, title=units)
#cb=pl.colorbar(orientation='vertical')
#cb=pl.colorbar(cax=cx,orientation='vertical',extend='both')
p2.method=pl.contour
del p2.copts['cmap']
del p2.copts['norm']
p2.copts.update(levels=clevs[0::2], colors='black')
p3=ptrs.Plotter2d(copts=dict(levels=clevs[0::2], colors='black'), 
                #cbar_opts=dict(orientation='horizontal', pad=0.2),
                 method=pl.contour,
                 cbar_opts=None,
                 axes=('time', 'lev'))
#p(exp.gam)
p3(exp.lower)

#ax=p.axis
for l in ax2.get_xticklabels():
    l.set_rotation(30)
ax2.xaxis.set_major_locator(dates.YearLocator())
#ax.yaxis.set_major_locator(ticker.FixedLocator((200,400,600,800,1000,1200,1400,1600,1800,2000,3000,4000,5000)))
ax2.invert_yaxis(); 
ax2.set_ylabel('depth, m')
ax2.yaxis.set_label_coords(-0.075,1.0)
ax2.set_xlabel('time, yr')
#cb.set_label(units)
#ax.set_yscale('bilinear',threshold=1000.,multiplier=5.)
#pl.grid(); 
pl.tight_layout()
pl.show()
pl.savefig(path+'/temp_gm.png')


