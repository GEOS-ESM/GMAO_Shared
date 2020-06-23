#!/bin/env python

import os
import scipy as sp
import matplotlib.pyplot as pl
from mpl_toolkits.basemap.cm import sstanom, s3pcpn_l
from matplotlib import dates
from g5lib import field

# Read validation data set
obs={}
path=os.environ['NOBACKUP']+'/verification/stress_mon_clim'
execfile(path+'/ctl.py')
obs['ctl']=ctl

tx=ctl.fromfile('taux',kind=0).clim(12); tx.data*=10
ty=ctl.fromfile('tauy',kind=0).clim(12); ty.data*=10

tx.shiftgrid(30.);
tx.grid['lon']=sp.where(tx.grid['lon']<29.,tx.grid['lon']+360,\
                          tx.grid['lon'])

ty.shiftgrid(30.);
ty.grid['lon']=sp.where(ty.grid['lon']<29.,ty.grid['lon']+360,\
                          ty.grid['lon'])

var=field.cmplx(tx,ty); var.name=ctl.name+' TAU'

ind=[0,1,11]; obs['djf']=var.subset(tind=ind).ave(0); obs['djf'].name+=', DJF'
ind=[5,6,7]; obs['jja']=var.subset(tind=ind).ave(0); obs['jja'].name+=', JJA'
obs['am']=var.ave(0); obs['am'].name+=', Annual Mean'

# Calculate equatorial profile
lonind=sp.logical_and(var.grid['lon'][0]>=130.0,var.grid['lon'][0]<=280.0)
latind=sp.logical_and(var.grid['lat'][:,0]>=-2.1,var.grid['lat'][:,0]<=2.0)
obs['eqprof']=obs['am'].subset(iind=lonind,jind=latind).ave(2)

# Equatorial Annual Cycle
obs['eqac']=var.subset(iind=lonind,jind=latind).ave(2)
obs['eqac'].data-=obs['eqac'].ave(0).data
obs['eqac'].name=var.name+', Eq. Annual Cycle'

# Plots
path=os.environ['NOBACKUP']+'/verification/stress_mon_clim/pics'

copts1={}
copts1['levels']=(0.,0.2,0.4,0.6,0.8,1.,1.5,2,2.5,3)
copts1['cmap']=s3pcpn_l

def plot_map(figure,F,copts):
    Nq=10
    x=field.absolute(F)
    pl.figure(figure); pl.clf()
    x.copts=copts
    x.plot_map()
    F.plot_quiver(Nq)
    pl.show()

# DJF
season='djf'
plot_map(1,obs[season],copts1)
pl.savefig(path+'/tau_'+season+'_qscat.png')

# JJA
season='jja'
plot_map(1,obs[season],copts1)
pl.savefig(path+'/tau_'+season+'_qscat.png')

# AM
season='am'
plot_map(1,obs[season],copts1)
pl.savefig(path+'/tau_'+season+'_qscat.png')

# Plot Equatorial Annual Cycle
pl.figure(2);pl.clf()
obs['eqac'].copts={'levels': sp.arange(-0.2,0.21,0.02),\
                   'cmap' : sstanom,\
                   'timefmt': dates.DateFormatter('%b')}
obs['eqac'].plot2d()
obs['eqac'].copts={'func': pl.contour,\
                   'colors': 'black',\
                   'levels': sp.arange(-0.2,0.21,0.04),\
                   'timefmt': dates.DateFormatter('%b')}
obs['eqac'].plot2d()
ax=pl.gca(); ax.yaxis.set_major_locator(dates.MonthLocator())
ax.set_title(obs['ctl'].name+' Eq. Annual cycle')
pl.grid(); pl.show()
pl.savefig(path+'/taux_eq_ac_qscat.png')
    
