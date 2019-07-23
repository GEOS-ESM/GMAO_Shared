#!/bin/env python

import os
import scipy as sp
import matplotlib.pyplot as pl
import my_utils as utl
import my_plots as mpl

print 'Running SSS...'

# Read Levitus SSS climatology
path=os.environ['NOBACKUP']+'/verification/levitus'
execfile(path+'/ctl.py')

obs={}
obs['name']='Levitus SSS'
obs['ctl']=ctl

obs['clim']=obs['ctl'].fromfile('salt',kind=0); obs['clim'].name=obs['name']
obs['clim'].shiftgrid(30.)
obs['clim'].grid['lon']=sp.where(obs['clim'].grid['lon']<29.,obs['clim'].grid['lon']+360.,obs['clim'].grid['lon'])

# Calculate DJF, JJA and annual mean
obs['djf']=obs['clim'].subset(tind=[0,1,11]).ave(0); obs['djf'].name=obs['djf'].name+', DJF'
obs['jja']=obs['clim'].subset(tind=[5,6,7]).ave(0); obs['jja'].name=obs['jja'].name+', JJA'
obs['am']=obs['clim'].subset(tind=slice(0,None)).ave(0)
obs['am'].name=obs['am'].name+', Annual Mean'


###################### Do plots #######################################################
path=os.environ['NOBACKUP']+'/verification/levitus/pics'

def plot_field(field, fig):
    pl.figure(fig); pl.clf()
    field.copts={'levels': sp.arange(32.0,38.1,0.4)}
    field.plot_map(); field.copts.clear()
    field.copts['func']=field.map.contour
    field.copts['levels']=sp.arange(32.0,38.1,0.8)
    field.copts['colors']='black'
    field.plot_map()
    mean,std=field.aave(ret_std=True)
    mpl.draw_stat((mean.data.squeeze(),std.data.squeeze()))
    pl.show()

# DJF
plot_field(obs['djf'],1)
pl.savefig(path+'/sss_lev_djf.png')

# JJA
plot_field(obs['jja'],2)
pl.savefig(path+'/sss_lev_jja.png')

# Annual Mean
plot_field(obs['am'],3)
pl.savefig(path+'/sss_lev_am.png')

print '...done'
