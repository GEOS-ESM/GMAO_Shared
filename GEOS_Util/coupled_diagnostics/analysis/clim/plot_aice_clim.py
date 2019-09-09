#! /usr/bin/env python

from netCDF4 import Dataset
import matplotlib
matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)
matplotlib.rc('lines', linewidth=3)
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import array
#import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import cmocean as cm
import glob
import struct
import datetime
import time
from importlib import import_module
import sys
import os
import re
host=os.environ['HOST']
if re.match(r"^pfe", host):
   sys.path.append('/home6/bzhao/python_utils')
elif re.match(r"^discover", host):
   sys.path.append('/home/bzhao/python_utils')
else:
   sys.path.append('/home6/bzhao/python_utils')
import read_utils
import data_utils
import plot_utils
import math_utils
#import get_info
#from pylab import *


def plot_pole_new(LON,LAT,VAR,levels,setcolor,setnorm,titlestr,POLE,ptype,MER):

    if  (POLE=='N'):
        m = Basemap(projection='npstere',lon_0=0,boundinglat=45)
    if  (POLE=='S'):
        m = Basemap(projection='spstere',lon_0=180,boundinglat=-45)
    m.drawcoastlines()
    m.fillcontinents()
    m.drawcountries()
    plt.title(titlestr)
    x, y =m(LON,LAT)
    if ptype=='cont':
        m.contourf(x,y,VAR, levels, origin='lower',cmap=setcolor, norm=setnorm, extend='both')
        #plt.colorbar(orientation='vertical',extend='both',shrink=0.4)
    if ptype=='plot':
        plt.plot(x,y,'.')#marker='.',color='k')
    if ptype=='scatter':
        #plt.plot(x,y,'.')#marker='.',color='k')
        VAR[abs(VAR)>999]='nan'
        print 'min='+str(np.nanmin(VAR))
        print 'max='+str(np.nanmax(VAR))
        valmin=min(levels)
        valmax=max(levels)

        plt.scatter(x,y,8*VAR/VAR,VAR,marker='o',vmin=valmin,vmax=valmax,cmap='jet',linewidths=0)

    m.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(np.arange(0.,420.,30.),labels=MER) # draw meridians

def plot_pole_2(LON,LAT,VAR, LON1,LAT1,VAR1, levels, levels1, setcolor,setnorm,titlestr,POLE,ptype,MER):

    if  (POLE=='N'):
        m = Basemap(projection='npstere',lon_0=0,boundinglat=45)
    if  (POLE=='S'):
        m = Basemap(projection='spstere',lon_0=180,boundinglat=-45)
    m.drawcoastlines()
    m.fillcontinents()
    m.drawcountries()
    plt.title(titlestr)
    x, y =m(LON1,LAT1)
    cs=plt.contour(x,y,VAR1, levels1, origin='lower',colors='magenta', linewidths=2)
    lines, legs = cs.legend_elements()
    x, y =m(LON,LAT)
    plt.contourf(x,y,VAR, levels, origin='lower',cmap=setcolor, norm=setnorm, extend='both')

    m.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(np.arange(0.,420.,30.),labels=MER) # draw meridians
    return lines,legs

#POLE='N'
POLE=sys.argv[2]
fig_index=1
cmp = cm.cm.ice
#fig = plt.figure(num=fig_index, figsize=(8,5), facecolor='w')
fig = plt.figure(num=fig_index, figsize=(8,12), facecolor='w')
if POLE=='N':
	fbot_levels = np.array([0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0,  5.0]) 
	#fbot_levels = np.arange(0, 3.75, 0.25) 
else:
	fbot_levels = np.array([0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0]) 
cmap2,norm=plot_utils.rescaled_cmap(fbot_levels, cmap=cmp)
aice_levels = np.array([0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 
                     0.8, 0.9, 0.95, 0.99]) 
cmap1,norm1=plot_utils.rescaled_cmap(aice_levels, cmap=cmp)
#cmap1.set_under('w')
#cmap2.set_under('w')

line_levs = np.array([0.15])

#SEASON='M08'
#YEAR='1973'
#YEAR='2010'
if POLE == 'N':
  SEASON=['M03', 'M09']
  MONTH=['MAR', 'SEP']
else:
  SEASON=['M09', 'M02']
  MONTH=['SEP', 'FEB']


try:
    exp=import_module(sys.argv[1])
    EXPDIR=exp.data_path 
    HOMDIR=os.environ['HOMDIR'] 
    EXPID=exp.expid
    PLOT_PATH=exp.plot_path
    try:
       os.makedirs(PLOT_PATH)
    except OSError:
       pass
    pngname = 'aice_clim_'+POLE
except ImportError:
    EXPDIR=sys.argv[1]
    HOMDIR=EXPDIR
    EXPID=EXPDIR.split('/')[-1]
    PLOT_PATH = './'
    pngname = EXPID+'_AICE_CLIM_'+POLE
COLLECTION='geosgcm_seaice'
#EXPDIR=sys.argv[1]
#EXPID=EXPDIR.split('/')[-1]
i = 1
for sea,mon in zip(SEASON,MONTH):
  fname=EXPDIR+'/'+COLLECTION+'/'+EXPID+'.'+COLLECTION+'.monthly.clim.'+sea+'.nc4' 
  print fname
  if os.path.isfile(fname):
     ncfile = Dataset(fname, 'r', format='NETCDF4')
     fbot=ncfile.variables['HICE'][:]
     LON=ncfile.variables['LON'][:]
     LAT=ncfile.variables['LAT'][:]
     lon = LON
     lat = LAT
     aice=ncfile.variables['AICE'][:]
     tmask=ncfile.variables['TMASK'][:]
     ncfile.close()
  else:
     files = glob.glob(EXPDIR+'/'+COLLECTION+'/*monthly.????'+sea[-2:]+'.nc4')
     files.sort() 
     ncfile = Dataset(files[0], 'r', format='NETCDF4')
     LON=ncfile.variables['LON'][:]
     LAT=ncfile.variables['LAT'][:]
     lon = LON
     lat = LAT
     tmask=ncfile.variables['TMASK'][:]
     ncfile.close()
     aice=np.zeros((1, tmask.shape[1], tmask.shape[2]))
     fbot=np.zeros((1, tmask.shape[1], tmask.shape[2]))
     for f in files:
       ncfile = Dataset(f, 'r', format='NETCDF4')
       hi=ncfile.variables['HICE'][:]
       ai=ncfile.variables['AICE'][:]
       ncfile.close()
       fbot += hi
       aice += ai
     fbot /= float(len(files))
     aice /= float(len(files))
     


#print aice.shape
  print LON.shape
  aicem = ma.masked_where(tmask<0.5, aice)
#aicem = ma.masked_where(aice>=9999.0, aice)
  print aicem.max()
  idx=math_utils.find_max_index(aicem)
  print idx

#fbotm = ma.masked_where(fbot>=9999.0, fbot)
  fbotm = ma.masked_where(tmask<0.5, fbot)

  month = int(sea[-2:])
  hlon,hlat,hsic=data_utils.hadisst1_ice_data(month)


  titlestr=mon
#ax1 = plt.axes([-0.05, 0.225, 0.6, 0.6])
  ax=plt.subplot(2,1,i)
  meridians=[1,0,1,1]
#plot_utils.plot_pole(lon,lat,aicem[0,:,:],aice_levels,'',POLE,'cont',meridians)
#plot_pole_new(lon,lat,aicem[0,:,:],aice_levels,cmap1,norm1,'',POLE,'cont',meridians)
  artists, legends = plot_pole_2(lon,lat,aicem[0,:,:], hlon, hlat, hsic, aice_levels, line_levs, cmap1,norm1,'',POLE,'cont',meridians)
  plt.title(titlestr,y=1.05,size=23)
  legends = [u'HADISST-1']
  leg=ax.legend(artists, legends, loc=(0.0, 0.05), fancybox=True, framealpha=0.0, fontsize=18) 
  #leg.get_frame().set_alpha(0.0)
  for text in leg.get_texts():
    plt.setp(text, color = 'm')
  i += 1
coloraxis = [0.85, 0.1, 0.05, 0.8]
cx = fig.add_axes(coloraxis, label='', title='')
cbar=plt.colorbar(cax=cx,orientation='vertical',ticks=list(aice_levels),extend='both')

#pngname=EXPID+'_AICE_CLIM_'+POLE
print pngname
plt.savefig(PLOT_PATH+'/'+pngname)
#pcolor(lon,lat,aicem[0,:,:])
#colorbar()
