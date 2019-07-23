#! /usr/bin/env python

from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import array
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
#import cmocean as cm
import glob
import struct
from importlib import import_module
import datetime
import time
import sys
import os
import re
host=os.environ['HOST']
if re.match(r"^pfe", host):
   sys.path.append('/home6/bzhao/python_utils')
   NOBACKUP='/nobackup/bzhao'
elif re.match(r"^discover", host):
   sys.path.append('/home/bzhao/python_utils')
   NOBACKUP='/discover/nobackup/bzhao'
else:
   sys.path.append('/home/bzhao/python_utils')
   NOBACKUP='/nobackup/bzhao'
import read_utils
import data_utils
import plot_utils
import math_utils
#import get_info
#from pylab import *


def plot_pole_new(LON,LAT,VAR,levels,setcolor,setnorm,titlestr,POLE,ptype,MER):

    if  (POLE=='N'):
        m = Basemap(projection='npstere',lon_0=0,boundinglat=65)
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
    plt.contour(x,y,VAR1, levels1, origin='lower',colors='magenta', linewidths=2)
    x, y =m(LON,LAT)
    plt.contourf(x,y,VAR, levels, origin='lower',cmap=setcolor, norm=setnorm, extend='both')

    m.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(np.arange(0.,420.,30.),labels=MER) # draw meridians


POLE='N'
fig_index=1
#cmp = cm.cm.ice
#fig = plt.figure(num=fig_index, figsize=(8,5), facecolor='w')
fig = plt.figure(num=fig_index, figsize=(14,14), facecolor='w')
if POLE=='N':
	fbot_levels = np.array([0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0,  5.0]) 
	#fbot_levels = np.arange(0, 3.75, 0.25) 
else:
	fbot_levels = np.array([0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0]) 
cmap2,norm=plot_utils.rescaled_cmap(fbot_levels)
aice_levels = np.array([0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 
                     0.8, 0.9, 0.95, 0.99]) 
cmap1,norm1=plot_utils.rescaled_cmap(aice_levels)
#cmap1.set_under('w')
#cmap2.set_under('w')

line_levs = np.array([0.15])

is_fm_yrs=['04', '05', '06', '07', '08']
is_on_yrs=['03','04', '05', '06', '07']

is_N = 19600
is_x = 140
is_y = 140


isfm=np.zeros((5,is_N,5))
ison=np.zeros((5,is_N,5))
is_dir=NOBACKUP+'/ObservationData/ICESat/'
for n in range(1,len(is_fm_yrs)+1,1):
	is_name=is_dir+'icesat_icethk_fm'+is_fm_yrs[n-1]+'_filled.dat'
	isfm[n-1]=np.loadtxt(is_name)
for n in range(1,len(is_on_yrs)+1,1):
	is_name=is_dir+'icesat_icethk_on'+is_on_yrs[n-1]+'_filled.dat'
	ison[n-1]=np.loadtxt(is_name)

#print isfm.shape, ison.shape
#print isfm[:,0,:]

isfmo = ma.masked_where(isfm==9999.0,isfm)
isono = ma.masked_where(ison==9999.0,ison)
isfmo = ma.masked_where(isfmo==-1.0,isfmo)
isono = ma.masked_where(isono==-1.0,isono)

#isfmo[isfmo==-1.0]=0.0
#isono[isono==-1.0]=0.0

#print isfmo.shape
#print isfmo[:,0,:]
isfm_m=np.mean(isfmo,axis=0)
ison_m=np.mean(isono,axis=0)

#print isfm_m.shape, ison_m.shape
#print isfm_m[0,:]
#SEASON='M08'
#YEAR='1973'
#YEAR='2010'

isfmg = np.reshape(isfm_m[:,-1], (is_x, is_y))
isfmlon = np.reshape(isfm_m[:,1], (is_x, is_y))
isfmlat = np.reshape(isfm_m[:,0], (is_x, is_y))

isfmg *= 0.01

isong = np.reshape(ison_m[:,-1], (is_x, is_y))
isonlon = np.reshape(ison_m[:,1], (is_x, is_y))
isonlat = np.reshape(ison_m[:,0], (is_x, is_y))

isong *= 0.01

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
    pngname = 'hice_icesat'
except ImportError:
    EXPDIR=sys.argv[1]
    HOMDIR=EXPDIR
    EXPID=EXPDIR.split('/')[-1]
    PLOT_PATH = './'
    pngname = EXPID+'_HICE_ICESAT'
COLLECTION='geosgcm_seaice'

#EXPDIR=sys.argv[1]
#EXPID=EXPDIR.split('/')[-1]
SEASON='M03'
fname=EXPDIR+'/'+COLLECTION+'/'+EXPID+'.'+COLLECTION+'.monthly.clim.'+SEASON+'.nc4' 
print fname
if os.path.isfile(fname):
   ncfile = Dataset(fname, 'r', format='NETCDF4')
   hi03=ncfile.variables['HICE'][0]
   LON=ncfile.variables['LON'][:]
   LAT=ncfile.variables['LAT'][:]
   lon = LON
   lat = LAT
   tmask=ncfile.variables['TMASK'][0]
   ncfile.close()

SEASON='M10'
fname=EXPDIR+'/'+COLLECTION+'/'+EXPID+'.'+COLLECTION+'.monthly.clim.'+SEASON+'.nc4' 
print fname
if os.path.isfile(fname):
   ncfile = Dataset(fname, 'r', format='NETCDF4')
   hi10=ncfile.variables['HICE'][0]
   ncfile.close()

SEASON='M11'
fname=EXPDIR+'/'+COLLECTION+'/'+EXPID+'.'+COLLECTION+'.monthly.clim.'+SEASON+'.nc4' 
print fname
if os.path.isfile(fname):
   ncfile = Dataset(fname, 'r', format='NETCDF4')
   hi11=ncfile.variables['HICE'][0]
   ncfile.close()


hifall=(hi10*31.0+hi11*30.0)/(31.0+30.0)
hispr=hi03

print LON.shape

hifall = ma.masked_where(tmask<0.5, hifall)
hispr = ma.masked_where(tmask<0.5, hispr)



titlestr=EXPID+' Feb-Mar'
#ax1 = plt.axes([-0.05, 0.225, 0.6, 0.6])
plt.subplot(2,2,1)
meridians=[1,0,1,1]
#plot_utils.plot_pole(lon,lat,aicem[0,:,:],aice_levels,'',POLE,'cont',meridians)
plot_pole_new(lon,lat,hispr,fbot_levels,cmap2,norm,'',POLE,'cont',meridians)
plt.title(titlestr,y=1.1,size=20)
#coloraxis = [0.05, 0.1, 0.4, 0.035]
#cx = fig.add_axes(coloraxis, label='m', title='1')
cbar=plt.colorbar(orientation='vertical',ticks=list(fbot_levels),extend='both',shrink=0.8)

titlestr=EXPID+' Oct-Nov'
#ax2 = plt.axes([0.425, 0.225, 0.6, 0.6])
plt.subplot(2,2,2)
meridians=[1,0,1,1]
#plot_utils.plot_pole(lon,lat,aicem[0,:,:],aice_levels,'',POLE,'cont',meridians)
plot_pole_new(lon,lat,hifall,fbot_levels,cmap2,norm,'',POLE,'cont',meridians)
plt.title(titlestr,y=1.1,size=20)

#plt.suptitle(EXPID,y=0.96,fontsize=25,fontweight='bold')
coloraxis = [0.5, 0.1, 0.5, 0.035]
#cx = fig.add_axes(coloraxis, label='m', title='m')
#cbar=plt.colorbar(cax=cx,orientation='horizontal',ticks=list(fbot_levels),extend='both')
cbar=plt.colorbar(orientation='vertical',ticks=list(fbot_levels),extend='both',shrink=0.8)

titlestr='ICESat Feb-Mar 2004-2008'
#ax1 = plt.axes([-0.05, 0.225, 0.6, 0.6])
plt.subplot(2,2,3)
meridians=[1,0,1,1]
#plot_utils.plot_pole(lon,lat,aicem[0,:,:],aice_levels,'',POLE,'cont',meridians)
plot_pole_new(isfmlon,isfmlat,isfmg,fbot_levels,cmap2,norm,'',POLE,'cont',meridians)
plt.title(titlestr,y=1.1,size=20)
#coloraxis = [0.05, 0.1, 0.4, 0.035]
#cx = fig.add_axes(coloraxis, label='m', title='1')
cbar=plt.colorbar(orientation='vertical',ticks=list(fbot_levels),extend='both',shrink=0.8)

titlestr='ICESat Oct-Nov 2003-2007'
#ax1 = plt.axes([-0.05, 0.225, 0.6, 0.6])
plt.subplot(2,2,4)
meridians=[1,0,1,1]
#plot_utils.plot_pole(lon,lat,aicem[0,:,:],aice_levels,'',POLE,'cont',meridians)
plot_pole_new(isonlon,isonlat,isong,fbot_levels,cmap2,norm,'',POLE,'cont',meridians)
plt.title(titlestr,y=1.1,size=20)
#coloraxis = [0.05, 0.1, 0.4, 0.035]
#cx = fig.add_axes(coloraxis, label='m', title='1')
cbar=plt.colorbar(orientation='vertical',ticks=list(fbot_levels),extend='both',shrink=0.8)
#plt.suptitle(EXPID,y=0.96,fontsize=16,fontweight='bold')
#pngname=EXPID+'_HICE_ICESat'
#print pngname
plt.savefig(PLOT_PATH+'/'+pngname)
#pcolor(lon,lat,aicem[0,:,:])
#colorbar()
