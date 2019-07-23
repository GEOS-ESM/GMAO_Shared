#!/usr/bin/env python

from netCDF4 import Dataset
import datetime
import numpy as np
import matplotlib
matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)
matplotlib.rc('lines', linewidth=3)
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.mlab as mlab
import matplotlib.cbook as cbook
from pylab import *
import numpy.ma as ma
from scipy.spatial import cKDTree
from importlib import import_module
import os
import sys
import glob
import subprocess
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
import data_utils

def computeMeanExtent(filename):
    print filename
    allfields=np.genfromtxt(filename)
    extent=allfields[1:-1,-2]
#print allfields.shape
#print extent
    extent=ma.masked_where(extent==-9999,extent)
    mextent=np.mean(extent)
    return mextent

def compute_extent(fname,area,mask): 
   ncfile = Dataset(fname, 'r', format='NETCDF4')
   aice=ncfile.variables['AICE'][0]
   vel=ncfile.variables['VEL'][0]
   tmask=ncfile.variables['TMASK'][0]
   ncfile.close()
   atemp=area*vel
   atemp1=area.copy()
   #atemp=ma.masked_where(tmask<0.5, atemp)
   atemp=ma.masked_where(mask<0.5, atemp)
   atemp=ma.masked_where(aice<1.e-3, atemp)
   #atemp1=ma.masked_where(tmask<0.5, atemp1)
   atemp1=ma.masked_where(mask<0.5, atemp1)
   atemp1=ma.masked_where(aice<1.e-3, atemp1)
   return np.sum(atemp)/np.sum(atemp1) 

def compute_extent_mon(EXPDIR, COLLECTION, year0, year1, mon, area, POLE):
    asum = np.zeros(year1-year0+1)
    asum[:] = -9999.0
    for year in range(year0, year1+1):
        filename=EXPID+'.'+COLLECTION+'.monthly.'+str(year)+mon+'.nc4'
        fname=EXPDIR+'/'+COLLECTION+'/'+filename
        print fname
        if os.path.isfile(fname): 
            asum[year-year0] = compute_extent(fname,area,POLE)
    asum = ma.masked_where(asum==-9999.0, asum)
    return np.mean(asum)  


def arctic_mask(mask, lon, lat, tmask):
   xx = mask.copy()
   lon[lon<0.0] = lon[lon<0.0] + 360.0 
   xx[np.logical_and(np.logical_and(lon<236., lon>103.), lat>68.0)] = 1.0 
   xx[np.logical_and(np.logical_not(np.logical_and(lon<236., lon>103.)), lat>79.0)] = 1.0 
   xx[tmask < 0.5] = 0.0 
   return xx

def maskout_coastal(mask, lon, lat, tmask):
   xx = mask.copy()
   lon[lon<0.0] = lon[lon<0.0] + 360.0 
   lon_in = lon[tmask<0.5]
   lat_in = lat[tmask<0.5]
   lon_out = lon[np.logical_and(tmask>0.5, mask>0.5)]
   lat_out = lat[np.logical_and(tmask>0.5, mask>0.5)]
   zout = mask[np.logical_and(tmask>0.5, mask>0.5)].copy().flatten()
   #lon_out = lon.copy()
   #lat_out = lat.copy()
   xs, ys, zs = lon_lat_to_cartesian(lon_in.flatten(), lat_in.flatten())
   xt, yt, zt = lon_lat_to_cartesian(lon_out.flatten(), lat_out.flatten())
   tree = cKDTree(zip(xs, ys, zs))
   rs=tree.query_ball_point(zip(xt, yt, zt), 150.0/6371.0)
   for i,r in enumerate(rs):
       if r:
          zout[i] = 0.0
   zout.shape = lon_out.shape
   xx[np.logical_and(tmask>0.5, mask>0.5)] = zout
   return xx

      
def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z
   

nsidc_total_extent = np.zeros(12)
#allfields=np.genfromtxt('N_04_area.txt', dtype=(int, int, string, string, float, float))

dims={360:'1', 720:'05', 1440:'025'}

start_year=None
end_year=None

try:
    start_year=int(sys.argv[2])
except:
    pass

try:
    end_year=int(sys.argv[3])
except:
    pass
#year = start_year

#our_years=np.arange(start_year,end_year+1,1) 
our_months=np.arange(1,13,1)
total_extent = np.zeros(12)

years    = mdates.YearLocator()   # every year
months   = mdates.MonthLocator(range(1,13), bymonthday=15, interval=1)  # every month
yearsFmt = mdates.DateFormatter('%Y')
monsFmt = mdates.DateFormatter('%b')

x_axis = np.zeros(12)
i=1990
for j in our_months:
	x_axis[j-1]=mdates.date2num(datetime.date(i,j,15)) 
print x_axis
fig = plt.figure(figsize=(8,6), facecolor='w')
ax = fig.add_subplot(111)

linespec=['o-','o--']

#EXPDIR='/discover/nobackup/bzhao/'+EXPID
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
    pngname = 'drift_speed_anncycle'
except ImportError:
    EXPDIR=sys.argv[1]
    HOMDIR=EXPDIR
    EXPID=EXPDIR.split('/')[-1]
    PLOT_PATH = './'
    pngname = EXPID+'_DRIFT_SPEED_ANNCYCLE'
COLLECTION='geosgcm_seaice'

aa=subprocess.check_output(['grep', 'OGCM_IM', HOMDIR+'/AGCM.rc'])
im=int(aa.rstrip().decode('utf-8').split(':')[-1])
aa=subprocess.check_output(['grep', 'OGCM_JM', HOMDIR+'/AGCM.rc'])
jm=int(aa.rstrip().decode('utf-8').split(':')[-1])

area1,area,tmask=data_utils.get_area_from_cice(dims[im])
lon,lat = data_utils.get_latlon_from_cice(dims[im])

mask = np.zeros((jm, im))

mask = arctic_mask(mask,lon,lat, tmask[0])

mask = maskout_coastal(mask, lon, lat, tmask[0])

'''
m = Basemap(projection='npstere',lon_0=0,boundinglat=45, resolution='l')
m.drawcoastlines()
m.fillcontinents()
m.drawcountries()

x, y =m(lon,lat)
m.pcolormesh(x,y,mask,vmin=0.0, vmax=1.0)
m.drawparallels(np.arange(-90.,120.,15.),labels=[0,0,0,0]) # draw parallels
m.drawmeridians(np.arange(0.,420.,30.),labels=[0,0,0,0]) # draw meridians
plt.colorbar(orientation='vertical',extend='both',shrink=0.8)
plt.show()
'''

if len(glob.glob(EXPDIR+'/'+COLLECTION+'/'+'*.monthly.clim.*')) > 0 and start_year is None:
   for mon in range(1,13,1):
       SEASON='M'+str(mon)
       if mon < 10:
           SEASON='M0'+str(mon)
       filename=EXPID+'.'+COLLECTION+'.monthly.clim.'+SEASON+'.nc4'
       fname=EXPDIR+'/'+COLLECTION+'/'+filename
       print fname
       total_extent[mon-1] = compute_extent(fname, area, mask) 
else:
   flist = glob.glob(EXPDIR+'/'+COLLECTION+'/'+'*.monthly.[0-9]*')
   flist.sort()  
   if start_year is None:
      start_year = int(flist[0].split('/')[-1][-10:-6])
   if end_year is None:
      end_year = int(flist[-1].split('/')[-1][-10:-6])
   for mon in range(1,13,1):
       SEASON=str(mon)
       if mon < 10:
           SEASON='0'+str(mon)
       total_extent[mon-1] = compute_extent_mon(EXPDIR, COLLECTION, start_year, end_year, SEASON, area, POLE) 

total_extent = total_extent*86400*1.e-3  # m/s -> km/day
print total_extent


l1=ax.plot(x_axis, total_extent, 'b'+linespec[0] )

ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monsFmt)
ax.xaxis.set_minor_locator(months)

datemin=datetime.date(1990,1,1)
datemax=datetime.date(1990,12,31)

ax.set_xlim(datemin, datemax)
ax.set_ylim(0.0, 14.0)

ax.format_xdata = mdates.DateFormatter('%m-%d')
'''
if len(EXPID) > 30: 
    pos = (0.35, 0.11)
elif len(EXPID) > 20:
    pos = (0.45, 0.11)
else:
    pos = (0.6, 0.11)
fig.legend((l1, l2), (EXPID, 'SSM/I'), pos, fontsize=15)
'''
ax.grid(True)
fig.autofmt_xdate(bottom=0.1, rotation=0, ha='center')
title('Sea Ice Drift Speed ',fontsize=25)
ax.set_ylabel(r'km/day',fontsize=20)
#fig.legend((l1, l2), ('GEOS5', 'ICESat'), 'upper right')
#ax.set_size('large')
#plt.show()
if start_year and end_year:
   pngname=pngname+'_'+str(start_year)+'-'+str(end_year)
plt.savefig(PLOT_PATH+'/'+pngname)




