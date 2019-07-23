#!/usr/bin/env python

from netCDF4 import Dataset
import datetime
import numpy as np
import matplotlib
matplotlib.rc('lines', linewidth=3)
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.mlab as mlab
import matplotlib.cbook as cbook
from pylab import *
import numpy.ma as ma
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

def compute_extent(fname,area,POLE): 
   ncfile = Dataset(fname, 'r', format='NETCDF4')
   aice=ncfile.variables['AICE'][:]
   hice=ncfile.variables['HICE'][:]
   tmask=ncfile.variables['TMASK'][:]
   LAT=ncfile.variables['LAT'][:]
   lat = LAT
   ncfile.close()
   atemp=area
   atemp=ma.masked_where(tmask[0,:,:]<0.5, atemp)
   #atemp=ma.masked_where(aice[0,:,:]>1.e10, atemp)
   atemp=ma.masked_where(aice[0,:,:]<0.15,  atemp)
   #atemp=ma.masked_where(hice[0,:,:]<0.06,  atemp)
   if POLE=='N':
       atemp=ma.masked_where(lat<0.0, atemp)
   if POLE=='S':
       atemp=ma.masked_where(lat>0.0, atemp)
   return np.sum(atemp) 

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
    pngname = 'extent_anncycle'
except ImportError:
    EXPDIR=sys.argv[1]
    HOMDIR=EXPDIR
    EXPID=EXPDIR.split('/')[-1]
    PLOT_PATH = './'
    pngname = EXPID+'_EXTENT_ANNCYCLE'
COLLECTION='geosgcm_seaice'

aa=subprocess.check_output(['grep', 'OGCM_IM', HOMDIR+'/AGCM.rc'])
im=int(aa.rstrip().decode('utf-8').split(':')[-1])

area1,area,tmask=data_utils.get_area_from_cice(dims[im])


for k,POLE in enumerate(['N', 'S']):
  for mon in range(1,13,1):
    SEASON=str(mon)
    if mon < 10:
      SEASON='0'+str(mon)
    fname=NOBACKUP+'/ObservationData/NSIDC/'+POLE+'_'+SEASON+'_area.txt'
    nsidc_total_extent[mon-1]=computeMeanExtent(fname)
  print nsidc_total_extent
  if len(glob.glob(EXPDIR+'/'+COLLECTION+'/'+'*.monthly.clim.*')) > 0 and start_year is None:
      for mon in range(1,13,1):
          SEASON='M'+str(mon)
          if mon < 10:
              SEASON='M0'+str(mon)
          filename=EXPID+'.'+COLLECTION+'.monthly.clim.'+SEASON+'.nc4'
          fname=EXPDIR+'/'+COLLECTION+'/'+filename
          print fname
          total_extent[mon-1] = compute_extent(fname, area, POLE) 
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

  total_extent = total_extent*1.e-6
  print total_extent

  departure=total_extent-nsidc_total_extent
  rms = np.sqrt(sum(departure*departure)/len(departure))
  print 'RMS deviations from SSMI is: ',rms


  l1,l2=ax.plot(x_axis, total_extent, 'b'+linespec[k], x_axis, nsidc_total_extent, 'r'+linespec[k])

  ax.xaxis.set_major_locator(months)
  ax.xaxis.set_major_formatter(monsFmt)
  ax.xaxis.set_minor_locator(months)

  datemin=datetime.date(1990,1,1)
  datemax=datetime.date(1990,12,31)

  ax.set_xlim(datemin, datemax)
  ax.set_ylim(0.0, 25.0)

  ax.format_xdata = mdates.DateFormatter('%m-%d')
  if len(EXPID) > 30: 
      pos = (0.35, 0.11)
  elif len(EXPID) > 20:
      pos = (0.45, 0.11)
  else:
      pos = (0.6, 0.11)
  fig.legend((l1, l2), (EXPID, 'SSM/I'), pos, fontsize=15)
  ax.grid(True)
fig.autofmt_xdate(bottom=0.1, rotation=0, ha='center')
title('Sea Ice Extent ',fontsize=25)
ax.set_ylabel(r'$10^6km^2$',fontsize=20)
#fig.legend((l1, l2), ('GEOS5', 'ICESat'), 'upper right')
#ax.set_size('large')
#plt.show()
if start_year and end_year:
   pngname=pngname+'_'+str(start_year)+'-'+str(end_year)
plt.savefig(PLOT_PATH+'/'+pngname)





