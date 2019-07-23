#! /usr/bin/env python

from netCDF4 import Dataset
import matplotlib
matplotlib.rc('xtick', labelsize=15) 
matplotlib.rc('ytick', labelsize=15) 
matplotlib.rc('lines', linewidth=3)
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import array
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap as bm
import matplotlib.dates as mdates
import glob
from importlib import import_module
import struct
import datetime
import time
import sys
import os
import subprocess
import re
host=os.environ['HOST']
if re.match(r"^pfe", host):
   sys.path.append('/home6/bzhao/python_utils')
elif re.match(r"^discover", host):
   sys.path.append('/home/bzhao/python_utils')
else:
   sys.path.append('/home6/bzhao/python_utils')
#import read_utils
import plot_utils
import data_utils
#import get_info
from pylab import *

def get_geos_monthly_aice(filename):
	ncfile = Dataset(filename, 'r', format='NETCDF4')
	aice=ncfile.variables['AICE'][0,:,:]
	ncfile.close()
	return aice

def get_field_ts(EXPDIR, EXPID, COLLECTION, FIXYEAR, year, Year, ind_i, ind_j, ind_k, fld, fldname):
   for mon in np.arange(1,13,1):
	month = str(mon)
	if mon < 10:
		month='0'+str(mon)
	filename=EXPID+'.'+COLLECTION+'.monthly.'+FIXYEAR+month+'.nc4'
	fname=EXPDIR+'/'+COLLECTION+'/'+Year+'/'+filename
        print fname
        if os.path.isfile(fname): 
		ncfile = Dataset(fname, 'r', format='NETCDF4')
		if ind_k >= 0:
			atemp=ncfile.variables[fldname][0,ind_k,:,:]
		else:
			atemp=ncfile.variables[fldname][0,:,:]
		fld[(year-1)*12+mon-1]=atemp[ind_j,ind_i]
		ncfile.close()
	else:
		fld[(year-1)*12+mon-1]=-9999.0

def get_extent_ts(EXPDIR, EXPID, COLLECTION, FIXYEAR, year, Year, start_year, POLE, ar, fld):
   for mon in np.arange(1,13,1):
	month = str(mon)
	if mon < 10:
		month='0'+str(mon)
	#filename=EXPID+'.'+COLLECTION+'.monthly.'+str(year)+month+'.nc4'
	filename=EXPID+'.'+COLLECTION+'.monthly.'+Year+month+'.nc4'
	fname=EXPDIR+'/'+COLLECTION+'/'+filename
        print fname
        if os.path.isfile(fname): 
		ncfile = Dataset(fname, 'r', format='NETCDF4')
		aice=ncfile.variables['AICE'][:]
		hice=ncfile.variables['HICE'][:]
		tmask=ncfile.variables['TMASK'][:]
		LAT=ncfile.variables['LAT'][:]
		lat = LAT[:,:]
		ncfile.close()
		atemp=ar
       		atemp=ma.masked_where(tmask[0,:,:]<0.5, atemp)
       		atemp=ma.masked_where(aice[0,:,:]>1.e10, atemp)
       		atemp=ma.masked_where(aice[0,:,:]<0.15,  atemp)
       		#atemp=ma.masked_where(hice[0,:,:]<0.06,  atemp)
       		if POLE=='N':
       			atemp=ma.masked_where(lat<0.0, atemp)
       		if POLE=='S':
       			atemp=ma.masked_where(lat>0.0, atemp)
		fld[(year-start_year)*12+mon-1] = np.sum(atemp) 
	else:
		fld[(year-start_year)*12+mon-1]=-9999.0

def get_volume_ts(EXPDIR, EXPID, COLLECTION, FIXYEAR, year, Year, start_year, POLE, ar, fld, ice):
   for mon in np.arange(1,13,1):
	month = str(mon)
	if mon < 10:
		month='0'+str(mon)
	#filename=EXPID+'.'+COLLECTION+'.monthly.'+str(year)+month+'.nc4'
	filename=EXPID+'.'+COLLECTION+'.monthly.'+Year+month+'.nc4'
	fname=EXPDIR+'/'+COLLECTION+'/'+filename
        print fname
        if os.path.isfile(fname): 
		ncfile = Dataset(fname, 'r', format='NETCDF4')
		aice=ncfile.variables['AICE'][0,:,:]
		if ice==1:
			hice=ncfile.variables['HICE'][0,:,:]
		else:
			hice=ncfile.variables['HSNO'][0,:,:]
		tmask=ncfile.variables['TMASK'][:]
		LAT=ncfile.variables['LAT'][:]
		lat = LAT[:,:]
		ncfile.close()
		atemp=ar*hice
		a1=atemp.copy()
		ar1=ar.copy()	
       		atemp=ma.masked_where(tmask[0,:,:]<0.5, atemp)
       		a1=ma.masked_where(tmask[0,:,:]<0.5, a1)
       		atemp=ma.masked_where(aice>1.e10, atemp)
       		a1=ma.masked_where(aice>1.e10, a1)
       		if POLE=='N':
       			atemp=ma.masked_where(lat<0.0, atemp)
       			a1=ma.masked_where(lat<0.0, a1)
       			ar1=ma.masked_where(lat<0.0, ar1)
       		if POLE=='S':
       			atemp=ma.masked_where(lat>0.0, atemp)
       			a1=ma.masked_where(lat>0.0, a1)
       			ar1=ma.masked_where(lat>0.0, ar1)
       		a1=ma.masked_where(hice<1.e-11, a1)
       		ar1=ma.masked_where(hice<1.e-11, ar1)
		fld[(year-start_year)*12+mon-1] = np.sum(atemp) 
		if ice==0:
			print np.sum(a1)/np.sum(ar1)
	else:
		fld[(year-start_year)*12+mon-1]=-9999.0

def get_year(year):
   if year < 10:
      Year = '000'+str(year)
   elif year < 100:
      Year = '00'+str(year)
   elif year < 1000:
      Year = '0'+str(year)
   else:
      Year = str(year)
   return Year
#fig_index=1
#fig = plt.figure(num=fig_index, figsize=(8,5), facecolor='w')
#fig = plt.figure(num=fig_index, figsize=(12,7), facecolor='w')

#ind_i = 500
#ind_j = 398

#ind_i = 527
#ind_j = 393


dims={360:'1', 720:'05', 1440:'025'}

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
    pngname = 'extent_ts'
except ImportError:
    EXPDIR=sys.argv[1]
    HOMDIR=EXPDIR
    EXPID=EXPDIR.split('/')[-1]
    POLE=sys.argv[2]
    PLOT_PATH = './'
    pngname = EXPID+'_EXTENT_Monthly'

COLLECTION='geosgcm_seaice'

files = glob.glob(EXPDIR+'/'+COLLECTION+'/*monthly.??????.nc4')
files.sort()

start_yr = int(files[0].split('/')[-1].split('.')[-2][:4])
end_yr =  int(files[-1].split('/')[-1].split('.')[-2][:4])+1
#start_yr = 1950
#end_yr = 1974
start_year = start_yr
end_year = end_yr

#nsidc_data=data_utils.nsidc_ice_data(POLE, 0, 1)
nsidc_data_n=data_utils.nsidc_ice_data('N', 0, 1, start_yr, end_year)
nsidc_data_s=data_utils.nsidc_ice_data('S', 0, 1, start_yr, end_year)
#print nsidc_data.shape

#get grid cell area
aa=subprocess.check_output(['grep', 'OGCM_IM', HOMDIR+'/AGCM.rc'])
im=int(aa.rstrip().decode('utf-8').split(':')[-1])
area1,area,tmask=data_utils.get_area_from_cice(dims[im])

piomass_data=data_utils.piomass_volume_data(start_yr, end_yr-1)
print piomass_data
piomass_data = piomass_data * 0.1

num_yrs = end_yr-start_yr
fld  = np.zeros(num_yrs*12)
fld1 = np.zeros(num_yrs*12)
fld2 = np.zeros(num_yrs*12)
fld3 = np.zeros(num_yrs*12)
fld4 = np.zeros(num_yrs*12)
fld5 = np.zeros(num_yrs*12)
FIXYEAR='1950'

POLE = 'N'
for year in range(start_yr, end_yr):
    Year = get_year(year)
    get_extent_ts(EXPDIR, EXPID, COLLECTION,  FIXYEAR, year, Year, start_yr, POLE, area,  fld)
    get_volume_ts(EXPDIR, EXPID, COLLECTION,  FIXYEAR, year, Year, start_yr, POLE, area1, fld1, 1)
    get_volume_ts(EXPDIR, EXPID, COLLECTION,  FIXYEAR, year, Year, start_yr, POLE, area1, fld2, 0)

POLE = 'S'
for year in range(start_yr, end_yr):
    Year = get_year(year)
    get_extent_ts(EXPDIR, EXPID, COLLECTION,  FIXYEAR, year, Year, start_yr, POLE, area,  fld3)
    get_volume_ts(EXPDIR, EXPID, COLLECTION,  FIXYEAR, year, Year, start_yr, POLE, area1, fld4, 1)
    get_volume_ts(EXPDIR, EXPID, COLLECTION,  FIXYEAR, year, Year, start_yr, POLE, area1, fld5, 0)


fld = ma.masked_where(fld==-9999.0, fld) 
fld1 = ma.masked_where(fld1==-9999.0, fld1) 
fld2 = ma.masked_where(fld2==-9999.0, fld2) 
fld3 = ma.masked_where(fld3==-9999.0, fld3) 
fld4 = ma.masked_where(fld4==-9999.0, fld4) 
fld5 = ma.masked_where(fld5==-9999.0, fld5) 
#fld2 = fld2 - 273.15

#fld3[:] = -9999.0
#fld3[-12:-1]=nsidc_data[0:11]
#fld3[-1] = nsidc_data[-1]
#fld3 = ma.masked_where(fld3==-9999.0, fld3) 

nsidc_data_n=np.tile(nsidc_data_n, num_yrs)
nsidc_data_s=np.tile(nsidc_data_s, num_yrs)
#fld3=nsidc_data

total_extent_n = fld*1.e-6
total_volume_n = fld1*1.e-13
total_volume_s_n = fld2*1.e-13
total_extent_s = fld3*1.e-6
total_volume_s = fld4*1.e-13
total_volume_s_s = fld5*1.e-13
#print total_extent
#print total_volume


years    = mdates.YearLocator(1)   # every year
months   = mdates.MonthLocator(range(1,13), bymonthday=1, interval=6)  # every month
#months   = mdates.YearLocator()
yearsFmt = mdates.DateFormatter('%Y')
monsFmt = mdates.DateFormatter('%Y %b')

#dstart=datetime.date(start_year,1,15)
#dend=datetime.date(end_year,12,15)
#timedelt = datetime.timedelta(days=30)
#x_axis=mdates.drange(dstart,dend,timedelt) 

our_years=np.arange(start_year,end_year,1)
our_months=np.arange(1,13,1)

datemin=datetime.date(start_year,1,1)
datemax=datetime.date(end_year-1,12,31)

x_axis = np.zeros(num_yrs*12)
for i in our_years:
        for j in our_months:
                x_axis[(i-start_year)*12+j-1]=mdates.date2num(datetime.date(i,j,15))

for i in our_years:
        for j in our_months:
		#if j==4 or j==5:
		if j==12:
			#print i, j, total_extent[(i-start_year)*12+j-1], total_volume[(i-start_year)*12+j-1] 
			print i, j, total_extent_n[(i-start_year)*12+j-1], nsidc_data_n[(i-start_year)*12+j-1] 

extent_min = 0.0
extent_max = 25.0
volume_min = 0.0
#	volume_max = 2.5
#	extent_min = 1.0
#	extent_max = 20.0
#	volume_min = 0.0
volume_max = 5.0
	#extent_min = 0.0
	#extent_max = 40.0
	#volume_min = 0.0
	#volume_max = 15.0



panelsx=3
panelsy=1
#print x_axis
fig = plt.figure(figsize=(16,10), facecolor='w')

ax = fig.add_subplot(panelsx,panelsy,1)
l1=ax.plot(x_axis, total_extent_n, 'b-',  label='NH')
l1=ax.plot(x_axis, total_extent_s, 'b--', label='SH')
l1=ax.plot(x_axis, nsidc_data_n, 'k-',  label='NSIDC-NH')
l1=ax.plot(x_axis, nsidc_data_s, 'k--', label='NSIDC-SH')
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
#ax.xaxis.set_minor_locator(months)
ax.set_xlim(datemin, datemax)
ax.set_ylim(extent_min, extent_max)
ax.legend(fontsize=17)
ax.grid()
ax.format_xdata = mdates.DateFormatter('%Y')
#ax.set_ylabel(r'$10^6km^2$',fontsize=20)
ax.set_ylabel(r'10**6 km**2',fontsize=20)
title(EXPID+' EXTENT',fontsize=25)


print x_axis.shape
#print piomass_data.shape

ax = fig.add_subplot(panelsx,panelsy,2)
l1=ax.plot(x_axis, total_volume_n, 'b-',label='NH')
l1=ax.plot(x_axis, total_volume_s, 'b--',label='SH')
#print x_axis.shape
#print piomass_data.shape
l2=ax.plot(x_axis, piomass_data, 'k-', label='PIOMASS')
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
#ax.xaxis.set_minor_locator(months)
ax.set_xlim(datemin, datemax)
ax.set_ylim(volume_min, volume_max)
#ax.legend((l2), ('PIOMASS'), 'upper right')
ax.legend(fontsize=17)
ax.grid()
ax.format_xdata = mdates.DateFormatter('%Y')
#ax.set_ylabel(r'$10^{13}m^3$',fontsize=20)
ax.set_ylabel(r'10**13 m**3',fontsize=20)
title(EXPID+' ICE VOLUME',fontsize=25)

ax = fig.add_subplot(panelsx,panelsy,3)
l1=ax.plot(x_axis, total_volume_s_n, 'b-')
l1=ax.plot(x_axis, total_volume_s_s, 'b--')
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
#ax.xaxis.set_minor_locator(months)
ax.set_xlim(datemin, datemax)
ax.grid()
#ax.set_ylim(volume_min, volume_max)
ax.format_xdata = mdates.DateFormatter('%Y')
#ax.set_ylabel(r'$10^{13}m^3$',fontsize=20)
ax.set_ylabel(r'10*13 m**3',fontsize=20)
title(EXPID+' SNOW VOLUME',fontsize=25)


fig.autofmt_xdate()
#fig.autofmt_xdate(bottom=0.1, rotation=0, ha='center')
plt.savefig(PLOT_PATH+'/'+pngname)
