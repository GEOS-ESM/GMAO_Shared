#!/usr/bin/env python

'''
to plot the SST and SIC boundary conditions

SA, Jul 2018, Aug 2020, Feb 2022.
'''
#-----------------------------------------------------------------

import argparse

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
#-----------------------------------------------------------------

def main():

   comm_args       = parse_args()

   data_path       = comm_args['data_path']

   year_           = comm_args['year']
   month_          = comm_args['month']
   day_            = comm_args['day']

   date_ = datetime( year_, month_, day_, 00, 0, 0)
   file_name = data_path + 'sst_fraci_' + date_.strftime('%Y%m%d') + '.nc4'

	# read SST and ICEC
   ds = xr.open_dataset(file_name)
   sst=ds.SST.squeeze('time')
   ice=ds.FRACI.squeeze('time')
   #-----------------------------------------------------------------

   print('Plotting for [%s] using\n[%s]'%(date_.strftime('%Y%m%d'), file_name))
	
	# Plot
   fig = plt.figure(figsize=(16,4))

	# plot sst
   ax1 = plt.subplot(121, projection=ccrs.PlateCarree( central_longitude=0))

   ax1.coastlines()
   ax1.add_feature(cfeature.LAND,  facecolor='white', zorder=100)
   ax1.add_feature(cfeature.LAKES, edgecolor='gray',  zorder=50)

   im1 = ax1.pcolormesh(ds['lon'], ds['lat'],sst,\
         transform=ccrs.PlateCarree(), cmap=plt.cm.jet)#, vmin=270., vmax=310.)

   gl = ax1.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='0.85', alpha=0.5, linestyle='-', draw_labels=True)
   gl.top_labels = False
   cbar=fig.colorbar(im1, extend='both', shrink=0.5, ax=ax1)
   cbar.set_label(r'SST')
   ax1.set_title(r'[mean=%.2f, sdev=%2.f]'%(np.nanmean(sst, dtype=np.float64), np.nanstd(sst, dtype=np.float64)))

	# plot sic
   ax2 = plt.subplot(122, projection=ccrs.PlateCarree())
   ax2.coastlines()
   ax2.add_feature(cfeature.LAND,  facecolor='white', zorder=100)
   ax2.add_feature(cfeature.LAKES, edgecolor='gray',  zorder=50)

   im2 = ax2.pcolormesh(ds['lon'], ds['lat'],ice,\
         transform=ccrs.PlateCarree(), cmap=plt.cm.jet)#, vmin=0., vmax=1.)

   gl = ax2.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='0.85', alpha=0.5, linestyle='-', draw_labels=True)
   gl.top_labels = False
   cbar=fig.colorbar(im2, extend='both', shrink=0.5, ax=ax2)
   cbar.set_label(r'FRACI')
   ax2.set_title(r'[mean=%.2f, sdev=%2.f]'%(np.nanmean(ice, dtype=np.float64), np.nanstd(ice, dtype=np.float64)))
#  -----------------------------------------------------------------

   fName_fig = 'SST_SIC_%s'%(date_.strftime('%Y%m%d'))
   plt.savefig(fName_fig + '.png', dpi=120)	
   plt.close('all')
#-----------------------------------------------------------------

def parse_args():
   p = argparse.ArgumentParser(description = \
       'A script to plot SST and sea ice concentration used as lower boundary conditions by the GEOS AGCM')

   p.add_argument('-data_path',   type=str, help= 'path to the data',  default='/discover/nobackup/sakella/BCs/test_ocean_ext/prepare_ocnExtData/')

   p.add_argument('-year',   type=int, help= 'year',  default=2020)
   p.add_argument('-month',  type=int, help= 'month', default=11)
   p.add_argument('-day',    type=int, help= 'day',   default=28)

   return vars(p.parse_args())
#-----------------------------------------------------------------

if __name__ == '__main__':

        main()
