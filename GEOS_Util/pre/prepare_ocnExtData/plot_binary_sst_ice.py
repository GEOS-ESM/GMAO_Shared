#!/usr/bin/env python

'''
to plot the SST and SIC binary formatted boundary conditions

SA, Jul 2018, Aug 2020, Feb 2022.
'''
#-----------------------------------------------------------------

import  argparse

import  numpy                as      np
import  matplotlib.pyplot    as      plt

from    datetime             import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# read a binary file, uses f2py, see file for compile instructions... change as you may please!
from read_sst_ice_bcs import read_bin_data
#-----------------------------------------------------------------

def file_names(data_path, date_in):

	sst_file = data_path + 'Ostia_sst_' + date_in.strftime('%Y%m%d') + '.bin'
	sic_file = data_path + 'Ostia_ice_' + date_in.strftime('%Y%m%d') + '.bin'

	return sst_file, sic_file
#-----------------------------------------------------------------

def main():

   comm_args       = parse_args()

   data_path       = comm_args['data_path']

   year_           = comm_args['year']
   month_          = comm_args['month']
   day_            = comm_args['day']

   date_ = datetime( year_, month_, day_, 00, 0, 0)
   [sst_file, sic_file] = file_names( data_path, date_)

	# read binary SST and SIC. Both are of dimension(1440, 2880)
   [date_ofSST, nlon, nlat, lon, lat, sst_] = read_bin_data(sst_file, date_.strftime('%Y%m%d')) # SST  is in K
   [date_ofSIC, nlon, nlat, lon, lat, sic_] = read_bin_data(sic_file, date_.strftime('%Y%m%d')) # SIC  is non-dimensional

   sst_[sst_ == -999.] = np.nan
   sic_[sic_ == -999.] = np.nan
   #-----------------------------------------------------------------

   print('Plotting for [%s] using\n[%s]\n[%s]'%(date_.strftime('%Y%m%d'), sst_file, sic_file))
	
	# Plot
   fig = plt.figure(figsize=(16,4))

	# plot sst
   ax1 = plt.subplot(121, projection=ccrs.PlateCarree())

   ax1.coastlines()
   ax1.add_feature(cfeature.LAND,  facecolor='white')
   ax1.add_feature(cfeature.LAKES, edgecolor='gray')

   im1 = ax1.pcolormesh(lon, lat, sst_, transform=ccrs.PlateCarree(),\
         cmap=plt.cm.jet)#, vmin=270., vmax=310.)

   plt.colorbar(im1, pad=0.01, shrink=0.75)
   plt.title(r'SST (deg K) for %s'%(date_.strftime('%Y%m%d')))
   plt.axis('off')

	# plot sic
   ax2 = plt.subplot(122, projection=ccrs.PlateCarree())
   ax2.coastlines()
   ax2.add_feature(cfeature.LAND,  facecolor='white')
   ax2.add_feature(cfeature.LAKES, edgecolor='gray')

   im2 = ax2.pcolormesh(lon, lat, sic_, transform=ccrs.PlateCarree(),\
         cmap=plt.cm.jet)#, vmin=0.0, vmax=1.0)

   plt.colorbar(im2, pad=0.01, shrink=0.75)
   plt.title(r'SIC for %s'%(date_.strftime('%Y%m%d')))
   plt.axis('off')
#  -----------------------------------------------------------------

   fName_fig = 'binary_SST_SIC_%s'%(date_.strftime('%Y%m%d'))
   plt.savefig(fName_fig + '.png', dpi=120)	
   plt.close('all')
#-----------------------------------------------------------------

def parse_args():
   p = argparse.ArgumentParser(description = \
       'Script to plot binary format SST and sea ice concentration used as lower boundary conditions by the GCM')

   p.add_argument('-data_path',   type=str, help= 'path to the data',  default='/discover/nobackup/sakella/BCs/test_ocean_ext/develop_15Feb2022/')

   p.add_argument('-year',   type=int, help= 'year',  default=2020)
   p.add_argument('-month',  type=int, help= 'month', default=11)
   p.add_argument('-day',    type=int, help= 'day',   default=28)

   return vars(p.parse_args())
#-----------------------------------------------------------------

if __name__ == '__main__':

        main()
#-----------------------------------------------------------------
