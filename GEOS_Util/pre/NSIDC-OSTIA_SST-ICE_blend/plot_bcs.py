#!/usr/bin/env python

'''
to plot the SST and SIC boundary conditions (binary files) that are used by the GCM

SA, Jul 2018, Aug 2020
'''
#-----------------------------------------------------------------

import  argparse

import  numpy                as      np
import  matplotlib.pyplot    as      plt

from    datetime             import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# read a binary file, uses f2py, see file for compile instructions... change as you may please!
from    read_ops_bcs         import read_bin_data
#-----------------------------------------------------------------

def file_names(data_path_, fndSST, date_in):

	data_path = data_path_+ fndSST + '/'

	sst_file_ = data_path + 'Ostia_sst_' + date_in.strftime('%Y%m%d') + '.bin'
	sic_file_ = data_path + 'Ostia_ice_' + date_in.strftime('%Y%m%d') + '.bin'
	mask_file_= data_path + 'mask'                                    + '.bin'

	return sst_file_, sic_file_, mask_file_
#-----------------------------------------------------------------

def main():

   comm_args       = parse_args()

   data_path       = comm_args['data_path']

   year_           = comm_args['year']
   month_          = comm_args['month']
   day_            = comm_args['day']
   whichFndSST     = comm_args['fndSST']

   date_ = datetime( year_, month_, day_, 00, 0, 0)	
   [sst_file, sic_file, mask_file] = file_names( data_path, whichFndSST, date_)

	# read the OPS SST & SIC. Both are of dimension(1440, 2880); see read_ops_bcs.pyf
   [date_ofSST, nlon, nlat, lon, lat, sst_] = read_bin_data(sst_file, date_.strftime('%Y%m%d')) # SST  is in K
   [date_ofSIC, nlon, nlat, lon, lat, sic_] = read_bin_data(sic_file, date_.strftime('%Y%m%d')) # SIC  is non-dimensional
   [date_ofMask,nlon, nlat, lon, lat, mask_]= read_bin_data(mask_file,date_.strftime('%Y%m%d')) # MASK is non-dimensional

   sst_[sst_ == -999.] = np.nan
   sic_[sic_ == -999.] = np.nan
   #-----------------------------------------------------------------

   print('Plotting for [%s] using\n[%s]\n[%s]'%(date_.strftime('%Y%m%d'), sst_file, sic_file))
	
	# Plot
   fig = plt.figure(figsize=(10,6))

	# plot sst
   ax1 = plt.subplot(221, projection=ccrs.PlateCarree())
#  ax1 = plt.subplot(221, projection=ccrs.PlateCarree( central_longitude=180))

#  ax1.set_global()
   ax1.coastlines()
#  ax1.stock_img()
   ax1.add_feature(cfeature.LAND,  facecolor='white')
#  ax1.add_feature(cfeature.LAKES, edgecolor='gray')

   im1 = ax1.pcolormesh(lon, lat, sst_, transform=ccrs.PlateCarree(),\
         cmap=plt.cm.jet)#, vmin=270., vmax=310.)

   plt.colorbar(im1, pad=0.01, shrink=0.75)
   plt.title(r'SST (deg K) for %s'%(date_.strftime('%Y%m%d')))
   plt.axis('off')

	# plot sic
   ax2 = plt.subplot(222, projection=ccrs.PlateCarree())
#  ax2.set_global()
   ax2.coastlines()
   ax2.add_feature(cfeature.LAND,  facecolor='white')
#  ax2.add_feature(cfeature.LAKES, edgecolor='gray')

   im2 = ax2.pcolormesh(lon, lat, sic_, transform=ccrs.PlateCarree(),\
         cmap=plt.cm.jet)#, vmin=0.0, vmax=1.0)

   plt.colorbar(im2, pad=0.01, shrink=0.75)
   plt.title(r'SIC for %s'%(date_.strftime('%Y%m%d')))
   plt.axis('off')
#  -----------------------------------------------------------------

#  plot mask
   ax3 = plt.subplot(223, projection=ccrs.PlateCarree())
#  ax3.set_extent([-180, 180, 20, 55], ccrs.PlateCarree())
   ax3.coastlines()
   ax3.add_feature(cfeature.LAND,  facecolor='white')
   ax3.add_feature(cfeature.LAKES, edgecolor='gray')

   im3 = ax3.pcolormesh(lon, lat, mask_, transform=ccrs.PlateCarree(),\
         cmap=plt.cm.jet)#, vmin=280., vmax=310.)

   plt.colorbar(im3, pad=0.01, shrink=0.75)
   plt.title(r'Mask (=0 or 1, latter over Great Lakes and Caspian Sea')
   plt.axis('off')
#  -----------------------------------------------------------------

   ax4 = plt.subplot(224, projection=ccrs.PlateCarree())
   ax4.set_extent([-180, 180, 10, 60], ccrs.PlateCarree())
   ax4.coastlines()
   ax4.add_feature(cfeature.LAND,  facecolor='white')
   ax4.add_feature(cfeature.LAKES, edgecolor='gray')

   im4 = ax4.pcolormesh(lon, lat, mask_, transform=ccrs.PlateCarree(),\
         cmap=plt.cm.jet)#, vmin=280., vmax=310.)

   plt.colorbar(im4, pad=0.01, shrink=0.75)
   plt.title(r' Zoom in mask')
   plt.axis('off')

#-----------------------------------------------------------------
   fName_fig = whichFndSST + '_bcs_SST_SIC_%s'%(date_.strftime('%Y%m%d'))
   plt.savefig(fName_fig + '.png', dpi=120)	
   plt.close('all')
#-----------------------------------------------------------------

def parse_args():
   p = argparse.ArgumentParser(description = \
       'Script to plot the SST and sea ice concentration used as lower boundary conditions by the GCM')

   p.add_argument('-data_path',   type=str, help= 'path to the data',  default='/discover/nobackup/sakella/lake_mask/')

   p.add_argument('-year',   type=int, help= 'year',  default=2006)
   p.add_argument('-month',  type=int, help= 'month', default=4)
   p.add_argument('-day',    type=int, help= 'day',   default=1)
   p.add_argument('-fndSST', type=str, help= 'which foundation (e.g., OSTIA1, OSTIA2)', default='OSTIA1')

   return vars(p.parse_args())
#-----------------------------------------------------------------

if __name__ == '__main__':

        main()
#-----------------------------------------------------------------
