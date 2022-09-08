#!/usr/bin/env python3

'''
A collection of functions for: prepare_SST_ICEC/
'''

import xesmf as xe
import xarray as xr
import numpy as np
import netCDF4 as nc4
import datetime as datetime

def gather_gmao2022_kpar(in_file):
  """
  Read binary KPAR file in GMAO operations path.
  This file is known to be a monthly climatology.

  Inputs:
  ------
  in_file:       string (file name to read)

  Returns:       
  ------
  clim_kpar:     array (rank-3) (12, lat, lon) <-- Monthly climatology.
  lat:           array (rank-1) (latitudes)
  lon:           array (rank-1) (longitudes)
  """

  ny, nx = [720, 1440] # dimensions
  dy, dx = [180./ny, 360./nx] # resolution

  # date-line edge, pole-edge regular lat-lon
  lat  =  np.arange(-90.+0.5*dy,90.,dy)
  lon  =  np.arange(-180.+0.5*dx,180.,dx)

  # binary data format
  # adapted from
  # /home/yvikhlia/nobackup/coupled/Forcings/MOM6/CF0090x6C_TM0360xTM0320/interp_kpar/interp_kpar.py
  dlist = []
  dlist.append(('head', 'f4', 16))
  dlist.append(('h1','i4'))
  dlist.append(('xx','f4',(ny,nx)))
  dlist.append(('h2','i4'))

  ds = np.memmap(in_file, dtype=dlist, mode='r') # read
  nRec = ds.shape[0] # number of records

  # this is what comes to reading the data sequentially
  clim_kpar = np.zeros((nRec,ny,nx))
  for i, x in enumerate(ds):
   print("reading record:\t{}".format(i))
   clim_kpar[i,:,:] = x['xx'][:]

#
# plt.figure(figsize=(16,12))
# for iRec in range(nRec):
#  plt.subplot(4,4,iRec+1)
#  plt.contourf(lon, lat, data[iRec,:,:], np.linspace(0.01, 0.26, 26), cmap=plt.cm.gist_ncar), plt.colorbar()
#  plt.title('%i'%(iRec))
#
#Notice from above plots that:
# record 12 is same as 0
# record 13 is same as 1
# plt.pcolormesh(lon, lat, data[0,:,:]-data[12,:,:], vmin=-0.01, vmax=0.01, cmap=plt.cm.bwr), plt.colorbar()  

  return clim_kpar[0:12,:,:], lat, lon
# ---------------------------

def regrid_data(in_data, in_lat, in_lon, out_ny, out_nx):
  """
  Regrids from one resolution to another using xesmf.
  Input and output grids are assumed to be regular.

  Inputs:
  ------
  in_data:       array (rank-3 or rank2: (time, lat, lon) or (lat, lon))
  in_lat:        array (rank-1) (input resolution latitudes)
  in_lon:        array (rank-1) (input resolution longitudes)
  out_ny:        output resolution along latitude
  out_nx:        output resolution along longitude

  Returns:       
  ------
  out_data:      array (same as `in_data`)
  out_lat:       array (rank-1) (output resolution latitudes)
  out_lon:       array (rank-1) (output resolution longitudes)
  """

  # define input grid as an xarray dataset
  ds_in = xr.Dataset(
    {"latitude":  (["latitude"],  in_lat),
     "longitude": (["longitude"], in_lon)
    })

  # do same for output grid
  dy, dx = [180./out_ny, 360./out_nx] # resolution

  # date-line edge, pole-edge regular lat-lon
  out_lat  =  np.arange(-90. +0.5*dy, 90.,  dy)
  out_lon  =  np.arange(-180.+0.5*dx, 180., dx)

  ds_out = xr.Dataset(
    {"latitude":  (["latitude"],  out_lat),
     "longitude": (["longitude"], out_lon)
    })

  # regridder
  print('Generating regridding weights...')
  regridder = xe.Regridder(ds_in, ds_out, "conservative", periodic=True)
  print('Done!')
  # save weights
  fn = regridder.to_netcdf()
  # to reuse weights
  # fn = 'conservative_720x1440_360x720.nc' # name of weights file
  # regridder = xe.Regridder(ds_in, ds_out, 'conservative', weights=fn) 

  print('Regridding data...')
  if in_data.ndim == 3: # input data is 3-dim: (time, lat, lon)
    out_data = np.zeros((in_data.shape[0], out_ny, out_nx))
    for nrec in range(in_data.shape[0]):
      out_data[nrec,:,:] = regridder(in_data[nrec,:,:])      
  else: # 2-d: (lat, lon)
    out_data = regridder(in_data)
  print('Done!')

# plt.figure(figsize=(16, 4))
# plt.subplot(121)
# plt.pcolormesh(lon, lat, data[0,:,:], vmin=0., vmax=0.25, cmap=plt.cm.gist_ncar), plt.colorbar()
# plt.title('input')

# plt.subplot(122)
# plt.pcolormesh(lon_out, lat_out, data_out, vmin=0., vmax=0.25, cmap=plt.cm.gist_ncar), plt.colorbar()
# plt.title('output')

  return out_data, out_lat, out_lon
# ---------------------------

def write_out(fName,\
              begin_date, time_units, time_calendar, time_in, lats_in, lons_in,\
              data_in, data_name, data_longName, data_units, data_FillValue, data_vMin, data_vMax):

  """
  Write out a file with specified arguments, one kind of `data` only, one record at a time.

  Inputs:
  ------
  fName:         string (file name to write)

  begin_date:    datetime (reference date for time in output)
  time_units:    string (minutes OR hours OR days OR months OR years)
  time_calendar: string (gregorian OR whatever!)
  time_in:       datetime object with dates (a bunch of dates)

  lats_in:       array (rank-1) (latitudes)
  lons_in:       array (rank-1) (longitudes)

  data_in:       array (rank-3) (dim: (time, lat, lon))
  data_name:     string (name of the data)
  data_lonName:  string (descriptive name of data)
  data_units:    string (SI Units)
  data_FillValue:float  (missing value)
  data_vMin:     float (valid minimum value)
  data_vMax:     float (valid maximum value)

  Returns:       
  ------
                 None; write out a file with above `fName`.
  """

  nlats = lats_in.shape[0]
  nlons = lons_in.shape[0]
  #
  
  ncid = nc4.Dataset(fName, mode='w', format='NETCDF4')

  # global attributes
  current_date = datetime.datetime.today()
  ncid.title       = 'ExtData_for_GEOS-ESM'
  ncid.institution = 'NASA Global Modeling and Assimilation Office'
  ncid.source      = 'https://github.com/GEOS-ESM/GMAO_Shared/tree/main/GEOS_Util/pre/prepare_ocnExtData/'
  ncid.references  = 'https://github.com/GEOS-ESM/GMAO_Shared/tree/main/GEOS_Util/pre/prepare_ocnExtData/README.md'
  ncid.history     = 'File created on ' + current_date.strftime("%B %d, %Y, %H:%M:%S")

  time = ncid.createDimension('time', None)
  lat  = ncid.createDimension('lat',  nlats)
  lon  = ncid.createDimension('lon',  nlons)

  times             = ncid.createVariable('time','f4', ('time',))
  times.long_name   = 'time'
  times.units       = time_units + begin_date.strftime(" since %Y-%m-%d %H:%M:%S")
  times.calendar    = time_calendar
  times[:]          = nc4.date2num(time_in, units=times.units, calendar=times.calendar)

  lats              = ncid.createVariable('lat', 'f8', ('lat',))
  lats.long_name    = 'latitude'
  lats.units        = 'degrees_north'
  lats[:]           = lats_in

  lons              = ncid.createVariable('lon', 'f8', ('lon',))
  lons.long_name    = 'longitude'
  lons.units        = 'degrees_east'
  lons[:]           = lons_in

  data              = ncid.createVariable(data_name,'f4',\
                      ('time','lat','lon',), fill_value=data_FillValue)
  data.long_name    = data_longName
  data.units        = data_units
  data.vmin         = data_vMin
  data.vmax         = data_vMax
  data[:]           = data_in

  ncid.close()
