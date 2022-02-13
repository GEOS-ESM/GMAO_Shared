#!/usr/bin/env python3

"""
Purpose:
 - Write out a climatological k-PAR (Photosynthetically Available Radiation) file
 - _Current_ input is from what is used in the NASA Global Modeling and Assimilation Office (GMAO) 
   [Weather Analysis and Prediction System](https://gmao.gsfc.nasa.gov/weather_prediction/) and 
   [MERRA-2 reanalysis](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/). This is a climatological field
   whose source was **not** documented, therefore unknow.
"""

from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import numpy as np

from utils import write_out, gather_gmao2022_kpar, regrid_data


# path to a file on NCCS discover (https://www.nccs.nasa.gov/) maintained by GMAO "operations group"
input_file = '/discover/nobackup/projects/gmao/share/dao_ops/fvInput/g5gcm/bcs/realtime/SST/1440x720/SEAWIFS_KPAR_mon_clim.1440x720'

# input resolution of above data (on regular lat-lon grid)
# output resolution = 0.125 <-- Ben, do we need a way to regrid? Or write it out at input resolution?
in_res, out_res = [0.25, 0.25]

# get the data and its lat-lon
[clim_kpar, lat, lon] = gather_gmao2022_kpar(input_file)

# write out (to nc4 format)
if in_res == out_res: # at same resolution as input.
  out_data = clim_kpar
  out_lat  = lat
  out_lon  = lon
else:
  ny_out, nx_out = [int(180./out_res), int(360./out_res)] # dimensions
  [out_data, out_lat, out_lon] = regrid_data(clim_kpar, lat, lon, ny_out, nx_out)
#

out_file_template = 'clim_kpar_GMAO2022'
begin_date = datetime(year=1, month=1, day=15) # climatological data, mid-month at 00UTC
out_frequency = 1   # ouput frequency in months
num_records   = 12  # monthly climatology
# one data record at a time
for n in range(num_records):
  dates = begin_date + n*relativedelta(months=out_frequency)
  print('writing out for:\t', dates)
  out_file = out_file_template + '_' + dates.strftime('%Y%m%d') + '.nc'
  write_out(out_file, begin_date, 'days', '365_day', dates, out_lat, out_lon,\
  out_data[n,:,:], 'kPAR', 'climatological photosynthetically available radiation', 'm-1', 1.e+15, 0., 1.)
