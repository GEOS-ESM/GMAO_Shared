#!/usr/bin/env python3

'''
A collection for scripts for: prepare_SST_ICEC/
'''

import netCDF4 as nc4

def write_out(fout='extData.nc'):
  ncid = nc4.Dataset(output_file,mode='w',format='NETCDF4')
  ncid.close()
