import os
import cftime
import xarray as xr

__all_=['ds']

oceanval=os.environ.get('OCEANVAL',
                        '/discover/nobackup/projects/gmao/oceanval/verification')

# Time units are wrong in the nc file, so no decode times
ds=xr.open_dataset(f'{oceanval}/GPCP/GPCP_2017jul.nc', decode_times=False)
ds.coords['time']=cftime.num2date(ds['time'], 
                                  'months since 1979-01-01', 
                                  calendar='360_day')
