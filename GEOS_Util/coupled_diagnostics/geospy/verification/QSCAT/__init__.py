import os
import xarray as xr

__all_=['ds']

oceanval=os.environ.get('OCEANVAL',
                        '/discover/nobackup/projects/gmao/oceanval/verification')

ds=xr.open_dataset(f'{oceanval}/QSCAT/qscat_clim.nc')
