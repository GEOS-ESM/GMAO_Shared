import os
import cftime
import xarray as xr

__all_=['ds_t','ds_s','ds_rho']

oceanval=os.environ.get('OCEANVAL',
                        '/discover/nobackup/projects/gmao/oceanval/verification')

mo=[str(i).zfill(2) for i in range(1,13)]

flist=[oceanval+'/WOA18/'+'woa18_decav_t'+mm+'_04.nc' for mm in mo]
ds_t=xr.open_mfdataset(flist,decode_times=False)
ds_t.coords['time']=cftime.num2date(ds_t.time,ds_t.time.units,calendar='360_day')

flist=[oceanval+'/WOA18/'+'woa18_decav_s'+mm+'_04.nc' for mm in mo]
ds_s=xr.open_mfdataset(flist,decode_times=False)
ds_s.coords['time']=cftime.num2date(ds_s.time,ds_s.time.units,calendar='360_day')
