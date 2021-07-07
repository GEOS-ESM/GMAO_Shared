'''
This modules provides utilities for reading different GEOS collections. 
'''

import os
import yaml, cftime
import numpy as np
import xarray as xr

R=6378e3 # Earth radius

def _load_geos(exp,collection):
    ds=xr.open_mfdataset(f'{exp["data_path"]}/{collection}/{exp["expid"]}.{collection}.monthly.??????.nc4')

    IM=ds.dims['lon']
    JM=ds.dims['lat']

    dx=np.radians(np.ones((JM,IM))*360./IM)*R
    dx*=np.cos(np.radians(ds['lat'].values[:,np.newaxis]))
    
    dy=np.radians(np.ones((JM,IM))*180./(JM-1))*R
    dy[0,:]*=0.5
    dy[-1,:]*=0.5
    
    area=np.ones((JM,IM))*dy
    area[1:-1,:]*=0.5*(0.5*dx[0:-2,:] + dx[1:-1,:] + 0.5*dx[2:,:])
    area[0,:]*=0.5*(dx[0,:]+dx[1,:])
    area[-1,:]*=0.5*(dx[-2,:]+dx[-1,:])
        
    ds.update({'dx': (('lat','lon'), dx),
               'dy': (('lat','lon'), dy),
               'area': (('lat','lon'), area)})
    
    return ds
        
def _load_tripolar(exp, collection):
    ds=xr.open_mfdataset(f'{exp["data_path"]}/{collection}/{exp["expid"]}.{collection}.monthly.??????.nc4')
    ds=ds.set_coords(['LON','LAT'])
    return ds

def _load_mom(exp, collection):
    ds=xr.open_mfdataset(f'{exp["data_path"]}/MOM_Output/{collection}.e??????01_00z.nc',decode_times=False)
#    ds.coords['time']=cftime.num2date(ds['time'],'days since 0001-01-01')

# Replace below fields with fields from MAPL_Tripolar.nc
#    ds_static=xr.open_dataset(exp["ocean_static"])
#    ds.coords.update({'geolon': ds_static['geolon'],
#                      'geolat': ds_static['geolat'],
#                      'geolon_c': ds_static['geolon_c'],
#                      'geolat_c': ds_static['geolat_c'],
#                      'geolon_u': ds_static['geolon_u'],
#                      'geolat_u': ds_static['geolat_u'],
#                      'geolon_v': ds_static['geolon_v'],
#                      'geolat_v': ds_static['geolat_v']})
    
#    ds.update({'dx': (('yh','xh'), ds_static['dxt']),
#               'dy': (('yh','xh'), ds_static['dyt']),
#               'area': (('yh','xh'), ds_static['area_t'])})

    return ds 

def _get_loader(type):
    loaders={'GEOS': _load_geos,
             'GEOSTripolar': _load_tripolar,
             'MOM': _load_mom}
    return loaders[type]

def load_exps(exp_conf):
    '''
    Returns a list of name spaces with meta data for experiment we want to plot and 
    all experiments we want to compare to.
    '''
    
    exps=[]
    with open(exp_conf) as ff:
        exp=yaml.safe_load(ff)
        exp.setdefault('cmpexp',[])
        exp.setdefault('dates',[None,None])
        exps.append(exp)
        
    for expconf in exps[0]['cmpexp']:
        with open(expconf,'r') as ff:
            exp=yaml.safe_load(ff)
            exp.setdefault('dates',[None,None])
            exps.append(exp)

    # Make directory for plots if it does not exists
    try:
        os.makedirs(exps[0]['plot_path'])
    except OSError:
        pass

    return exps

def load_collection(exps, colname, type='GEOS'):
    '''
    Loads a collection 'colname' for all experiments in 'exps' list.
    Returns a list of xarray data sets.
    '''

    dsets=[]
    for exp in exps:
        dsets.append(_get_loader(type)(exp, colname))

    return dsets

def close(dsets):
    for ds in dsets:
        ds.close()
