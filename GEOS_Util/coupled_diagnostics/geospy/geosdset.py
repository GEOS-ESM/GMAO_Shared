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

    ds=ds.rename_dims({'lon': 'x', 'lat': 'y'})

    IM=ds.dims['x']
    JM=ds.dims['y']

    dx=np.radians(np.ones((JM,IM))*360./IM)*R
    dx*=np.cos(np.radians(ds['lat'].values[:,np.newaxis]))
    
    dy=np.radians(np.ones((JM,IM))*180./(JM-1))*R
    dy[0,:]*=0.5
    dy[-1,:]*=0.5
    
    area=np.ones((JM,IM))*dy
    area[1:-1,:]*=0.5*(0.5*dx[0:-2,:] + dx[1:-1,:] + 0.5*dx[2:,:])
    area[0,:]*=0.5*(dx[0,:]+dx[1,:])
    area[-1,:]*=0.5*(dx[-2,:]+dx[-1,:])

    ds.update({'dx': (('y','x'), dx),
               'dy': (('y','x'), dy),
               'mask': (('y','x'), ds['MASKO'][0]),
               'area': (('y','x'), area)})
    
    return ds
        
def _load_tripolar(exp, collection):
    ds=xr.open_mfdataset(f'{exp["data_path"]}/{collection}/{exp["expid"]}.{collection}.monthly.??????.nc4')
    grid=xr.open_dataset(f'{exp["griddir"]}/MAPL_Tripolar.nc')
    dx=(grid['htn']+grid['hus'])*0.5
    dy=(grid['hte']+grid['huw'])*0.5
    ds=ds.assign_coords({'LON':(('lat','lon'),ds['LON'][0]), 
                         'LAT':(('lat','lon'),ds['LAT'][0])})
    # Rename coords for XESMF
    ds=ds.rename({'lon':'x', 'lat':'y'})
    ds=ds.rename({'LON': 'lon', 'LAT': 'lat'})
    ds.coords.update({'x': ds.lon[0,:], 'y': ds.lat[:,0]})

    ds.update({'dx': (('y','x'), dx), 
               'dy': (('y','x'), dy),
               'mask': (('y','x'), grid['mask']),
               'area': (('y','x'), grid['areat']),
               'basin': (('y','x'), grid['basin']),               
               'atl_mask': (('y','x'), grid['atl_mask'])})
    
    return ds

def _load_mom(exp, collection):
    ds=xr.open_mfdataset(f'{exp["data_path"]}/MOM_Output/{collection}.e??????01_00z.nc',decode_times=False)

    # Now need to decode times
    time=cftime.num2date(ds['time'],'days since 0001-01-01')
    # And convert to numpy.datetime64 (this is rediculous)
    ds.coords.update({'time': [np.datetime64(xx) for xx in time]})

    grid=xr.open_dataset(f'{exp["griddir"]}/MAPL_Tripolar.nc')
    ds.update({'mask': (('yh','xh'), grid['mask']),
               'area': (('yh','xh'), grid['areat']),
               'basin': (('yh','xh'), grid['basin']),               
               'atl_mask': (('yh','xh'), grid['atl_mask'])})

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
