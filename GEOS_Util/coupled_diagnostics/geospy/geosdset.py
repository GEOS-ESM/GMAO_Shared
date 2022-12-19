'''
This modules provides utilities for reading different GEOS collections. 
'''

import os
import yaml, cftime
import numpy as np
import xarray as xr

R=6378e3 # Earth radius

def _load_geos(exp,vardata):
    colname=vardata['colname']
    templ=vardata.get('template','{data_path}/{colname}/{expid}.{colname}.monthly.??????.nc4')
    fnames=templ.format(data_path=exp['data_path'], colname=vardata['colname'], expid=exp['expid'])
    ds=xr.open_mfdataset(fnames)

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
        
def _load_tripolar(exp, vardata):
    colname=vardata['colname']
    templ=vardata.get('template','{data_path}/{colname}/{expid}.{colname}.monthly.??????.nc4')
    fnames=templ.format(data_path=exp['data_path'], colname=vardata['colname'], expid=exp['expid'])
    ds=xr.open_mfdataset(fnames)
    monthly='Xdim' in ds.coords.keys()
    grid=xr.open_dataset(f'{exp["griddir"]}/MAPL_Tripolar.nc')
    dx=(grid['htn']+grid['hus'])*0.5
    dy=(grid['hte']+grid['huw'])*0.5

    # Rename coords for XESMF
    if monthly:
        xname='Xdim'; yname='Ydim'; lonname='lons'; latname='lats'
    else:
        xname='lon'; yname='lat'; lonname='LON'; latname='LAT'

    ds=ds.assign_coords({lonname:((yname,xname),ds[lonname][0].values), 
                         latname:((yname,xname),ds[latname][0].values)})
    ds=ds.rename({xname:'x', yname:'y'})
    ds=ds.rename({lonname: 'lon', latname: 'lat'})
    ds.coords.update({'x': ds.lon[0,:], 'y': ds.lat[:,0]})

    ds.update({'dx': (('y','x'), dx.values), 
               'dy': (('y','x'), dy.values),
               'mask': (('y','x'), grid['mask'].values),
               'area': (('y','x'), grid['areat'].values),
               'basin': (('y','x'), grid['basin'].values),               
               'atl_mask': (('y','x'), grid['atl_mask'].values)})
    
    return ds

def _load_mom(exp, vardata):
    colname=vardata['colname']
    templ=vardata.get('template','{data_path}/MOM_Output/{colname}.e??????01_00z.nc')
    fnames=templ.format(data_path=exp['data_path'], colname=vardata['colname'])
    ds=xr.open_mfdataset(fnames, decode_times=False)

    # Now need to decode times
    time=cftime.num2date(ds['time'],'days since 0001-01-01')
    # And convert to numpy.datetime64 (this is rediculous)
    ds.coords.update({'time': [np.datetime64(xx) for xx in time]})

    grid=xr.open_dataset(f'{exp["griddir"]}/MAPL_Tripolar.nc')
    ds.update({'mask': (('yh','xh'), grid['mask'].values),
               'area': (('yh','xh'), grid['areat'].values),
               'basin': (('yh','xh'), grid['basin'].values),               
               'atl_mask': (('yh','xh'), grid['atl_mask'].values)})

    return ds 

def _get_loader(vardata):
    coltype=vardata['coltype']
    loaders={'GEOS': _load_geos,
             'GEOSTripolar': _load_tripolar,
             'MOM': _load_mom}
    return loaders[coltype]

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

def load_collection(exps, colname, coltype='GEOS'):
    '''
    Loads a collection 'colname' for all experiments in 'exps' list.
    Returns a list of xarray data sets.
    '''

    dsets=[]
    for exp in exps:
        dsets.append(_get_loader(coltype)(exp, colname))

    return dsets

def load_data(exps, plotname, defaults=None):
    '''
    Load data for plot "plotname" from experiments "exps".
    Use metadata (collection and variable name in collection) from exp 
    config if defaults is not None.
    '''
    dsets=[]
    for exp in exps:
        vardata=exp['plots'].get(plotname,defaults)

        if vardata is not None:
            dsets.append(_get_loader(vardata)(exp, vardata))
        else:
            raise Exception(f'''
            Metadata for {plotname} should be eigther passed as argument to "load_data"
            or set as options in exp config file.
            ''')

    return dsets

def close(dsets):
    for ds in dsets:
        ds.close()
