import os
import pathlib, importlib
import xarray as xr

def load_ds(exp, collection, type='GEOS'):
    '''
    Returns an xarray dataset with GEOS collection.
    - exp: module with experiment metadata;
    - collection: collection name, string;
    - .
    '''
    
    if type=='GEOS':
        path=pathlib.Path(f'{exp.data_path}/{collection}')
        flist=list(path.glob(f'{exp.expid}.{collection}.monthly.?????[0-9].nc4')) 
        flist.sort()
        ds=xr.open_mfdataset(flist,combine='by_coords')
    elif type=='GEOSTripolar':
        path=pathlib.Path(f'{exp.data_path}/{collection}')
        flist=list(path.glob(f'{exp.expid}.{collection}.monthly.?????[0-9].nc4')) 
        flist.sort()
        ds=xr.open_mfdataset(flist,combine='by_coords')
        ds=ds.assign_coords({'LON': ds['LON'],
                             'LAT': ds['LAT']})
    elif type=='MOM':
        path=pathlib.Path(f'{exp.data_path}/MOM_Output')
        flist=list(path.glob(f'{collection}.e*_00z.nc'))
        flist.sort()
        ds=xr.open_mfdataset(flist,combine='by_coords')
        ds_static=xr.open_dataset(exp.ocean_static)
        ds=ds.assign_coords({'geolon': ds_static['geolon'],
                             'geolat': ds_static['geolat'],
                             'geolon_c': ds_static['geolon_c'],
                             'geolat_c': ds_static['geolat_c'],
                             'geolon_u': ds_static['geolon_u'],
                             'geolat_u': ds_static['geolat_u'],
                             'geolon_v': ds_static['geolon_v'],
                             'geolat_v': ds_static['geolat_v']})

    return ds    

def load_exps(expid):
    '''
    Returns a list of module withmeta data for experiment we want to plot (expid) and 
    all experiments we want to compare to.
    '''
    
    exps=[]
    exps.append(importlib.import_module(expid))
    for exp in exps[0].cmpexp:
        exps.append(importlib.import_module(exp))

    # Make directory for plots if it does not exists
    try:
        os.makedirs(exps[0].plot_path)
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
        dsets.append(load_ds(exp,colname,type))

    return dsets
