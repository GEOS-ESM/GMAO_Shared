import pathlib, importlib
import xarray as xr

def load_ds(exp, collection, type='GEOS'):
    '''
    Returns an xarray dataset with GEOS collection.
    - exp: module with experiment metadata;
    - collection: collection name, string.
    '''
    
    if type == 'GEOS':
        path=pathlib.Path(f'{exp.data_path}/{collection}')
        flist=list(path.glob(f'{exp.expid}.{collection}.monthly.?????[0-9].nc4')) 
    elif type=='MOM':
        path=pathlib.Path(f'{exp.data_path}/MOM_Output')
        flist=list(path.glob(f'{collection}.e*_00z.nc'))

    flist.sort()
    ds=xr.open_mfdataset(flist,combine='by_coords')

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

    return exps

def load_collection(exps, colname):
    '''
    Loads a collection 'colname' for all experiments in 'exps' list.
    Returns a list of xarray data sets.
    '''

    dsets=[]
    for exp in exps:
        dsets.append(load_ds(exp,colname))

    return dsets
