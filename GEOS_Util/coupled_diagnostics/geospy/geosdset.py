import importlib
import pathlib
import xarray as xr

def load_ds(expid, collection, type='GEOS'):
    '''
    Returns an xarray dataset with GEOS collection.
    - expid: experiment id;
    - collectiond: collection name, string.
    '''
    
    exp=importlib.import_module(expid)
    if type == 'GEOS':
        path=pathlib.Path(f'{exp.data_path}/{collection}')
        flist=list(path.glob(f'{expid}.{collection}.monthly.?????[0-9].nc4')) 
    elif type=='MOM':
        path=pathlib.Path(f'{exp.data_path}/MOM_Output')
        flist=list(path.glob(f'{collection}.e*_00z.nc'))

    flist.sort()
    ds=xr.open_mfdataset(flist,combine='by_coords')

    return ds,exp
    
