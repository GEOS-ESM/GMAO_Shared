import scipy as sp
from utils import scalar2slice

class Grid(sp.ndarray):
    """
    This class contains metadata for Grads data set
    self['lon'] - array of longitudes
    self['lat'] - array of latitudes
    self['lev'] - array of levels

    Usage:
    g=Grid(lon,lat,lev)
    """

    def __new__(cls, lon=sp.zeros((1,1)), lat=sp.zeros((1,1)), lev=sp.zeros(1)):

        if (lon.ndim == 1) and (lat.ndim == 1):
            lon,lat=sp.meshgrid(lon,lat)
        
        im=lon.shape
        jm=lat.shape
        km=lev.shape
        
        dtype=sp.dtype([('lev',sp.float32,km),\
                        ('lat',sp.float32,jm),\
                        ('lon',sp.float32,im)])
        
        obj=sp.ndarray.__new__(cls,(),dtype=dtype)
        obj['lon']=lon; obj['lat']=lat; obj['lev']=lev
        return obj

    def subset(self,iind=slice(None),jind=slice(None),kind=slice(None)):
        """
        Make a subset according to indices
        """
        i,j,k=scalar2slice(iind,jind,kind)
        lon=self['lon'][j][:,i]; lat=self['lat'][j][:,i]; lev=self['lev'][k];
        return Grid(lon,lat,lev)
    
    @property
    def dims(self):
        km=self['lev'].shape
        jm=self['lat'].shape
        return km+jm

    @property
    def domain(self):
        dom={'lons': (self['lon'].min(),self['lon'].max()),\
             'lats': (self['lat'].min(),self['lat'].max()),\
             'levs': (self['lev'].min(),self['lev'].max())}

        if dom['lons'][0]==dom['lons'][1]:
            dom['lons']=(dom['lons'][0],)

        if dom['lats'][0]==dom['lats'][1]:
            dom['lats']=(dom['lats'][0],)

        if dom['levs'][0]==dom['levs'][1]:
            dom['levs']=(dom['levs'][0],)

        return dom
    
