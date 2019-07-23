'''
Contains Field class and field related routines.
'''

import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as pl
from matplotlib import dates, mlab
import netCDF4 as nc
import utils as utl
from domain import Domain
import plotters
import pca
import grid


class Field(object):
    """
    4d field of geographical data.
    Should contain data array, grid, time
    plotting subroutines, and other utilities
    
    data - 4d masked array
    time - array of datetime objects
    grid - Grid object
    name - Field name
    """
    def __init__(self,data,time,grid,name):
        self.data=data
        self.time=time
        self.grid=grid
        self.name=name
        
    def tonetcdf(self, fname, format='NETCDF4'):
        '''
        Dumps the field to netCDF file with metadata
        '''
        f=nc.Dataset(fname,mode='w',format=format)
        sh=self.dims
        
        f.createDimension('lon',sh[3])
        f.createDimension('lat',sh[2])
        f.createDimension('lev',sh[1])
        f.createDimension('time',None)
        
        f.createVariable('lon', 'f8', ('lon',)); f.variables['lon'][:]=self.grid['lon'][0]
        f.createVariable('lat', 'f8', ('lat',)); f.variables['lat'][:]=self.grid['lat'][:,0]
        f.createVariable('lev', 'f8', ('lev',)); f.variables['lev'][:]=self.grid['lev']
        f.createVariable('time', 'i8', ('time',))
        for i, tt in enumerate(self.time):
            f.variables['time'][i]=tt.toordinal()

        f.createVariable('LON','f8',('lat','lon',)); f.variables['LON'][:]=self.grid['lon']
        f.createVariable('LAT','f8',('lat','lon',)); f.variables['LAT'][:]=self.grid['lat']
        
        f.createVariable(self.name, 'f4', ('time','lev','lat','lon',)); f.variables[self.name][:]=self.data
        f.close()

    def subset(self, iind=slice(None), jind=slice(None), kind=slice(None), tind=slice(None)):
        '''
        Returns slice of the field.
        Indices should be valid for array indexing.
        If no indices are given, returns a copy.
        '''
        g=self.grid.subset(iind,jind,kind)
        
        i,j,k,t=utl.scalar2slice(iind,jind,kind,tind)
        time=self.time[t].copy()
        data=self.data[t][:,k][:,:,j][:,:,:,i].copy()

        out= Field(data,time,g,self.name)
        return out
    
    @property
    def dims(self):
        return self.time.shape+self.grid.dims

    @property
    def domain(self):

        dom=self.grid.domain
        dom['dates']=(self.time[0],self.time[-1])

        if dom['dates'][0]==dom['dates'][1]:
            dom['dates']=(dom['dates'][0],)

        return dom

    def index(self,**dom):
        '''
        Returns indexes given domain
        '''
        dd=self.domain
        dd.update(dom)

        return Domain(**dd)(self.grid,self.time)

    def __call__(self,**dom):
        '''
        Returns a subset of the field from a domain "dom". 
        '''
        
        return self.subset(**self.index(**dom))

    def gint(self, axis=-1):
        """
        Area integral of a field over a given axis.
        """
        # If axis size is one, return itself
        sh=self.data.shape
        if sh[axis]==1:
            return self
        
        r=6380e3
        X=sp.deg2rad(self.grid['lon']);Y=sp.deg2rad(self.grid['lat'])
        X*=sp.cos(Y)

        x=[dates.date2num(self.time),self.grid['lev'],\
           r*Y[sp.newaxis,sp.newaxis,:],r*X[sp.newaxis,sp.newaxis,:]]
        y=self.data.view(sp.ma.MaskedArray)
        # New dimensions
        newsh=list(sh); newsh[axis]=1
    
        # New grid and time
        ind=[slice(None)]*4; ind[axis]=slice(0,1)
        time=self.time[ind[0]]
        g=self.grid.subset(kind=ind[1],jind=ind[2],iind=ind[3])
        
        data=utl.gint(y,x[axis],axis);  data=data.reshape(newsh)
        
        return Field(data,time,g,self.name)

    def aint(self, area=None):
        '''
        Integral over area. area is a field type.
        If area is given as an argument, do sum over self.data*area.data,
        otherwise do consecutive integrals over lon and lats
        '''

        if area is None:
            X=self.gint(3).gint(2)
        else:
            Y=self.subset()
            Y.data*=area.data
            X=Y.subset(iind=0,jind=0)
            X.data[:,:,0,0]=Y.data.sum(3).sum(2)

        return X
        
    def sum(self, axis=0):
        ind=[slice(None)]*4; ind[axis]=slice(0,1)
        x=self.subset(tind=ind[0],kind=ind[1],jind=ind[2],iind=ind[3])
        newsh=x.dims
        x.data=self.data.sum(axis=axis).reshape(newsh)
        return x

    def mean(self, axis=0, ret_std=False, anom=False):
        ind=[slice(None)]*4; ind[axis]=slice(0,1)
        mean=self.subset(tind=ind[0],kind=ind[1],jind=ind[2],iind=ind[3])
        newsh=mean.dims
        mean.data=self.data.mean(axis=axis).reshape(newsh)

        if anom: self.data-=mean.data
            
        if ret_std:
            if anom:
                std=Field((self.data)**2,\
                          self.time,self.grid,self.name).mean(axis)
            else:
                std=Field((self.data-mean.data)**2,\
                          self.time,self.grid,self.name).mean(axis)
            std.data=sp.sqrt(std.data)
            return mean,std 
        else:
            return mean

    def ave(self, axis=-1, ret_std=False, anom=False):
        """
        Average field over a given axis.
        If ret_std=True - return standard deviation.
        If anom=True - remove mean from self.data
        """
        
        I=Field(sp.ma.masked_array(sp.ones(self.data.shape),mask=self.data.mask),\
                self.time,self.grid,self.name)
        I=I.gint(axis)
        
        mean=self.gint(axis)
        mean.data=mean.data/I.data

        if anom: self.data-=mean.data
            
        if ret_std:
            if anom:
                std=Field((self.data)**2,\
                          self.time,self.grid,self.name).ave(axis)
            else:
                std=Field((self.data-mean.data)**2,\
                          self.time,self.grid,self.name).ave(axis)
            std.data=sp.sqrt(std.data)
            return mean,std 
        else:
            return mean
        
    def aave(self, ret_std=False, anom=False, area=None):
        """
        Average field over area.
        If area is given as an argument, do sum over self.data*area.data,
        otherwise do consecutive integrals over lon and lats.
        If ret_std=True - return standard deviation.
        If anom=True - remove mean from self.data
        """
        
        mean=self.aint(area=area)

        if area is None:
            I=Field(sp.ma.masked_array(sp.ones(self.data.shape),mask=self.data.mask),\
                        self.time,self.grid,'')
            I=I.aint(area=area)
            mean.data=mean.data/I.data
        else:
            mean.data[:,:,0,0]/=area.data.sum(3).sum(2)
        
        if anom: self.data-=mean.data
        
        if ret_std:
            if anom:
                std=Field((self.data)**2,\
                          self.time,self.grid,self.name).aave(area=area)
            else:
                std=Field((self.data-mean.data)**2,\
                          self.time,self.grid,self.name).aave(area=area)
            std.data=sp.sqrt(std.data)
            return mean,std
        else:
            return mean

    def running_mean(self,N=1):
        '''
        Running mean
        '''
        Nt=self.data.shape[0]
        x=self.data.copy()
        for i in xrange(Nt):
            self.data[i]=x[max(i-N,0):min(i+N,Nt)].mean(0)

    def shiftgrid(self,lon0,lonsin=None,cyclic=360.0):
        '''
        Adapted basemap.shiftgrid.
        Shifts a field to origin at longitude lon0.
        Lonsin is a vector of longitudes. Works only on cyclic grid.
        lon0 - starting longitude for shifted grid
        '''

        if lonsin is None: lonsin=self.grid['lon'][0]
        i0=utl.find_nearest_index(lonsin,lon0)
        i0_shift = len(lonsin)-i0

        dataout  = sp.ma.zeros(self.data.shape,self.data.dtype)
        lonsout = sp.zeros(self.grid['lon'].shape,self.grid['lon'].dtype)
        latsout = sp.zeros(self.grid['lat'].shape,self.grid['lat'].dtype)

        lonsout[:,0:i0_shift] = self.grid['lon'][:,i0:]
        latsout[:,0:i0_shift] = self.grid['lat'][:,i0:]
        dataout[:,:,:,0:i0_shift] = self.data[:,:,:,i0:]

        lonsout[:,i0_shift:] = self.grid['lon'][:,:i0]+cyclic
        latsout[:,i0_shift:] = self.grid['lat'][:,:i0]
        dataout[:,:,:,i0_shift:] = self.data[:,:,:,:i0]
        
        self.data=dataout
        self.grid['lon']=lonsout
        self.grid['lat']=latsout

    def regrid(self, newgrid, newmask=None, interp='nearest'):
        '''
        Regrids a field to a grid defined in newgrid

        interp is interpolation method.
        '''
        X=newgrid['lon'];Y=newgrid['lat']

        newshape=self.time.shape+self.grid['lev'].shape+newgrid['lat'].shape
        nt,km,jm,im=newshape

        x=self.data.reshape(-1,self.grid['lat'].shape[0],self.grid['lon'].shape[1])
        
        xnew=sp.ma.masked_array(sp.zeros((nt*km,jm,im)))

        for ii,xx in enumerate(x):
            mm=~sp.ma.getmaskarray(xx)
            xnew[ii]=interpolate.griddata(zip(self.grid['lon'][mm],self.grid['lat'][mm]),xx[mm],(X,Y),method=interp)

        xnew=xnew.reshape(newshape)
        if newmask is not None: 
            xnew.mask=newmask
        xnew[sp.isnan(xnew)]=sp.ma.masked

        self.data=xnew            
        self.grid=grid.Grid(newgrid['lon'], newgrid['lat'], self.grid['lev'])
        
    def g2g(self,ogrid,tiledata):
        '''
        Interpolates field from one horizontal grid onto another using 
        exchange grid (tiles).
        
        ogrid - output grid

        tiledata - array of records of size Nt (Nt - number of tiles).
        tiledata['iin'], tiledata['jin'], tiledata['iout'], tiledata['jout'], 
        tiledata['frac']  - 
        arrays of size Nt of indexes on the input grid, indexes on the output grid, 
        corresponding to a tile N, and fraction of tile on the output grid.
        '''
        outdims=ogrid.dims[-2:]

        self.data=utl.g2g(self.data,tiledata,outdims)
        self.grid=grid.Grid(ogrid['lon'],ogrid['lat'],self.grid['lev'])


    def vinterp(self,newgrid, revert_axis=False, newmask=None):
        '''
        Performs vertical interpolation from self.grid['lev'] to newgrid['lev'].
        Levels should be monotonicaly increasing. If levels are decreasing, i.g. pressure levels,
        use revert_axis=True.
        '''

        if revert_axis:
            ind=slice(None,None,-1)
        else:
            ind=slice(None)

        lev=self.grid['lev'][ind]
        newlev=newgrid['lev'][ind]
        
        nt,km,jm,im=self.dims
        kmout=newlev.size

        data=self.data[:,ind,:,:]
        mm=data.mean()
        data[data.mask]=mm
        f=interpolate.interp1d(lev,data,axis=1)
        self.data=sp.ma.masked_values(f(newlev)[:,ind,:,:],mm)
        if newmask is not None: 
            self.data.mask=newmask

        self.grid=grid.Grid(self.grid['lon'],self.grid['lat'],newgrid['lev'])
        
    def d(self):
        '''
        Display the field.
        '''

        p=plotters.plotter(self.dims)
        p(self)
        return p



    def clim(self,freq,anom=False):
        '''
        Calculate annual or daily climatology based on frquency freq.
        For example, for monthly data, freq=12 will give annual cycle.
        If anom=True, subtracts a climatology from the data.
        '''
        
        cl=self.subset(tind=range(freq))
        
        for i,x in enumerate(cl.data):
            cl.data[i]=self.subset(tind=slice(i,None,freq)).data.mean(0)
            if anom: self.data[i::freq]-=cl.data[i]
        return cl

    def time_mean(self, freq):
        '''
        Makes a time average of field with frequency freq.
        For example, freq=12 will return annual field given monthly data,
        freq=24 will return daily field given hourly data.
        '''
        if self.time.size%freq != 0:
            raise Exception('Size of time series must be multiple of '+str(freq))
        var=self.subset(tind=slice(0,None,freq))
        var.data=self.data.reshape(-1,12,*self.grid.dims).mean(1)
            
        return var

    def detrend(self):
        '''
        Removes a time trend from data.
        '''

        # self.ave is awareof missing data, but sp.linalg.lstsq is not
        # need to remove mean first and set missing data to zeros
        mask=sp.ma.getmaskarray(self.data).copy()
        self.data[mask]=0.0; self.data.mask=mask; self.ave(0,anom=True)
        
        nt,km,jm,im=self.data.shape
        a=sp.mat(sp.ones((nt,2))); a[:,0]=sp.arange(nt).reshape(nt,1)
        b=sp.linalg.lstsq(a,self.data.reshape((nt,km*im*jm)).view(sp.matrix))
        self.data-=(a*b[0]).view(sp.ndarray).reshape(self.data.shape)

        return b

    def composites(self,ind):
        '''
        Makes composites of self based on index

        ind - a field object, representing index

        Returns two composite fields for self> std and self<-std 
        '''

        # Check if ind has same time scale as self

        if not(sp.all(self.time==ind.time)):
            raise Exception('Field time scale does not correspond to index time scale')

        std=ind.ave(0,ret_std=True)[1].data
        tplus=sp.where(ind.data>std)[0]
        tminus=sp.where(ind.data<-std)[0]

        pos=self.subset(tind=0); pos.data=self.subset(tind=tplus).data.mean(0)[sp.newaxis]
        neg=self.subset(tind=0); neg.data=self.subset(tind=tminus).data.mean(0)[sp.newaxis]
        
        return pos,neg

    def pca(self, keep=None, center=False, weight=True):
        '''
        Performss principal component analysis on data field, and stores
        a PCA object. Please, remove climatology, detrend data etc before
        calling this method. 

        If center=True, PCA object will center data using mean and standard deviation.
        If weight=True, multiply data by area weights.
        '''

        nt,km,jm,im=self.data.shape

        # multiply data by area factor, reshape, return matrix
        if weight:
            factor=sp.cos(sp.deg2rad(self.grid['lat']))
            factor[factor<0.]=0.
            factor=sp.sqrt(factor)
        else:
            factor=sp.ones(self.grid['lat'].shape)
        mask=sp.ma.getmaskarray(self.data).copy()
        self.data[mask]=0.0
        self.data*=factor[sp.newaxis,sp.newaxis]
        X=self.data.reshape((nt,km*jm*im)).view(sp.ndarray) 

        self._pc=pca.PCA(X, center=center, keep=keep)

        self.data/=factor[sp.newaxis,sp.newaxis]
        self.data[mask]=self.data.fill_value
        self.data.mask=mask

    def from_array(self,a):
        x=self.subset(tind=1,kind=0)
        sh=x.grid['lat'].shape
        x.data[0,0]=sp.ma.masked_array(a.reshape(sh),mask=x.data[0,0].mask)
        return x

    def get_pc(self,num=1):
        '''
        Returns the field which contains #num PC and EOF of a current field.
        '''

        if getattr(self,'_pc',None) is None: self.pca()

        e=self.from_array(sp.asarray(self._pc.eofs[:,num-1]))
        e.name=self.name+' EOF'+str(num)+' '+str(round(self._pc.fracs[num-1]*10000)/100)+'%'

        tmask=self.aave().data.mask
        p=self.subset(kind=0,jind=0,iind=0)
        p.data[:,0,0,0]=sp.asarray(self._pc.pcs)[:,num-1]
        p.data.mask=tmask
        p.name=self.name+' PC'+str(num)
        
        return p,e
            
def cmplx(re,imag):
    """
    Makes a complex field out of two real fields
    """
    F=re.subset()
    F.data=re.data+1j*imag.data

    return F

def absolute(F):
    """
    Returns absolute of complex field
    """
    out=F.subset()
    out.data=sp.absolute(F.data)
    out.data.set_fill_value(abs(F.data.fill_value))
    return out

