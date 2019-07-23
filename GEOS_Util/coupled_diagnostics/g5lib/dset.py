import scipy as sp
from numpy.core.records import fromarrays
import netCDF4 as nc
from grid import Grid
from field import Field
from domain import Domain
from utils import scalar2slice


class Dset(object):
    """
    Base class for data set subclasses

    flist - list of files

    vlist - either list of variables, or list of tuples sutable for construction
    numpy.ndarray.dtype

    grid - dimensions of data set

    time - a time span of data set, array of datetime dates

    undef - undefined value

    name - name of dataset
    
    """
    def __init__(self,flist=[],vlist=[],grid=Grid(), time=sp.array(0,dtype='|O'),\
                 undef=1e15, name=''):

        print 'opening Dset'
        
        self.flist=flist
        self.vlist=vlist
        self.grid=grid
        self.name=name
        self.undef=undef
        self.time=time

    @property
    def domain(self):
        
        dom=self.grid.domain
        dom['dates']=(self.time[0],self.time[-1])

        if dom['dates'][0]==dom['dates'][1]:
            dom['dates']=(dom['dates'][0],)

        return dom

    def subset(self,iind=slice(None),\
               jind=slice(None),kind=slice(None),tind=slice(None), dtype=sp.float32):
        '''
        Creates empty field from dataset with dimentions (tind,kind,jind,iind)
        '''

        g=self.grid.subset(iind=iind,jind=jind,kind=kind)
        t,=scalar2slice(tind)
        time=self.time[t]

        sh=time.shape+g.dims
        data=sp.ma.masked_array(sp.zeros(sh),dtype=dtype,fill_value=self.undef)
        return Field(data,time,g, self.name)

    def fromfile(self,varname,iind=slice(None),\
                 jind=slice(None),kind=slice(None),tind=slice(None), dtype=sp.float32):
        '''
        Returns an empty field. To be overloaded in the derived class.
        '''

        var=self.subset(iind,jind,kind,tind)
        var.name=self.name+' '+varname

        return var
    
    def index(self,**dom):
        '''
        Returns indexes given domain
        '''
        dd=self.domain
        dd.update(dom)

        return Domain(**dd)(self.grid,self.time)

    def __call__(self,varname,fromfile=True, dtype=sp.float32, **dom):
        '''
        Returns a field "varname" from a domain "dom". If fromfile is true,
        actually reads the field, otherwise returns empty field. 
        '''
        
        ind=self.index(**dom)
        if fromfile:
            return self.fromfile(varname,dtype=dtype, **ind)
        else:
            return self.subset(dtype=dtype, **ind)

    def tave(self,varname,tind=slice(None)):
        '''
        Returns a field averaged over time interval defined by tind. Saves memory.
        Need to define fromfile functions in derived class first.

        (Needs to be remade as ave(varname,iind,jind,kind,tind,axis))
        
        DO NOT USE THIS
        '''

        t=sp.arange(self.time.size)[tind]
        
        var=self.subset(tind=t[0]);
        var.name+=' '+varname
        I=sp.ma.masked_array(sp.zeros(var.data.shape))
        for tt in t:
            x=self.fromfile(varname,tind=tt)
            var.data[0]=sp.ma.concatenate((var.data,x.data),0).sum(0)
            Ix=sp.ma.masked_array(sp.ones(var.data.shape),mask=x.data.mask)
            I.data[0]=sp.ma.concatenate((I.data,Ix),0).sum(0)

        var.data/=I

        return var

    def time_mean(self,varname,tind=slice(None)):
        '''
        Returns a field averaged over time interval defined by tind. 
        Saves memory compared to self.fromfile(varname).ave(0).
        '''

        t=sp.arange(self.time.size)[tind]
        var=self.subset(tind=t[0])
        var.data[:]=0.0
        for tt in t:
            var.data+=self.fromfile(varname,tind=tt).data
        var.data/=t.size
        return var
            
    def __del__(self):
        print 'deleting Dset'
    

class GADset(Dset):
    """
    This class contains metadata for multiple files data set in Grads format

    flist - list of file names

    vlist - either list of variables, or list of tuples sutable for construction
    numpy.ndarray.dtype

    grid - dimensions of data set

    time - a time span of data set, array of datetime dates
    
    """
    def __init__(self,flist=[],vlist=[],grid=Grid(), time=sp.array(0,dtype='|O'),\
                 undef=1e15, name='',mode='r'):

        print 'opening GADset'
        
        super(GADset,self).__init__(flist,vlist,grid,time,undef,name)

        # File map
        self.fmap=sp.zeros((time.size,2),dtype=sp.int32)
        self.files=sp.empty(len(flist),dtype='|O')

        srec=0
        rectot=0
        for i,fname in enumerate(flist):
            try:
                self.files[i]=sp.memmap(fname,dtype=vlist, mode=mode)
                nrec=self.files[i].shape[0]
            except IOError, err:
                print fname
                print err
                nrec=1
                
            erec=srec+nrec
            self.fmap[srec:erec,0]=i
            self.fmap[srec:erec,1]=sp.arange(srec,erec)-rectot
            srec=erec
            rectot+=nrec

    def fromfile(self,varname,iind=slice(None),\
                 jind=slice(None),kind=slice(None),tind=slice(None),dtype=sp.float32):
        '''
        Reads subset of the field. To be overloaded in the derived class.
        '''
        ii,jj,kk,tt=scalar2slice(iind,jind,kind,tind)

        var=self.subset(iind,jind,kind,tind)
        var.name=self.name+' '+varname

        for i,mm in enumerate(self.fmap[tt]):
            find=mm[0];  recind=scalar2slice(mm[1])
            print 'Reading '+varname+' from ',self.flist[find]
            if self.files[find] is not None:
                var.data[i]=sp.ma.masked_values(
                    self.files[find][varname][recind][:,kk][:,:,jj][:,:,:,ii],
                    self.undef)
            else:
                print 'File does not exist, values are set to missing'
                var.data[i]=self.undef; var.data[i]=sp.ma.masked

        return var

    def __del__(self):
        print 'deleting GADset'
        for f in self.files:
            del f
        del self.files
        super(GADset,self).__del__()
        
        
class NCDset(Dset):
    """
    This class contains metadata for multiple files netCDF data set

    Additional parameters:

    lonname, latname, levname, timename - names for longitude, latitude, level and time dimensions
    nrecs - container of same size as flist . Containes number of time records in each file. If missing,
    number of records will be read during initialization according to size of time dimension (can be slow
    on some data sets).
    """

    def __init__(self,flist=[],lonname='lon',latname='lat',levname='lev',\
                 timename='time',time=sp.array(0,dtype='|O'), undef=1e15, name='', nrecs=None):

        print 'opening NCDset'
        
        # File map
        self.fmap=sp.zeros((time.size,2),dtype=sp.int32)
        
        # Read grid info
        for fname in flist:
            try:
                ff=nc.Dataset(fname)
                break
            except (RuntimeError, IOError):
                ff=None

        if ff is not None:
            vlist={kk:{an:av for (an,av) in zip(vv.ncattrs(),[vv.getncattr(att) for att in vv.ncattrs()])} 
                   for (kk,vv) in ff.variables.items()}
            lat=ff.variables[latname][:]
            lon=ff.variables[lonname][:]
            try:
                lev=ff.variables[levname][:]
            except KeyError:
                lev=sp.zeros(1)
            ff.close()
        else:
            raise Exception('NCDset: none of files provided in file list exists.')

        # Fill file map
        srec=0
        rectot=0
        for i,fname in enumerate(flist):
            if nrecs is None:
                try:
                    f=nc.Dataset(fname)
                    nrec=f.variables[timename].shape[0]
                    f.close()
                except (RuntimeError, IOError):
                    nrec=1
            else:
                nrec=int(nrecs[i])
                
            erec=srec+nrec
            self.fmap[srec:erec,0]=i
            self.fmap[srec:erec,1]=sp.arange(srec,erec)-rectot
            srec=erec
            rectot+=nrec
            
        grid=Grid(lon,lat,lev)
            
        super(NCDset,self).__init__(flist,vlist,grid,time,undef,name)



    def fromfile(self,varname,iind=slice(None),jind=slice(None),kind=slice(None),\
                 tind=slice(None),maskandscale=True, dtype=sp.float32):
        '''
        Reads subset of the field. To be overloaded in the derived class.
        '''
        ii,jj,kk,tt=scalar2slice(iind,jind,kind,tind)
        
        var=self.subset(iind,jind,kind,tind,dtype=dtype)
        var.name=self.name+' '+varname
        
        for i,mm in enumerate(self.fmap[tt]):
            find=mm[0];  recind,=scalar2slice(mm[1])
            fname=self.flist[find]
            print 'Reading '+varname+' from ', fname
            try:
                f=nc.Dataset(fname)
                x=f.variables[varname]
                x.set_auto_maskandscale(maskandscale)
                x=x[recind]
                f.close()
            except (RuntimeError, IOError), err:
                print err
                x=sp.ma.masked_array(sp.ones(self.grid.dims)*self.undef,mask=True)
            
            if len(x.shape) < 4:
                x=x[:]; x=x[sp.newaxis,...]            

            var.data[i]=x[:,kk][:,:,jj][:,:,:,ii]

        return var

    def __del__(self):
        print 'deleting NCDset'
        super(NCDset,self).__del__()
