from importlib import import_module
import datetime
import dateutil.rrule as rrule
import scipy as sp 
import netCDF4 as nc
import dset

class Ctl(dset.NCDset):
    def __init__(self,exp,collection,freq='MONTHLY'):

        fmt=getattr(exp,'fmt',{}).get(collection,{})

        if freq=='MONTHLY':
            r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(exp.start_year,1,15),until=datetime.date(exp.end_year,12,15))
            fmt=fmt.get('MONTHLY',
                        '{data_path}/{collection}/{expid}.{collection}.monthly.{date:%Y%m}.nc4')
        elif freq=='DAILY':
            r=rrule.rrule(rrule.DAILY,dtstart=datetime.date(exp.start_year,1,1),until=datetime.date(exp.end_year,12,31))
            fmt=fmt.get('DAILY',
                        '{data_path}/{collection}/{expid}.{collection}.monthly.{date:%Y%m%d}_1500z.nc4')
        else:
            raise Exception('Unsupported frequency '+freq)

        flist=sp.array([fmt.format(data_path=exp.data_path,collection=collection,expid=exp.expid,date=dd) for dd in r[:]])
        time=sp.array(r[:],dtype='|O')
        super(Ctl,self).__init__(flist,time=time,name=exp.expid, nrecs=sp.ones(flist.size))

class CtlTripolar(Ctl):
    def __init__(self,exp,collection,freq='MONTHLY'):

        super(CtlTripolar,self).__init__(exp,collection,freq)

        for fname in self.flist:
            try:
                f=nc.Dataset(fname)
                vlist=f.variables
                lat=vlist['LAT'][:]
                lon=vlist['LON'][:]
                f.close()
                break
            except RuntimeError, err:
                print fname
                print err
        sh=lon.shape; ind=[0]*len(sh); ind[-1]=ind[-2]=slice(None)
        self.grid['lon']=lon[ind]
        self.grid['lat']=lat[ind]

class CtlMOM(dset.NCDset):
    def __init__(self,exp):

        collection='MOM_Output'
        lonname='xt_ocean'
        latname='yt_ocean'
        levname='st_ocean'
        
        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(exp.start_year,1,15),until=datetime.date(exp.end_year,12,15))
        fmt=getattr(exp,'fmt',{}).get(collection,
                                    '{data_path}/{collection}/ocean_month.e{date:%Y%m}01_00z.nc')
        flist=sp.array([fmt.format(data_path=exp.data_path,collection=collection,date=dd) for dd in r[:]])

        time=sp.array(r[:],dtype='|O')
        super(CtlMOM,self).__init__(flist,lonname=lonname,latname=latname,levname=levname,\
                                 time=time,name=exp.expid, nrecs=sp.ones(flist.size))
        self.basin_mask=exp.basin_mask
        self.grid_spec=exp.grid_spec
        f=nc.Dataset(self.grid_spec)

        vlist=f.variables
        if vlist.has_key('x_T'):
            lon=vlist['x_T'][:]
            lat=vlist['y_T'][:]
        elif vlist.has_key('geolon_t'):
            lon=vlist['geolon_t'][:]
            lat=vlist['geolat_t'][:]
        else:
            raise Exception("Can't find latitudes and longitudes")

        self.grid['lon']=lon
        self.grid['lat']=lat
        f.close()

def read_exp(expid):
    return import_module(expid)
