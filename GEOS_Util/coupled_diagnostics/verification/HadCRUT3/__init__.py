import os
import netCDF4 as nc
import scipy as sp
from g5lib import dset

__all__=['ctl']

class Ctl(dset.NCDset):

    def __init__(self):
        name='HadCRUT3'
        
        flist=['/discover/nobackup/projects/gmao/oceanval/verification/HadCRUT3/HadCRUT3.nc']
        f=nc.Dataset(flist[0])
        tt=f.variables['t']
        time=nc.num2date(tt[:],tt.units)
        undef=f.variables['temp'].missing_value
        f.close()

        super(Ctl,self).__init__(flist,time=time,undef=undef,name=name,
                                 lonname='longitude',
                                 latname='latitude',
                                 levname='unspecified',
                                 timename='t')
        
    def get_hadcrut3gl(self,annual=False):
        '''
        Gets HadCRUT3 globaly averaged T anomalies
        '''
        if not annual:
            # Return monthly data
            return sp.loadtxt('/discover/nobackup/projects/gmao/oceanval/verification/HadCRUT3/hadcrut3gl.txt', usecols=range(1,13))[0::2].flatten()
        else:
            # Return annual means
            return sp.loadtxt('/discover/nobackup/projects/gmao/oceanval/verification/HadCRUT3/hadcrut3gl.txt', usecols=range(1,13))[0::2].mean(1)
    

ctl=Ctl()
