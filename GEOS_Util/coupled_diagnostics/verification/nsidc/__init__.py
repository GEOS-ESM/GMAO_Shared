import os
import scipy as sp

__all__=['anh','ash','enh','esh','anh_clim','ash_clim','enh_clim','esh_clim']

Nyr=38
xx=sp.ma.zeros((2,12,Nyr,2))
fmt='{path}/{pole}_{mm}_extent_v2.1.csv'
oceanval=os.environ.get('OCEANVAL',
                        '/discover/nobackup/projects/gmao/oceanval/verification')
path=oceanval+'/NSIDC/RAW/MONTHLY_TXT'

pole='N'
flist=[fmt.format(path=path, pole=pole, mm=str(mm).zfill(2)) for mm in xrange(1,13)]
for ii,ff in enumerate(flist):
    xx[0,ii]=sp.loadtxt(ff, skiprows=1, delimiter=',', usecols=(4,5))[:Nyr]

pole='S'
flist=[fmt.format(path=path, pole=pole, mm=str(mm).zfill(2)) for mm in xrange(1,13)]
for ii,ff in enumerate(flist):
    xx[1,ii]=sp.loadtxt(ff, skiprows=1, delimiter=',', usecols=(4,5))[:Nyr]
    
xx[xx==-9999.]=sp.ma.masked

enh=xx[0,:,:,0].transpose().flatten(); enh_clim=xx[0,:,:,0].mean(axis=1)
anh=xx[0,:,:,1].transpose().flatten(); anh_clim=xx[0,:,:,1].mean(axis=1)
esh=xx[1,:,:,0].transpose().flatten(); esh_clim=xx[1,:,:,0].mean(axis=1)
ash=xx[1,:,:,1].transpose().flatten(); ash_clim=xx[1,:,:,1].mean(axis=1)
