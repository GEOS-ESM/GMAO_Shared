import os
import scipy as sp
from scipy import interpolate
import glob
import netCDF4 as nc

# Input data

imi,jmi=1440,720
dlon,dlat=360./imi,180./jmi
lati,loni=sp.mgrid[-90+dlat/2:90:dlat,-180+dlon/2:180:dlon]
loni[loni>79.875]-=360.

itype=[]
itype.append(('head','f4',16))
itype.append(('h1','i4',1))
itype.append(('var','f4',(jmi,imi)))
itype.append(('h2','i4',1))

var='SST'
idir='/discover/nobackup/projects/gmao/share/dao_ops/fvInput/g5gcm/bcs/SST/'+str(imi)+'x'+str(jmi)+'/'
iflist=glob.glob1(idir,'dataoceanfile_MERRA2_'+var+'*data'); iflist.sort()

mask=sp.fromfile('/gpfsm/dnb02/projects/p43/bzhao/oceanmask_merra2_1440x720.bin', dtype=('float32',(jmi,imi)))[0]>0

# Output data

with nc.Dataset('/discover/nobackup/yvikhlia/coupled/Forcings/a90x540_o360x200/INPUT/grid_spec.nc') as ogfile:
    lono=ogfile.variables['geolon_t'][:]; lato=ogfile.variables['geolat_t'][:]
    jmo,imo=lono.shape
    wet=ogfile.variables['wet'][:]

otype=[]
otype.append(('head','f4',16))
otype.append(('h1','i4',1))
otype.append(('var','f4',(jmo,imo)))
otype.append(('h2','i4',1))

odir=os.environ['NOBACKUP']+'/workdir/SST'
try:
    os.makedirs(odir)
except OSError:
    pass

for ff in iflist:
    of=ff.replace(str(imi)+'x'+str(jmi),str(imo)+'x'+str(jmo))
    print ff+' --> '+of
    aa=sp.memmap(idir+'/'+ff,dtype=itype,mode='r')
    zout=sp.zeros(aa.size,dtype=otype)
    zout['h1']=zout['h2']=imo*jmo*4
    zout['head']=aa['head']
    zout['head'][:,13]=imo; zout['head'][:,14]=jmo
    for ii,zin in enumerate(aa['var']):
        zout['var'][ii]=sp.where(wet==0,1e15,
                                 interpolate.griddata(zip(loni[mask],lati[mask]),zin[mask],(lono,lato),method='nearest',fill_value=1e15))

    zout.tofile(odir+'/'+of)
