import os
import scipy as sp
import netCDF4 as nc

dir=os.environ['NOBACKUP']+'/workdir/move_runoff'

# Read ocean land mask
print 'Reading ocean land mask'
gfile=dir+'/grid_spec.nc'
ff=nc.Dataset(gfile)
lon=ff.variables['x_T'][:]
lat=ff.variables['y_T'][:]
wet=ff.variables['wet'][:]
sh=wet.shape
tmp=sp.zeros((sh[0]+2,sh[1]+2))
tmp[1:-1,1:-1]=wet
nwet=tmp[1:-1,1:-1]+tmp[:-2,:-2]+tmp[1:-1,:-2]\
    +tmp[2:,:-2]+tmp[2:,1:-1]+tmp[2:,2:]\
    +tmp[1:-1,2:]+tmp[:-2,2:]+tmp[:-2,1:-1]

# Make slmask
slmask=wet.copy()
slmask[nwet<6.0]=0.0

# Read routing table
print 'Reading routing data'
runoff_table=dir+'/runoff.bin'

Nt=sp.fromfile(runoff_table,dtype='i4',count=3)[1]

vlist=[('','i4'),('Nt','i4'),('','i4'),
       ('','i4'),('sind','i4',Nt),('','i4'),
       ('','i4'),('dind','i4',Nt),('','i4'),
       ('','i4'),('weight','f4',Nt),('','i4')]

route=sp.fromfile(runoff_table,dtype=vlist); route['sind']-=1; route['dind']-=1 # indises start with 0 

# Read ocean i,j corresponding to tmp from tile file
print 'Reading tile data '
tfile=dir+'/tile.dat'
tmp=sp.loadtxt(tfile, skiprows=8,dtype='f4'); tdata=sp.hstack((tmp,sp.array([sp.arange(tmp.shape[0])]).transpose()))
ijs=tdata[:,9:7:-1].astype(int); ijs-=1
# Subset ocean only tiles 
otdata=tdata[sp.where(tdata[:,0]==0)[0]]
oijs=otdata[:,9:7:-1].astype(int); oijs-=1
# Remove tiles which are over ocean land (slmask[i,j]==0)
tslmask=sp.zeros(otdata.shape[0])
for ii,oij in enumerate(oijs):
    gind=(oij[0],oij[1])
    tslmask[ii]=slmask[gind]
otdata_noland=otdata[sp.where(tslmask > 0)[0]]

# Construct a pretty array with tile data
rec=[('area','f4'),('lon','f4'),('lat','f4'),
     ('ii','i4'),('jj','i4'),('tnum','i4'),
     ('sinlat','f4'), ('coslat_coslon','f4'), ('coslat_sinlon','f4')]
tiles=sp.zeros(otdata_noland.shape[0],dtype=rec)
tiles['area']=otdata_noland[:,1]; tiles['lon']=otdata_noland[:,2]; tiles['lat']=otdata_noland[:,3]; 
tiles['ii']=otdata_noland[:,8]; tiles['jj']=otdata_noland[:,9]
tiles['tnum']=otdata_noland[:,12]
#tmp=tdata[sp.where(tdata[:,0]==0)[0]]; otiles=sp.hstack((tmp,sp.zeros((tmp.shape[0],3))))
##sp.save(dir+'/sind.npy',sind)
##sind=sp.load(dir+'/sind.npy')
#
#
#for tind in route['dind'][0]:
#    gind=(tdata[tind][9]-1,tdata[tind][8]-1)
#    if not slmask[gind]:
        # Find another dind

# Move discharge.
#print 'Moving discharge'
#def move_discharge(ji,nwet):
#  sj,ej=max(ji[0]-2,0),min(ji[0]+3,sh[0])
#  si,ei=max(ji[1]-2,0),min(ji[1]+3,sh[1])
#  tmp=nwet[sj:ej,si:ei]
#  if tmp.max() >= 6.0:
#    j,i=sp.unravel_index(tmp.argmax(),tmp.shape)
#    ji[0],ji[1]=sj+j,si+i
#  return ji
#  
#
#dind=sind.copy()
#for i,ji in enumerate(dind):
#  if nwet[ji[0],ji[1]] < 6.0 or wet[ji[0],ji[1]]==0.0:
#    # move discharge
#    dind[i]=move_discharge(ji,nwet)
#
## Move discharge out of Uruguay Bay
#ind=sp.argmax(sp.all(sind==(350,891),axis=1)); dind[ind]=sp.array((345,896))
#ind=sp.argmax(sp.all(sind==(349,891),axis=1)); dind[ind]=sp.array((344,896))
#ind=sp.argmax(sp.all(sind==(350,890),axis=1)); dind[ind]=sp.array((345,895))
#ind=sp.argmax(sp.all(sind==(349,889),axis=1)); dind[ind]=sp.array((344,894))
#ind=sp.argmax(sp.all(sind==(349,887),axis=1)); dind[ind]=sp.array((344,892))
#ind=sp.argmax(sp.all(sind==(351,888),axis=1)); dind[ind]=sp.array((346,893))
#ind=sp.argmax(sp.all(sind==(352,887),axis=1)); dind[ind]=sp.array((345,897))
#ind=sp.argmax(sp.all(sind==(347,891),axis=1)); dind[ind]=sp.array((342,896))
#ind=sp.argmax(sp.all(sind==(349,893),axis=1)); dind[ind]=sp.array((344,898))
#
## Move discharge out of Amazon Bay
#ind=sp.argmax(sp.all(sind==(502,919),axis=1)); dind[ind]=sp.array((505,924))
#ind=sp.argmax(sp.all(sind==(500,918),axis=1)); dind[ind]=sp.array((500,924))
#ind=sp.argmax(sp.all(sind==(497,922),axis=1)); dind[ind]=sp.array((502,927))
#ind=sp.argmax(sp.all(sind==(498,921),axis=1)); dind[ind]=sp.array((498,926))
#ind=sp.argmax(sp.all(sind==(501,921),axis=1)); dind[ind]=sp.array((501,926))
#ind=sp.argmax(sp.all(sind==(499,923),axis=1)); dind[ind]=sp.array((499,928))
#
#disch_old=wet.copy(); disch_old[sind[:,0,],sind[:,1]]=2
#disch_new=wet.copy(); disch_new[dind[:,0,],dind[:,1]]=2
#
#ind=sp.unique(sp.where(sind != dind)[0])
#tmp=sp.concatenate((sp.roll(sind[ind]+1,1,1),sp.roll(dind[ind]+1,1,1)),1)
#sp. savetxt('/home/yvikhlia/nobackup/workdir/move_runoff/ocean_routing_1440x1080.txt',tmp,fmt='%4.1i')
