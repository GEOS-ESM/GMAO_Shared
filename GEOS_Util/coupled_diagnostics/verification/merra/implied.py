import os
from matplotlib import pyplot as pl
import scipy as sp
import netCDF4 as nc
import merra

# Read MERRA surface fluxes
lw=merra.Ctl('tavgM_2d_rad_Nx')('LWGNT').clim(12).mean(0)
sw=merra.Ctl('tavgM_2d_rad_Nx')('SWGNT').clim(12).mean(0)
lhf=merra.Ctl('tavgM_2d_flx_Nx')('EFLUX').clim(12).mean(0)
shf=merra.Ctl('tavgM_2d_flx_Nx')('HFLUX').clim(12).mean(0)
netheat=sw(); netheat.data+=(lw.data-lhf.data-shf.data)

# Read ocean mask
f=nc.Dataset('/discover/nobackup/dao_ops/scratch/d5_merra_jan79/stage/d5_merra_jan79.const_2d_asm_Nx.hdf')
frocean=sp.ma.masked_values(f.variables['FROCEAN'][0],0)
f.close()
netheat.data*=frocean
netheat.data-=netheat.aave().data

zonal=netheat.gint(3)
imptrans=zonal()
jm=imptrans.dims[2]-1
for j in xrange(jm-1,-1,-1):
    imptrans.data[:,:,j,:]=-zonal.subset(jind=slice(j,jm)).gint(2).data

imptrans.data/=1e15
GWlats=[47, 24, -19, -30];
GWerr=[0.1, 0.3, 0.6, 0.3];
GWoht=[0.6, 1.8, -0.8, -0.6];
trenb=sp.ma.masked_values(sp.genfromtxt('/home/yvikhlia/verification/implied/ANNUAL_TRANSPORTS_1985_1989.ascii',skiprows=1),-999.0)/100.

path=os.environ['NOBACKUP']+'/merra/plots'
try:
    os.makedirs(path)
except OSError:
    pass


pl.figure(); pl.clf()
imptrans.d();
pl.plot(trenb[:,0],trenb[:,6], linewidth=2,color='red')
pl.plot(GWlats,GWoht,'*',color='green')
pl.errorbar(GWlats,GWoht,yerr=GWerr,fmt='*',color='green')
pl.plot(trenb[:,0],trenb[:,6]+trenb[:,13],color='red')
pl.plot(trenb[:,0],trenb[:,6]-trenb[:,13],color='red')
ax=pl.gca()
ax.set_title('Implied Ocean Heat Transport')
ax.set_ylabel('PW')
ax.set_xlim(-90.,90)
ax.set_ylim(-3,3)
ax.legend(('MERRA','Trenberth-Caron','Ganachaud-Wunsch'),loc=4)
pl.grid()
pl.show()
pl.savefig(path+'/implied.png')
