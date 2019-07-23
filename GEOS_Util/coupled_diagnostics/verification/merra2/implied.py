import os
from matplotlib import pyplot as pl
import scipy as sp
import netCDF4 as nc
import merra2

# Read MERRA2 surface fluxes
lw=merra2.CtlClim('tavgM_2d_rad_Nx')('lwgnt').mean(0)
sw=merra2.CtlClim('tavgM_2d_rad_Nx')('swgnt').mean(0)
lhf=merra2.CtlClim('tavgM_2d_flx_Nx')('eflux').mean(0)
shf=merra2.CtlClim('tavgM_2d_flx_Nx')('hflux').mean(0)
netheat=sw(); netheat.data+=(lw.data-lhf.data-shf.data)

# Read ocean mask
f=nc.Dataset('/gpfsm/dnb02/projects/p53/merra2/scratch/MERRA2_100/stage/MERRA2_100.const_2d_asm_Nx.00000000.nc4')
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

path=os.environ['NOBACKUP']+'/verification/merra2/plots'
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
ax.legend(('MERRA-2','Trenberth-Caron','Ganachaud-Wunsch'),loc=4)
pl.grid()
pl.show()
pl.savefig(path+'/implied.png')
