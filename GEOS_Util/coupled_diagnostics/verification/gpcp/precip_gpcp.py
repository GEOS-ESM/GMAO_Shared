import os
import matplotlib.pyplot as pl
from mpl_toolkits.basemap.cm import s3pcpn_l
import my_plots as mpl

# GPCP precipitation
path=os.environ['NOBACKUP']+'/verification/gpcp'
execfile(path+'/ctl.py')

clim=ctl.fromfile('precip').clim(12)
clim.shiftgrid(30.)
clim.grid['lon']=sp.where(clim.grid['lon']<28.,clim.grid['lon']+360.,clim.grid['lon'])

ind=[0,1,11]; djf=clim.subset(tind=ind).ave(0); djf.name=ctl.name+' '+djf.name+' DJF'
ind=[5,6,7]; jja=clim.subset(tind=ind).ave(0); jja.name=ctl.name+' '+jja.name+' JJA'
am=clim.ave(0); am.name=ctl.name+' '+am.name+' Annual Mean'


path=os.environ['NOBACKUP']+'/verification/gpcp/pics'

def plot_field(field, fig):
    clevs=(.2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.3, 2.6, 3, \
           3.3, 3.6, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 11, 13, 15)
    pl.figure(fig); pl.clf()
    field.copts={'levels': clevs,\
                 'cmap': s3pcpn_l}
    field.plot_map()
    mean,std=field.aave(ret_std=True)
    mpl.draw_stat((mean.data.squeeze(),std.data.squeeze()))
    pl.show()
    
# DJF
plot_field(djf,1)
pl.savefig(path+'/precip_gpcp_djf.png')

# JJA
plot_field(jja,2)
pl.savefig(path+'/precip_gpcp_jja.png')

# Annual Mean
plot_field(am,3)
pl.savefig(path+'/precip_gpcp_am.png')


