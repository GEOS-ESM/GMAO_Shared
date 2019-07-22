import os
import scipy as sp

__all__=['vsh','vsh_clim']

oceanval=os.environ.get('OCEANVAL',
                        '/discover/nobackup/projects/gmao/oceanval/verification')
path=oceanval+'/PIOMAS/RAW'
fmt='{path}/piomass.dat'

xx=sp.ma.masked_array(sp.loadtxt(fmt.format(path=path))[:,1:])
xx[xx==-1.0]=sp.ma.masked

vnh=xx.flatten(); vnh_clim=xx.mean(axis=0) 

