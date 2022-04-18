import os
from g5lib import dset, grid 

class Ctl(Dset):
    def __init__(self):
        name='GPCP'
        undef=-1e10

        flist=[os.environ['NOBACKUP']+'/verification/coads/coads2x25_climate.dat']

        lon=sp.arange(0.0,360.,2.5); lat=sp.arange(-90.,91,2.); lev=sp.zeros(1)

        im,jm,km=144,91,1
        vlist=[('ac','>f4',(km,jm,im))]
        
        super(Ctl,self).__init__()
