import os
import scipy as sp
from g5lib import dset
import datetime
import dateutil.rrule as rrule

class Ctl(dset.NCDset):
    def __init__(self,collection='tavgM_2d_flx_Nx'):
        '''
        collection: const_2d_asm_Nx  tavgM_2d_ocn_Nx  tavgM_3d_qdt_Cp
                    instM_2d_int_Nx  tavgM_2d_rad_Nx  tavgM_3d_rad_Cp
                    instM_3d_asm_Cp  tavgM_2d_slv_Nx  tavgM_3d_tdt_Cp
                    instM_3d_ana_Np  tavgM_2d_chm_Fx  tavgM_3d_trb_Cp
                    tavgM_2d_flx_Nx  tavgM_3d_cld_Cp  tavgM_3d_udt_Cp
                    tavgM_2d_int_Nx  tavgM_3d_mst_Cp
                    tavgM_2d_lnd_Nx  tavgM_3d_odt_Cp

        '''
        expid='MERRA'
        expdir='/discover/nobackup/projects/gmao/share/dao_ops/verification/MERRA_MEANS'
        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(1979,1,1),
                      until=datetime.date(2014,12,1))

        # Create meta-data
        flist=sp.array([expdir+'/MERRA.prod.assim.'
                        +collection+'.'
                        +str(date.year)
                        +str(date.month).zfill(2)
                        +'.hdf' 
                        for date in r[:]])

        time=sp.array(r[:],dtype='|O')

        super(Ctl,self).__init__(flist,time=time,name=expid, nrecs=sp.ones(flist.size),
                                 lonname='XDim',latname='YDim',timename='Time')


