import os
import scipy as sp
from g5lib import dset
import datetime
import dateutil.rrule as rrule

class CtlClim(dset.NCDset):
    def __init__(self,collection='tavgM_2d_flx_Nx'):
        '''
        collection: instM_2d_asm_Nx  tavgM_2d_adg_Nx  tavgM_2d_lfo_Nx  
                    tavgM_3d_odt_Np  instM_2d_gas_Nx  tavgM_2d_aer_Nx  
                    tavgM_2d_lnd_Nx  tavgM_3d_qdt_Np  instM_2d_int_Nx  
                    tavgM_2d_chm_Nx  tavgM_2d_ocn_Nx  tavgM_3d_rad_Np
                    instM_2d_lfo_Nx  tavgM_2d_csp_Nx  tavgM_2d_rad_Nx  
                    tavgM_3d_tdt_Np  instM_3d_asm_Np  tavgM_2d_flx_Nx  
                    tavgM_2d_slv_Nx  tavgM_3d_trb_Np  tavgM_2d_glc_Nx  
                    tavgM_3d_cld_Np  tavgM_3d_udt_Np  statM_2d_slv_Nx  
                    tavgM_2d_int_Nx  tavgM_3d_mst_Np
        '''
        expid='MERRA-2'
        expdir='/discover/nobackup/mbosilov/MERRA2/Climate/'+collection
        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(1980,1,1),
                      count=12)

        # Create meta-data
        flist=sp.array([expdir+'/MERRA2.clim.'+collection+'.nc'])

        time=sp.array(r[:],dtype='|O')

        super(CtlClim,self).__init__(flist,time=time,name=expid,
                                     lonname='longitude',latname='latitude')


class Ctl(dset.NCDset):
    def __init__(self,collection='tavg1_2d_flx_Nx'):
        '''
        collection: same as in CtlClim

        '''
        expid='MERRA-2'
        expdir='/discover/nobackup/projects/gmao/share/dao_ops/verification/MERRA2_MEANS/'+collection
        
        r=rrule.rrule(rrule.MONTHLY,
                      dtstart=datetime.date(1980,1,1),
                      until=datetime.date(2016,12,1))
        flist=sp.array([expdir+'/'+expid+'.'+collection+'.monthly.'
                        +str(date.year)+str(date.month).zfill(2)+'.nc4'
                        for date in r[:]
                        ])
        time=sp.array(r[:],dtype='|O')

        super(Ctl,self).__init__(flist,time=time,name=expid,nrecs=sp.ones(flist.size))
