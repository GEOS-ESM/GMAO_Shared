#!/usr/bin/env python3

'''
This module runs different plotting modules. Can be used as script.
'''
import sys
import geosdset

def main(exp_conf):
    
    # Load metadata for experiments to plot
    exps=geosdset.load_exps(exp_conf)

    # Load ocean 2d data
    ocn2d=geosdset.load_collection(exps, 'geosgcm_ocn2d')

    # Plot SST
    import sst
    sst.mkplots(exps, ocn2d)

    # Tropical Pacific SST plots
    import sst_tp
    sst_tp.mkplots(exps, ocn2d)    

    # Plot SSS
    import sss
    sss.mkplots(exps, ocn2d)

    # Plot surface wind stress
    import surf_stress
    surf_stress.mkplots(exps, ocn2d)

    geosdset.close(ocn2d)

    # Load surface data
    surf=geosdset.load_collection(exps, 'geosgcm_surf')

    # Plot precip
    import precip
    precip.mkplots(exps, surf)

    geosdset.close(surf)

    # Load ocean 3d data
    try:
        ocn3d=geosdset.load_collection(exps,'geosgcm_ocn3d')
    except OSError:
        ocn3d=geosdset.load_collection(exps,'prog_z',type='MOM')

    # Plot T profiles
    import temp
    temp.mkplots(exps, ocn3d)

    # Plot S profiles
    import salt
    salt.mkplots(exps, ocn3d)

    geosdset.close(ocn3d)

    # Add more plots....
    
if __name__ == "__main__":
    main(sys.argv[1])
