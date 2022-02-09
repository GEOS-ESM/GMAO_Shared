#!/usr/bin/env python3

'''
This module runs different plotting modules. Can be used as script.
'''
import sys, importlib
import geosdset

def main(exp_conf):
    
    # Load metadata for experiments to plot
    exps=geosdset.load_exps(exp_conf)

    plots=['sst', 'sst_tp', 'sss', 'surf_stress']

    for plot in plots:
        print(f'\nPlotting {plot}\n'.upper())
        importlib.import_module(plot).main(exps)

#    # Plot surface wind stress
#    import surf_stress
#    surf_stress.mkplots(exps, ocn2d)
#
#    # Plot surface wind stress in tropical Pacific
#    import zonal_stress_tp
#    zonal_stress_tp.mkplots(exps, ocn2d)
#
#    # Plot surface currents
#    import surf_current
#    surf_current.mkplots(exps, ocn2d)
#
#    geosdset.close(ocn2d)
#
#    # Load surface data
#    surf=geosdset.load_collection(exps, 'geosgcm_surf')
#
#    # Plot precip
#    import precip
#    precip.mkplots(exps, surf)
#
#    geosdset.close(surf)
#
#    # Load ocean 3d data
#    try:
#        ocn3d=geosdset.load_collection(exps,'geosgcm_ocn3d')
#    except OSError:
#        ocn3d=geosdset.load_collection(exps,'prog_z',type='MOM')
#
#    # Plot T profiles
#    import temp
#    temp.mkplots(exps, ocn3d)
#
#    # Plot S profiles
#    import salt
#    salt.mkplots(exps, ocn3d)
#
#    geosdset.close(ocn3d)
#
#    # Load MOM diagnostics
#    ocean_month=geosdset.load_collection(exps,'ocean_month',type='MOM')
#
#    # Plot AMOC stream function
#    import amoc
#    amoc.mkplots(exps,ocean_month)
#
#    geosdset.close(ocean_month)

    # Add more plots....
    
if __name__ == "__main__":
    main(sys.argv[1])
