#!/usr/bin/env python3

'''
This module runs different plotting modules. Can be used as script.
'''
import sys, importlib
import geosdset

def main(exp_conf):
    
    # Load metadata for experiments to plot
    exps=geosdset.load_exps(exp_conf)

    plot_modules=['sst', 'sst_tp', 'sss', 'surf_stress', 
                  'zonal_stress_tp', 'precip', 'temp', 'salt']

    for module in plot_modules:
        print(f'\nPlotting {module}\n'.upper())
        importlib.import_module(module).main(exps)

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
