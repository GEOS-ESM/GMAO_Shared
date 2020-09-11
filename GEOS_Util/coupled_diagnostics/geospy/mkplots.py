#!/usr/bin/env python3

'''
This module runs different plotting modules. Can be used as script.
'''

import geosdset

def main(expid):
    
    # Load metadata for experiments to plot
    exps=geosdset.load_exps(expid)

    # Load ocean 2d data
    ocn2d=geosdset.load_collection(exps, 'geosgcm_ocn2d')

    # Plot SST
    import sst
    sst.plots(exps, ocn2d)

    # Load ocean 3d data
    ocn3d=geosdset.load_collection(exps, 'prog_z', type='MOM')

    # Plot T profiles
    import temp
    temp.plots(exps, ocn3d)

    # Add more plots....
    
if __name__ == "__main__":
    import sys
    main(sys.argv[1])
