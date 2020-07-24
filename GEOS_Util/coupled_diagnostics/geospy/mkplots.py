import geosdset

def main(expid):
    
    # Load metadata for obs and experiments to plot
    exps=geosdset.load_exps(expid)

    # Load ocean 2d data
    ocn2d=geosdset.load_collection(exps, 'geosgcm_ocn2d')

    # Plot SST
    import sst
    sst.plots(exps, ocn2d)
    
if __name__ == "__main__":
    import sys
    main(sys.argv[1])
