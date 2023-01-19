#!/usr/bin/env python3

'''
Docstring
'''
import sys
import numpy as np
import matplotlib.pyplot as pl
import geosdset

def mkplots(exps,dsets):
    pass

def main(exps):
    dsets=geosdset.load_data(exps, plotname, defaults)
    mkplots(exps,dsets)
    geosdset.close(dsets)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    main(exps)
