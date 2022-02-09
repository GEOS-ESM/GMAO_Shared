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

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    mkplots(exps,dsets)
    geosdset.close(dsets)
