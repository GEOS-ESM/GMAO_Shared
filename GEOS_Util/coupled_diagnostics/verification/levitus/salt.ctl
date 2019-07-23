dset /gpfsm/dnb42/projects/p16/ssd/ocean/kovach/odas-2/val/LEVITUS_GRD/salt/salt.%m2
undef -100

options template sequential 

xdef 360 linear 0.5 1.0
ydef 180 linear -89.5 1.0
zdef 33 levels 0.    10.    20.    30.    50.    75.   100.   125. 150.   200.   250.   300.   400.   500.   600.   700. 800.   900.  1000.  1100.  1200.  1300.  1400.  1500. 1750.  2000.  2500.  3000.  3500.  4000.  4500.  5000. 5500.
tdef 12 linear jan2009 1mo

vars 1
salt 33 0 salinity
endvars
