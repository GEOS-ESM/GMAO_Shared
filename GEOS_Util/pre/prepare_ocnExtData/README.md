This directory contains scripts, fortran code to write out [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) formatted file(s) that could be
used as external dataset(s) in various applications of [GEOS-ESM](https://github.com/GEOS-ESM), for example:
- [GEOSgcm](https://github.com/GEOS-ESM/GEOSgcm)
- [GEOSadas](https://github.com/GEOS-ESM/GEOSadas)

It also has a few notebooks, see `notebooks/` that could be used to plot the data.

1. `write_kPAR.py`: A script to write out a climatological k-PAR (Photosynthetically Available Radiation) file.

**Note**:
 - For all python scripts, python3 or later version is needed.
