#!/usr/bin/env bash
# Modules for ioda-converters
# ---------------------------
#jedisrc=/discover/nobackup/projects/gmao/advda/swell/JediBundles/soca_ioda_SLES15
jedisrc=/discover/nobackup/projects/gmao/advda/swell/JediBundles/fv3_soca_SLES15_02062025
jedibuild=$jedisrc/build-intel-release

module purge
source $jedibuild/modules
#module load jedi_bundle
# Ioda Python
export LD_LIBRARY_PATH=$jedibuild/lib:$LD_LIBRARY_PATH

#what worked for us
export PYTHONPATH=$jedibuild/lib/python3.11/:$PYTHONPATH
export PYTHONPATH=$jedibuild/lib/python3.11/pyioda:$PYTHONPATH
export PYTHONPATH=$jedibuild/lib/python3.11/pyiodaconv:$PYTHONPATH
export PYTHONPATH=$jedibuild/lib/python3.11/pyiodautils:$PYTHONPATH

#instructions on the ioda-conv webpage (not working)
#export PYTHONPATH=$jedibuild/lib/pyiodaconv:$PYTHONPATH
#export PYTHONPATH=$jedisrc/iodaconv/src:$PYTHONPATH

export PATH=$jedibuild/bin/gmao_obs2ioda.py:$PATH
# List modules
ml

# Convert
#python $jedibuild/bin/proc_gsi_ncdiag.py -o gsi/ ioda
