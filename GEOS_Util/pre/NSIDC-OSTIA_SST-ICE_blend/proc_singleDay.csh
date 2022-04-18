#!/bin/csh

set year      = $1
set month     = $2
set day       = $3
set whichSST=$4
set toPlot    = $5

set build_dir = /discover/nobackup/sakella/geosMom6/develop_24Jul2020/GEOSgcm/build-Debug/bin/
# ---------------------------------------------------------

# generate input 
if ( -e input.txt) /bin/rm -f input.txt
gen_input.csh ${year} ${month} ${day} ${whichSST}

# generate SST and SIC files
${build_dir}/lake_sst_sic_EIGTHdeg.x input.txt

# clean up
/bin/rm -f  ${PWD}/${whichSST}/*.nc
if (${whichSST} == 'OSTIA') then
  /bin/mv    Ostia_*.bin ${PWD}/${whichSST}/
  /bin/mv    mask*.bin   ${PWD}/${whichSST}/
else
 echo "Unknown option for whichSST:" ${whichSST} ". Exiting"
 exit 1
endif

# plot if needed
if (${toPlot} == 'yes') then
  plot_bcs.py -year ${year} -month ${month} -day ${day} -fndSST ${whichSST}
  /bin/mv *.png ${PWD}/${whichSST}/
endif
# ---------------------------------------------------------

# example usage:
# -------------
# proc_singleDay.csh 2020 3 1 OSTIA 'yes'  <-- after 03/01 is okay, file with prefix exist: 'oisst-avhrr-v02r01.'
# proc_singleDay.csh 2020 7 1 OSTIA 'yes'

# gen_input.csh      2020 7 1 OSTIA
