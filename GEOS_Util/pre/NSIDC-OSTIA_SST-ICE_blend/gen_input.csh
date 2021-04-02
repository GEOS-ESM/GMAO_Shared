#!/bin/csh

set year    = $1
set month_  = $2
set today_  = $3
set whichSST=$4

if ($#argv < 4 ) then
  echo " "
  echo " Not all inputs are specified. Exiting!"
  echo " "
  exit 1
endif
# ---------------------------------------------------------

# run specific settings

# resolution
set NLAT      = 1440
set NLON      = 2880

# if NLAT=720, NLON=1440, i.e., MERRA-2 BCs resolution, then iMerra can be set to 1.
set iMerra          = 0

# make sure we do not have sea ice if SST > SST_Thr (set below).
set iAdjust_SST_SIC = 1
set SST_Thr         = 283.15
# ---------------------------------------------------------

set one       = 1
set tomorrow_ = `expr ${today_} + $one`

# form the dates
set month     =                `echo ${month_}    | awk '{printf "%02d", $1}'`
set today     = ${year}${month}`echo ${today_}    | awk '{printf "%02d", $1}'`
set tomorrow  = ${year}${month}`echo ${tomorrow_} | awk '{printf "%02d", $1}'`
# ---------------------------------------------------------

# get foundation SST file to scratch dir
set scDir = ${PWD}/${whichSST}/
if (! -e ${scDir}) /bin/mkdir -p ${scDir}

# form source file names and fetch
if (${whichSST} == 'OSTIA') then

  set ostia_data_loc = '/discover/nobackup/dao_ops/intermediate/flk/stage/ostia/'
  set ostia_fPref    = ''
  set ostia_fSuff_   = '-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc.bz2'
  set ostia_fSuff    = '-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc'

  set ostia_file_ops = `/bin/ls ${ostia_data_loc}/${ostia_fPref}${today}${ostia_fSuff_}`
  set fndSST_file    = ${scDir}/${ostia_fPref}${today}${ostia_fSuff}

  if (! -e ${fndSST_file}) then
#   /usr/bin/dmget ${ostia_file_ops}
    /bin/cp ${ostia_file_ops}    ${scDir}
    /usr/bin/bunzip2 ${scDir}/*.bz2
  endif

else
  echo "Unknown option for whichSST:" ${whichSST} ". Exiting"
  exit 99
endif

set reynolds_data_loc = '/gpfsm/dnb42/projects/p17/production/GEOS5odas-5.00/RC/OBS/REYN/2020/'
set reynolds_fPref    = 'oisst-avhrr-v02r01.'
set reynolds_fSuff    = '.nc'
#set reynolds_fSuff2  = '_preliminary.nc'
set reynolds_file     = `/bin/ls ${reynolds_data_loc}/${reynolds_fPref}${today}${reynolds_fSuff}`
# ---------------------------------------------------------

# generate input for processing
if (! -e input.txt) then
  echo ${today}           >input.txt
  echo ${tomorrow}       >>input.txt
  echo ${reynolds_file}  >>input.txt
  echo ${fndSST_file}    >>input.txt
  echo ${NLAT}           >>input.txt
  echo ${NLON}           >>input.txt
  echo ${iMerra}         >>input.txt
  echo ${iAdjust_SST_SIC}>>input.txt
  echo ${SST_Thr}        >>input.txt
else
  echo "input.txt already exists. Nothing to do. Exiting!"
endif
