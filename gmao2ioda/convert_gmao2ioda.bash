#!/bin/bash

####################################################################################
# This script convert gmao ocean observation to IODA/JEDI format for each platform #
# usage: ./convert_gmao2ioda.bash YYYYMMDDHH                                       #
####################################################################################


# load multiple modules (may updated by Doruk 
source py_modules_310_li.sh
source ocean_das_config
idate=$1
if [ $# -ne 1 ]; then
  echo "e.g., ./convert_gmao2ioda.bash 2022051600"
  exit 1
fi

# input format '+%Y%m%d%H`

year=$(echo $idate | cut -c1-4)
month=$(echo $idate | cut -c5-6)
day=$(echo $idate | cut -c7-8)
hour=$(echo $idate | cut -c9-10)

# print the extracted values
echo "Year:  $year"
echo "Month: $month"
echo "Day:   $day"
echo "Hour:  $hour"

## para
scripts='/discover/nobackup/lren1/jedi_obs/scripts'
ioda_data='/discover/nobackup/lren1/jedi_obs/ioda_ocean_obs'
gmao_data='/discover/nobackup/lren1/jedi_obs/gmao_ocean_obs'
     

platforms=("Argo" "CTD" "XBT" "TAO" "PIRATA" "RAMA" "Jason-1" "Jason-2" "Jason-3" "Saral" "Sentinel-3a" "Sentinel-6a" "Sentinel-3b" "ERS-1" "ERS-2" "TOPEX" "GEOSAT-2" "Envisat" "HY-2A" "CryoSat-2" "TOPEX-N" "Envisat-N" "Jason-1-N" "Jason-1G" "Jason-2-N" "Jason-3N" "CryoSat-2-N" "SWOTN" "SMOS" "AQUARIUS" "SMAPV5.0")

platforms=('Argo')
cd $scripts
for platform in "${platforms[@]}"; do

    echo '############'
    echo pf: $platform
# extracing gmao ocean observation in every 6-hour 
    ./ocean_obs_4ioda.py $year $month $day $hour $platform
    rm *.png

# converting the 6-hour gmao ocean observation to IODA/JEDI format
    InputName=gmao-obs-$idate-$platform
    OutName=ioda-obs-$idate-$platform
    ifn=$gmao_data/$InputName.nc
    ofn=$ioda_data/$OutName.nc
    echo "./gmao_obs2ioda.py -i $ifn -o $ofn"
    ./gmao_obs2ioda.py -i $ifn -o $ofn

done 

 

