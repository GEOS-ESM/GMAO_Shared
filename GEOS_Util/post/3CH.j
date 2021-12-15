#!/bin/csh -f

unsetenv LD_LIBRARY_PATH

setenv  ARCH `uname`
setenv  GEOSBIN /discover/nobackup/ltakacs/TAGS/Jason-4_0/GEOSagcm/Linux/bin
source $GEOSBIN/g5_modules
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BASEDIR}/${ARCH}/lib

setenv RUN_CMD "$GEOSBIN/esma_mpirun -np "

if ($?SLURM_NTASKS) then
        set NCPUS = $SLURM_NTASKS
else
        set NCPUS = 1
endif

# ----------------------------------------------------
# Define Experiment IDs, Locations, and Relevant Files
# ----------------------------------------------------

set nymdb = `cat 3CH.rc | grep -v \# | grep nymdb: | cut -d: -f2`
set nhmsb = `cat 3CH.rc | grep -v \# | grep nhmsb: | cut -d: -f2`
set nymde = `cat 3CH.rc | grep -v \# | grep nymde: | cut -d: -f2`
set nhmse = `cat 3CH.rc | grep -v \# | grep nhmse: | cut -d: -f2`
set ndt   = `cat 3CH.rc | grep -v \# | grep   ndt: | cut -d: -f2`

if( $ndt == '' ) set ndt = 21600

@ numexps = `cat 3CH.rc | grep -v \# | grep expid_ | wc -l`
set command_line = " "

@ n = 1
while( $n <= $numexps )
set expid = `cat 3CH.rc | grep -v \# | grep expid_$n | cut -d: -f2`
set files = `cat 3CH.rc | grep -v \# | grep files_$n | cut -d: -f2`

/bin/rm -f exp_files.$n
touch      exp_files.$n

set command_line = "$command_line -exp"$n" $expid "

set nymd = $nymdb
set nhms = 000000
while( $nymd <= $nymde )
        set year  = `echo $nymd | cut -b1-4`
        set month = `echo $nymd | cut -b5-6`
        set day   = `echo $nymd | cut -b7-8`
        set hour  = `echo $nhms | cut -b1-2`
        set dir = `echo $files | sed -e "s/%y4/$year/g" | sed -e "s/%m2/$month/g" | sed -e "s/%d2/$day/g" | sed -e "s/%h2/$hour/g"`
        echo Looking for $dir ...
        set list = `/bin/ls -1 ${dir}* | grep -v daily_mean`
        echo $list >> exp_files.$n
        set date = `tick $nymd $nhms $ndt`
        set nymd = $date[1]
        set nhms = $date[2]
end
set dum = `cat exp_files.$n`

set command_line = "$command_line $dum "

echo " "
@ n = $n + 1
end

/bin/rm -f exp_files.*

set im     = `cat 3CH.rc | grep -v \# | grep  im:     | cut -d: -f2` ; if( "$im"     != '' ) set im     = "-im     $im"
set jm     = `cat 3CH.rc | grep -v \# | grep  jm:     | cut -d: -f2` ; if( "$jm"     != '' ) set jm     = "-jm     $jm"
set levs   = `cat 3CH.rc | grep -v \# | grep  levs:   | cut -d: -f2` ; if( "$levs"   != '' ) set levs   = "-levs   $levs"
set fields = `cat 3CH.rc | grep -v \# | grep  fields: | cut -d: -f2` ; if( "$fields" != '' ) set fields = "-fields $fields"

set method = `cat 3CH.rc | grep -v \# | grep  interpolation_method: | cut -d: -f2`

if($method != '' ) then
   if( $method == 'BL' ) set method = "-method 1"
   if( $method == 'BC' ) set method = "-method 2"
   if( $method == 'BA' ) set method = "-method 3"
endif

# ----------------------------------------------------
# Perform 3CH Calculation
# ----------------------------------------------------

$RUN_CMD $NCPUS $GEOSBIN/3CH_mpi.x $command_line $im $jm $levs $fields $method \
                                   -nymdb $nymdb -nhmsb $nhmsb -nymde $nymde -nhmse $nhmse -ndt $ndt  

# ----------------------------------------------------
# Convert grads to nc4 file
# ----------------------------------------------------

set ctl_files = `/bin/ls -1t *ctl`
set ctl_file  = $ctl_files[1]
set length    = `echo $ctl_file | awk '{print length($0)}'`
@ loc = $length - 4
set  newfile = `echo $ctl_file | cut -b1-$loc`
echo Newfile = $newfile

$GEOSBIN/flat2hdf.x -flat $newfile.data -ctl $newfile.ctl -nymd $nymdb -nhms $nhmsb -ndt $ndt

