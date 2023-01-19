#!/usr/bin/env tcsh

# use -db flag to keep outdir from being deleted at end of run
#-------------------------------------------------------------
set flag = ""
if ($#argv > 0) then
    set flag = $1
endif

# run each test found
#--------------------
@ stat = 0
foreach file (test*.py)
    set echo
    $file -v $flag
    @ stat += $status
    unset echo
    sleep 2
end

# check status
#-------------
if ($stat > 0) then
    set msg = "FAILURE"
else
    set msg = "SUCCESS"
endif

echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "status = $stat"
echo $msg
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo
