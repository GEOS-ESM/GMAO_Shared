#!/bin/csh
    if ( $#argv < 3 ) goto usage

    setenv OPENGRADS $SHARE/dasilva/opengrads
     setenv OPENGRADS_VERSION 1.9-rc1-gmao
    set EXEC_DIR = /home/ravi/source_motel/GEOSdas-2_1_2-osse-01/GEOSadas/Linux/bin
    set OUT_FILE = ec2fv.ana.eta.nc4
    set RES = d

#   ----------------------------------------------------------
#      Extract User give arguments.
#   ----------------------------------------------------------

   foreach token ( $argv )
     
     if ( "$token" == "-o" ) then
      shift
      if ( $#argv < 1 ) goto usage
       set OUT_FILE = $1
       shift
     else if ( "$token" == "-res" ) then
      shift
      if ( $#argv < 1 ) goto usage
       set RES = $1
       shift
     endif
   end

    if ( $#argv < 3 ) goto usage
     set prog3d_file    = $1
     set prog2d_file    = $2
     set geos5_eta_file = $3
    endif

    if ( -e ${prog3d_file} ) then
      set prog3d_ctlfile = `basename $prog3d_file`
      /home/dasilva/bin/grib2ctl.pl -verf -no_suffix ${prog3d_file} > $prog3d_ctlfile.ctl

      if ( -e ${prog3d_ctlfile} ) then
       $SHARE/dasilva/opengrads/Linux-1.9-rc1-gmao/gribmap -i $prog3d_ctlfile.ctl
       $SHARE/dasilva/opengrads/Linux-1.9-rc1-gmao/gribmap -i $prog3d_ctlfile.ctl

      else
        echo " ----------------------------------------"
        echo "    $prog3d_ctlfile DOES NOT EXIST       "
        echo "    Something wrong with the grib3ctl.pl "
        echo " ----------------------------------------"
          exit
      endif
     else
        echo " ----------------------------------------"
        echo "    $prog3d_file DOES NOT EXIST          "
        echo " ----------------------------------------"
          exit
      endif

      if ( $RES == a ) then
       set HGRID = '72x46'
      else if (  $RES == b ) then
       set HGRID = '144x91'
      else if (  $RES == c ) then
       set HGRID = '288x181'
      else if (  $RES == d ) then
       set HGRID = '540x361'
      endif

    if ( -e ${prog2d_file} ) then
      set prog2d_ctlfile = `basename $prog2d_file`
      /home/dasilva/bin/grib2ctl.pl -verf -no_suffix ${prog2d_file} > $prog2d_ctlfile.ctl

      if ( -e ${prog2d_ctlfile} ) then
       $SHARE/dasilva/opengrads/Linux-1.9-rc1-gmao/gribmap -i $prog2d_ctlfile.ctl
       $SHARE/dasilva/opengrads/Linux-1.9-rc1-gmao/gribmap -i $prog2d_ctlfile.ctl
      else
        echo " ----------------------------------------"
        echo "    $prog2d_ctlfile DOES NOT EXIST       "
        echo "    Something wrong with the grib2ctl.pl "
        echo " ----------------------------------------"
          exit
      endif
    else
        echo " ----------------------------------------"
        echo "    $prog2d_file DOES NOT EXIST          "
        echo " ----------------------------------------"
          exit
    endif

    $SHARE/dasilva/opengrads/Linux-1.9-rc1-gmao/lats4d -o prog.eta -i $prog3d_ctlfile.ctl -vars u v t q o3
    $SHARE/dasilva/opengrads/Linux-1.9-rc1-gmao/lats4d -o prog.sfc -i $prog2d_ctlfile.ctl

    if ( -e prog.eta.nc ) then
      time $EXEC_DIR/GFIO_remap.x -res $RES -o prog.${HGRID}.eta.nc4 prog.eta.nc
    else
        echo " ----------------------------------------"
        echo "     prog.eta.nc DOES NOT EXIST          "
                   ls -lat prog.${HGRID}.eta.nc             
        echo "    Something wrong with the GFIO_remap  "
        echo " ----------------------------------------"
          exit
    endif
      
    if ( -e prog.sfc.nc ) then
      time $EXEC_DIR/GFIO_remap.x -res $RES -o prog.${HGRID}.sfc.nc4 prog.sfc.nc
    else
        echo " ----------------------------------------"
        echo "     prog.sfc.nc DOES NOT EXIST          "
                   ls -lat prog.${HGRID}.sfc.nc             
        echo "    Something wrong with the GFIO_remap  "
        echo " ----------------------------------------"
          exit
    endif

    if ( -e prog.${HGRID}.eta.nc4 && -e prog.${HGRID}.sfc.nc4 && -e  $geos5_eta_file ) then
     $EXEC_DIR/ec2fv.x   -o $OUT_FILE -anafile $geos5_eta_file -uafile $geos5_eta_file -phis prog.${HGRID}.sfc.nc4 prog.${HGRID}.eta.nc4
     if ( $status == 0 ) then
#      /bin/rm prog.${HGRID}.sfc.nc4 prog.${HGRID}.eta.nc4 *.nc *.ctl *.pdef *.idx
    else
        echo " ------------------------------------------------"
        echo "     prog.540x360.eta.nc4, prog.${HGRID}.sfc.nc4  "
        echo "     and $geos5_eta_file                         "
        echo "     DOES NOT EXIST.                             "
        echo "    Something wrong with the GFIO_remap.         "
                   ls -lat prog.${HGRID}.eta.nc4         
                   ls -lat prog.${HGRID}.sfc.nc4         
                   ls -lat $geos5_eta_file 
        echo " ------------------------------------------------"
          exit
    endif
 exit

usage:
clear
cat << ---//---

NAME

         ec2fv.csh -  Maps  Grib Nature 91 level model file
                      to GEOS-5 eta file.
        
SYNOPSIS

         ec2fv.csh   [-o] [-res] ecmwf_prog3d_file ecmwf_prog2d_file geos5_eta_file

DESCRIPTION


           Required scripts  get_leap,get_julday
           executable : ec2fv.x
           
          get_leap     :    Checks the year for leap year.
          get_julday   :    Converts the given date to day of the year.
          ec2fv.x      :    Reads Nature run model level prog and surface hdf files
	                    converts to fv mapped bkg.eta files.

OPTIONS

          -o           Output file ( default ec2fv.nature.ana.eta)
          -res         Horizontal resolution [ a        72x46
                                               b       144x91
                                               c       288x181
                                               d       540x361   default]

Example:

          ec2fv.csh -o OUT_FILE ecmwf_nature_prog_file ecmwf_nature_prog_sfc_file geos05_ana_eta_file
----//---
