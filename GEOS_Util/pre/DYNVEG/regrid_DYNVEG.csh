#!/bin/csh -f

# check whether GEOSldas or GCM
###############################

if( -e LDAS.rc ) set JOB_FILE = lenkf.j
if( -e AGCM.rc ) set JOB_FILE = gcm_run.j
    
set EXPID = `grep EXPID $JOB_FILE | head -1 | rev | cut -d' ' -f1 | rev`
set LineNo1 = `grep -n -m 1 "EXPDIR" $JOB_FILE | cut -d':' -f1`
set EXPDIR = `grep EXPDIR $JOB_FILE | head -1 | rev | cut -d' ' -f1 | rev`
set EXPDIR = `echo $EXPDIR | sed 's|$EXPID|'"$EXPID"'|g'`
set HOMDIR = `grep HOMDIR $JOB_FILE | head -1 | rev | cut -d' ' -f1 | rev`
set HOMDIR = `echo $HOMDIR | sed 's|$EXPDIR|'"$EXPDIR"'|g'`
set SPONSORID = `getsponsor | tail -3 | head -1 | cut -d'|' -f1`

if( -e LDAS.rc ) then
    set GEOSBIN = $EXPDIR/build/bin/
    set MODELING_SYSTEM = GEOSldas
    set BCSPATH = `ls -d $EXPDIR/output/*/rc_out/`
    set  BCSDIR = `grep BCS_PATH $EXPDIR/run/* | cut -d':' -f3`
    set TILFILE = `ls $EXPDIR/output/*/rc_out/*tilecoord.bin`
    set PRESCRIBE_DVG = `grep PRESCRIBE_DVG LDAS.rc | cut -d':' -f2`
endif

if( -e AGCM.rc ) then
    set GEOSBIN = `grep GEOSBIN $JOB_FILE | head -1 | rev | cut -d' ' -f2 | rev`
    set MODELING_SYSTEM = GCM
    set BCSDIR = `grep BCSDIR $JOB_FILE | head -1 | rev | cut -d' ' -f1 | rev`
    set AGCM_IM = `grep AGCM_IM AGCM.rc | head -1 | cut -d':' -f2`
    set AGCM_JM = `grep AGCM_JM AGCM.rc | head -1 | cut -d':' -f2`
    set OGCM_IM = `grep OGCM.IM_WORLD AGCM.rc | head -1 | cut -d':' -f2`
    set OGCM_JM = `grep OGCM.JM_WORLD AGCM.rc | head -1 | cut -d':' -f2` 
    set BCSPATH = `echo $BCSDIR/*${AGCM_IM}*${OGCM_IM}*${OGCM_JM}`
    set BCSDIR  = $BCSPATH
    set TILFILE = `ls $BCSPATH/*til | head -1 | rev | cut -d'/' -f1 | rev`
    set PRESCRIBE_DVG = `grep PRESCRIBE_DVG RC/GEOS_LandGridComp.rc | cut -d':' -f2`
endif

if( $PRESCRIBE_DVG == 1 ) set ACTUAL_DV = 'TRUE'   
if( $PRESCRIBE_DVG >= 2 ) set ACTUAL_DV = 'FALSE'  

/bin/mkdir -p ${EXPDIR}/VEGDATA/scratch
set PWD = `pwd`
set YEARB = `head -1 $HOMDIR/cap_restart | cut -c1-4`
set YEARE = `head -n 15 $HOMDIR/CAP.rc | grep END_DATE | cut -d':' -f2 | cut -d' ' -f2 | cut -c1-4`
@ YEARB--
@ YEARE++

cd ${EXPDIR}/VEGDATA/scratch

if($ACTUAL_DV == TRUE) then

    while ($YEARB <= $YEARE)
    
	tar -xvzf /discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/DYNVEG/CF0180x6C_TM0720xTM0410/CF180_DVG2_C4.tavg1_tile_veg.${YEARB}.tar.gz
	/bin/rm *monthly*
	
	@ YEARB++
	
    end

    set YEARB = `head -1 $HOMDIR/cap_restart | cut -c1-4`
    set YEARE = `head -n 15 $HOMDIR/CAP.rc | grep END_DATE | cut -d':' -f2 | cut -d' ' -f2 | cut -c1-4`

else
    tar -xvzf /discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/DYNVEG/CF0180x6C_TM0720xTM0410/CF180_DVG2_C4.tavg1_tile_veg.YYYY.tar.gz
endif

cat << _EOF_ > ${EXPDIR}/VEGDATA/scratch/dynveg_run.j
#!/bin/csh -f

#PBS -l walltime=12:00:00
#SBATCH --ntasks=56
#SBATCH --constraint=hasw
#SBATCH --account=$SPONSORID 
#PBS -N dynveg_RUN

setenv EXPID $EXPID
setenv EXPDIR $EXPDIR
setenv HOMDIR $HOMDIR
setenv GEOSBIN $GEOSBIN
source $GEOSBIN/g5_modules
setenv MODELING_SYSTEM $MODELING_SYSTEM
setenv BCSDIR $BCSDIR
setenv BCSPATH $BCSPATH
setenv TILFILE $TILFILE
setenv OMPI_MCA_shmem_mmap_enable_nfs_warning 0
setenv YEARB $YEARB
setenv YEARE $YEARE
setenv ACTUAL_DV $ACTUAL_DV

limit stacksize unlimited

cd ${EXPDIR}/VEGDATA/scratch

##mpirun -map-by core --mca btl ^vader -np 56 $GEOSBIN/regrid_DYNVEG.x
mpirun -np 56 $GEOSBIN/regrid_DYNVEG.x
sleep 300
/bin/rm *.nc *.nc4

_EOF_

chmod 755 ${EXPDIR}/VEGDATA/scratch/dynveg_run.j
cd ${EXPDIR}/VEGDATA/scratch 
sbatch dynveg_run.j

exit



