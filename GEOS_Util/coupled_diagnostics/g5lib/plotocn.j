#!/bin/csh 
#SBATCH --job-name=@PLOT_T
#SBATCH --ntasks=20
#SBATCH --time=12:00:00 
#SBATCH -o plotocn.out

if ( $HOSTNAME =~ "discover"* ) then
    source /usr/share/modules/init/csh
    module purge
    module load other/comp/gcc-5.3-sp3
    module load other/SSSO_Ana-PyD/SApd_4.2.0_py2.7_gcc-5.3-sp3
    set SCRDIR=/home/yvikhlia
    set G5LIBDIR=$SCRDIR/python
    set ANADIR=$SCRDIR/geos5/ana
    set VERIFICATION=$SCRDIR/verification
    setenv OCEANVAL /discover/nobackup/projects/gmao/oceanval/verification
else if ( $HOSTNAME =~ "pfe"* ) then
    module purge
    module load python/GEOSpyD/Ana2018.12_py2.7
    set SCRDIR=/nobackup/gmao_SIteam/Utilities/ocean-scripting
    set G5LIBDIR=$SCRDIR/python
    set ANADIR=$SCRDIR/geos5/ana
    set VERIFICATION=$SCRDIR/verification
    setenv OCEANVAL /nobackup/gmao_SIteam/ModelData/oceanval/verificarion
endif 

set EXPID=@EXPID
set HOMDIR=@HOMDIR
set EXPDIR=@EXPDIR
set WORKDIR=$EXPDIR/plot

setenv PYTHONPATH ${G5LIBDIR}:${VERIFICATION}:${HOMDIR}/..

if (! -e $WORKDIR  ) mkdir -p $WORKDIR
cd $WORKDIR 

sed -e "s|@expid|${EXPID}|g" -e "s|@anadir|${ANADIR}|g" $G5LIBDIR/g5lib/plotocn > plotocn

if ( $HOSTNAME =~ "discover"* ) then
    if($?PBS_JOBID) then
        /usr/local/other/pods/pods.sh $WORKDIR/plotocn
    else
        /bin/csh $WORKDIR/plotocn
    endif
else
    /bin/csh $WORKDIR/plotocn
endif

exit 0 
