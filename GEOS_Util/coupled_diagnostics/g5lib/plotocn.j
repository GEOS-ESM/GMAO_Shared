#!/bin/csh 
#SBATCH --job-name=@PLOT_T
#SBATCH --ntasks=20
#SBATCH --time=12:00:00 
#SBATCH -o plotocn.out

source /usr/share/modules/init/csh
module purge
module load python/GEOSpyD/Ana2019.10_py2.7
set SCRDIR=@GEOSDIR/coupled_diagnostics
set G5LIBDIR=$SCRDIR
set ANADIR=$SCRDIR/analysis
set VERIFICATION=$SCRDIR/verification

if (($HOSTNAME =~ discover*) || ($HOSTNAME =~ borg*)) then
    setenv OCEANVAL /discover/nobackup/projects/gmao/oceanval/verification
else if ( $HOSTNAME =~ "pfe"* ) then
    setenv OCEANVAL /nobackup/gmao_SIteam/ModelData/oceanval/verification
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
