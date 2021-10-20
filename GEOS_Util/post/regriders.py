#!/usr/bin/env python
#
import os
import subprocess
import shutil
import glob

tags_to_bcs ={"Icarus_Reynolds": "/discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/bcs/Icarus_Updated/Icarus_Reynolds/", 
              "Icarus-NLv3_MERRA-2": "/discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/bcs/Icarus-NLv3/Icarus-NLv3_MERRA-2/"}

def grid_to_string(grid):
  if(grid[0] == 'C'):
    n = int(grid[1:])
    return f'CF{n:04d}x6C' 
  else:
    xy = grid.upper().split('X')
    x = int(xy[0])
    y = int(xy[1])
    if (y/x == 6) :
      return f'CF{x:04d}x6C'
    else:
      return f'DE{x:04d}xPE{y:04d}'

class regrider:
  def __init__(self, config):
     self.common_in  = config['options_in']['COMMON']
     self.common_out = config['options_out']['COMMON']
     self.restarts_in = config['restarts_in']
     self.slurm_options = config['slurm_options']

     # deduce the options
     agrid_str = grid_to_string(self.common_in['agrid'])
     ogrid_str = grid_to_string(self.common_in['ogrid'])
     ao_str = agrid_str + '_' + ogrid_str
     self.in_bcsdir = tags_to_bcs[self.common_in['tag']] + ao_str
     self.in_til    = self.in_bcsdir + '/'+ao_str + '-Pfafstetter.til' 
     
     agrid_str = grid_to_string(self.common_out['agrid'])
     ogrid_str = grid_to_string(self.common_out['ogrid'])
     ao_str = agrid_str + '_' + ogrid_str
     self.out_bcsdir = tags_to_bcs[self.common_out['tag']] + ao_str
     self.out_til    = self.out_bcsdir + '/'+ao_str + '-Pfafstetter.til' 

class upperair(regrider):
  def __init__(self, config):
     super().__init__(config)
  def regrid(self):
     ...

class surface(regrider):
  def __init__(self, config):
     super().__init__(config)
     self.restarts_in = self.restarts_in['SURFACE']
     self.surf_in  = config['options_in']['SURFACE']
     self.surf_out = config['options_out']['SURFACE']

     print("--ToDO--\n")
     print("split saltwater after some tags(need to rank tag)")
     print("if split, set self.saltwater None")
     print("\n")
     self.saltwater = self.restarts_in.get('saltwater')
     self.openwater = self.restarts_in.get('openwater')
     self.seaice  = self.restarts_in.get('seaice')
     self.landice = self.restarts_in.get('landice')
     self.lake  = self.restarts_in.get('lake')
     self.route = self.restarts_in.get('route')
     self.catch = self.restarts_in.get('catch')
     self.catchcnclm40 = self.restarts_in.get('catchcnclm40')
     self.catchcnclm45 = self.restarts_in.get('catchcnclm45')
         
  def regrid(self):
     bindir = os.getcwd()
     outdir = self.common_out['outdir']
     if not os.path.exists(outdir) : os.makedirs(outdir)
     os.chdir(self.common_out['outdir'])
     InData_dir = outdir+'/InData/'
     OutData_dir = outdir+'/OutData/'
     if os.path.exists(InData_dir) : subprocess.call('rm -rf '+ InData_dir, shell = True)
     os.makedirs(InData_dir)
     if os.path.exists(OutData_dir) : subprocess.call('rm -rf '+ OutData_dir, shell = True)
     os.makedirs(OutData_dir)

     for key, value in self.restarts_in.items():
        if not value : continue
        src  = self.common_in['rstdir']+'/'+value
        if not os.path.exists(src) :
           print( "Could not find " + src + "--remove it from yaml config file add to restart dir")
           exit()
        dest = InData_dir+'/'+value
        if os.path.exists(dest) : os.unlink(dest)
        os.symlink(src, dest)
     in_til = InData_dir+'/' + self.in_til.split('/')[-1]
     out_til = OutData_dir+'/'+self.out_til.split('/')[-1]
     if os.path.exists(in_til)  : os.unlink(in_til)
     if os.path.exists(out_til) : os.unlink(out_til)
     os.symlink(self.in_til, in_til)
     os.symlink(self.out_til, out_til)

     exe = bindir + '/mk_LakeLandiceSaltRestarts '
     if (self.saltwater):
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.saltwater + ' 0 ' + str(self.surf_in['zoom'])
       print(cmd)
       subprocess.call(cmd, shell= True)

     if (self.openwater):
       exe = bindir + '/mk_LakeLandiceSaltRestarts '
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.openwater + ' 0 ' + str(self.surf_in['zoom'])
       print(cmd)
       subprocess.call(cmd, shell= True)

     if (self.seaice):
       exe = bindir + '/mk_LakeLandiceSaltRestarts '
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.seaice + ' 0 ' + str(self.surf_in['zoom'])
       print(cmd)
       subprocess.call(cmd, shell= True)

     if (self.lake):
       exe = bindir + '/mk_LakeLandiceSaltRestarts '
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.lake + ' 19 ' + str(self.surf_in['zoom'])
       print(cmd)
       subprocess.call(cmd, shell= True)

     if (self.landice):
       exe = bindir + '/mk_LakeLandiceSaltRestarts '
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.landice + ' 20 ' + str(self.surf_in['zoom'])
       print(cmd)
       subprocess.call(cmd, shell= True)

     if (self.route):
       route = bindir + '/mk_RouteRestarts '
       cmd = route + 'OutData/*.til '+ self.common_in["yyyymmddhh"][0:6]
       print(cmd)
       subprocess.call(cmd, shell= True)
     if ( self.catchcnclm40 or self.catchcnclm45) :
       dirname = os.path.dirname(self.out_til)
       clsm = dirname+'/clsm'
       os.symlink(clsm, 'OutData/clsm')

     mk_catch_j_template = """#!/bin/csh -f
#SBATCH --account={account}
#SBATCH --ntasks=84
#SBATCH --time=1:00:00
#SBATCH --job-name=mk_catch
#SBATCH --qos=debug
#SBATCH --output={outdir}/{mk_catch_log}
#

source {Bin}/g5_modules
set echo

#limit stacksize unlimited
unlimit

set esma_mpirun_X = ( {Bin}/esma_mpirun -np 84 )
set mk_CatchRestarts_X   = ( $esma_mpirun_X {Bin}/mk_CatchRestarts )
set mk_CatchCNRestarts_X = ( $esma_mpirun_X {Bin}/mk_CatchCNRestarts )
set Scale_Catch_X   = {Bin}/Scale_Catch
set Scale_CatchCN_X = {Bin}/Scale_CatchCN

set OUT_til   = OutData/*.til
set IN_til    = InData/*.til

if ({catchFLG}) then
    set catchIN = InData/*catch_internal_rst*
    set params = ( $OUT_til $IN_til $catchIN {surflay} )
    $mk_CatchRestarts_X $params

    if ({rescale}) then
        set catch_regrid = OutData/$catchIN:t
        set catch_scaled = $catch_regrid.scaled
        set params = ( $catchIN $catch_regrid $catch_scaled {surflay} )
        set params = ( $params {wemin} {wemout} )
        $Scale_Catch_X $params

        mv $catch_regrid $catch_regrid.1
        mv $catch_scaled $catch_regrid
    endif
endif

if ({catchcnFLG}) then
    # WY notes: it was never used in gcm
    if ({fromGCM} ) then
       set catchcnIN = InData/*catchcn_internal_rst*
       set params = ( $OUT_til $IN_til $catchcnIN {surflay} {rsttime} )
       $mk_CatchCNRestarts_X $params
    else # from GEOSldas
       set OUT_til   = `ls OutData/*.til | cut -d '/' -f2`
       /bin/cp OutData/*.til OutData/OutTileFile
       /bin/cp OutData/*.til InData/OutTileFile
       set RESTART_short = {RESTART_PATH}/{RESTART_ID}/output/{RESTART_DOMAIN}/
       set YYYY = `echo {rsttime} | cut -c1-4`
       set MM   = `echo {rsttime} | cut -c5-6`
       echo $MM, $YYYY
       set PARAM_FILE = `ls $RESTART_short/rc_out/Y$YYYY/M$MM/*ldas_catparam* | head -1`
       set params = ( -b OutData/ -d {rsttime} -e {RESTART_ID} -m catchcn{CN_VERSION} -s {surflay} -j Y -r R -p $PARAM_FILE -l $RESTART_short)
       $mk_GEOSldasRestarts_X $params
    endif

    if ({rescale}) then
        set catchcn_regrid = OutData/$catchcnIN:t
        set catchcn_scaled = $catchcn_regrid.scaled
        set params = ( $catchcnIN $catchcn_regrid $catchcn_scaled {surflay} )
        set params = ( $params {wemin} {wemout} )
        $Scale_CatchCN_X $params

        mv $catchcn_regrid $catchcn_regrid.1
        mv $catchcn_scaled $catchcn_regrid
    endif
endif
exit
"""
     catchFLG = '1' if self.catch else '0'
     catchcnFLG = '1' if (self.catchcnclm40 or self.catchcnclm45) else '0'

     catch1script =  mk_catch_j_template.format(Bin = bindir, account = self.slurm_options['account'], \
                  outdir = outdir, mk_catch_log = 'mk_catch_log.1', surflay = self.surf_in['surflay'],  \
                  wemin = self.surf_in['wemin'], wemout = self.surf_out['wemout'] ,  \
                  fromGCM = '0', catchFLG = catchFLG, catchcnFLG = catchcnFLG, rescale = '0', rsttime = self.common_in['yyyymmddhh'], \
                  RESTART_ID = self.surf_in['restart_id'], RESTART_PATH = self.surf_in['restart_path'], \
                  RESTART_DOMAIN = self.surf_in['restart_domain'], CN_VERSION = self.surf_in['cn_version'] ) 
     print(os.getcwd())
     catch1 = open('mk_catch.j.1','wt')
     catch1.write(catch1script)
     catch1.close()

     subprocess.call('sbatch -W mk_catch.j.1', shell= True)
     # step 2
     if (self.surf_in['rescale'] and ( self.catch or self.catchcnclm40 or self.catchcnclm45)) :
       if os.path.exists('InData.step1') : subprocess.call('rm -rf InData.step1', shell = True)
       shutil.move('InData', 'InData.step1')
       os.makedirs('InData')
       for catchfile in glob.glob("OutData/*catch*"):
         shutil.move(catchfile,"InData/")
       os.symlink(self.out_til, in_til)
       
       if ( (not (self.catchcnclm40 or self.catchcnclm45)) and self.catch) :
          dirname = os.path.dirname(self.out_til)
          clsm = dirname+'/clsm'
          os.symlink(clsm, 'OutData/clsm')

       catch2script =  mk_catch_j_template.format(Bin = bindir, account = self.slurm_options['account'], \
                  outdir = outdir, mk_catch_log = 'mk_catch_log.2', surflay = self.surf_in['surflay'],  \
                  wemin = self.surf_in['wemin'], wemout = self.surf_out['wemout'] ,  \
                  fromGCM = '0', catchFLG = catchFLG, catchcnFLG = catchcnFLG, rescale = '1', rsttime = self.common_in['yyyymmddhh'], \
                  RESTART_ID = self.surf_in['restart_id'], RESTART_PATH = self.surf_in['restart_path'], \
                  RESTART_DOMAIN = self.surf_in['restart_domain'], CN_VERSION = self.surf_in['cn_version'] )
     print(os.getcwd())
     catch2 = open('mk_catch.j.2','wt')
     catch2.write(catch2script)
     catch2.close()

     subprocess.call('sbatch -W mk_catch.j.2', shell= True)
     