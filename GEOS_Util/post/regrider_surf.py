#!/usr/bin/env python3
#
import os
import subprocess
import shutil
import glob
from regrider_base import *

class surface(regrider):
  def __init__(self, config):
     #v3.0
     #super().__init__(config)
     #v2.7
     super(surface, self).__init__(config)
     self.restarts_in = self.restarts_in['SURFACE']
     self.surf_in  = config['input']['parameters']['SURFACE']
     self.surf_out = config['output']['parameters']['SURFACE']

     self.saltwater = None
     self.openwater = None
     self.seaice  = None
     self.landice = None
     self.lake  = None
     self.route = None
     self.catch = None
     self.catchcnclm40 = None
     self.catchcnclm45 = None
     for rst in self.restarts_in :
        f = os.path.basename(rst)
        if (f.find('saltwater') != -1):
          self.saltwater = f
          continue
        if (f.find('seaice') != -1):
          self.seaice = f
          continue
        if (f.find('landice') != -1):
          self.landice = f
          continue
        if (f.find('lake') != -1):
          self.lake = f
          continue
        if (f.find('route') != -1):
          self.route = f
          continue
        if (f.find('catch_') != -1):
          self.catch = f
          continue
        if (f.find('catchcnclm40') != -1):
          self.catchcnclm40 = f
          continue
        if (f.find('catchcnclm45') != -1):
          self.catchcnclm45 = f
          continue
     tagout = self.common_out['tag']
     ogrid = self.common_out['ogrid']
     bctag = self.get_bcTag(tagout, ogrid)
     tagrank = self.tagsRank[bctag]

     self.surf_out['surflay'] = 20.
     if tagrank >=12 :
       self.surf_out['surflay'] = 50.
     self.surf_in['rescale'] = 0 
     if tagrank > self.tagsRank["Fortuna-2_0"]:
       self.surf_in['rescale'] = 1 
         
  def regrid(self):
     bindir = os.getcwd()
     out_dir = self.common_out['out_dir']
     if not os.path.exists(out_dir) : os.makedirs(out_dir)
     os.chdir(self.common_out['out_dir'])
     InData_dir = out_dir+'/InData/'
     OutData_dir = out_dir+'/OutData/'
     if os.path.exists(InData_dir) : subprocess.call(['rm','-rf', InData_dir])
     os.makedirs(InData_dir)
     if os.path.exists(OutData_dir) : subprocess.call(['rm','-rf', OutData_dir])
     os.makedirs(OutData_dir)
    
     for rst in self.restarts_in:
        f = os.path.basename(rst)
        dest = InData_dir+'/'+f
        if os.path.exists(dest) : shutil.remove(dest)
        print('\nCopy ' + rst + ' to ' +dest)
        shutil.copy(rst,dest)

     in_til = InData_dir+'/' + self.in_til.split('/')[-1]
     out_til = OutData_dir+'/'+self.out_til.split('/')[-1]
     if os.path.exists(in_til)  : shutil.remove(in_til)
     if os.path.exists(out_til) : shutil.remove(out_til)
     print('\n Copy ' + self.in_til + ' to ' + in_til)
     shutil.copy(self.in_til, in_til)
     print('\n Copy ' + self.out_til + ' to ' + out_til)
     shutil.copy(self.out_til, out_til)

     ogrid = self.common_out['ogrid']
     tagout = self.common_out['tag']
     bctag  = self.get_bcTag(tagout, ogrid)
     tagrank = self.tagsRank[bctag]
     split_saltwaterFLG = '0'
     if  tagrank >= self.tagsRank["Icarus_Reynolds"] and self.saltwater:
       split_saltwaterFLG = '1'
      
     exe = bindir + '/mk_LakeLandiceSaltRestarts '
     if (self.saltwater):
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.saltwater + ' 0 ' + str(self.surf_in['zoom'])
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))
  
       # split Saltwater
       ogrid = self.common_out['ogrid']
       tagout = self.common_out['tag']
       bctag  = self.get_bcTag(tagout, ogrid)
       tagrank = self.tagsRank[bctag]
       if  split_saltwaterFLG == '1':
         print("\nSplitting Saltwater...\n")
         cmd = bindir+'/SaltIntSplitter ' + out_til + ' ' + 'OutData/'+self.saltwater
         print('\n'+cmd)
         subprocess.call(shlex.split(cmd))
         self.openwater = None
         self.seaice  = None

     if (self.openwater):
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.openwater + ' 0 ' + str(self.surf_in['zoom'])
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     if (self.seaice):
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.seaice + ' 0 ' + str(self.surf_in['zoom'])
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     if (self.lake):
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.lake + ' 19 ' + str(self.surf_in['zoom'])
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     if (self.landice):
       cmd = exe + 'OutData/*.til InData/*.til InData/'+self.landice + ' 20 ' + str(self.surf_in['zoom'])
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     if (self.route):
       route = bindir + '/mk_RouteRestarts '
       cmd = route + 'OutData/*.til '+ self.common_in["yyyymmddhh"][0:6]
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     if ( self.catchcnclm40 or self.catchcnclm45) :
       dirname = os.path.dirname(self.out_til)
       clsm = dirname+'/clsm'
       print( "\nsymbolic link clsm to OutData/clsm")
       os.symlink(clsm, 'OutData/clsm')

     mk_catch_j_template = """#!/bin/csh -f
#SBATCH --account={account}
#SBATCH --ntasks=84
#SBATCH --time=1:00:00
#SBATCH --job-name=mk_catch
#SBATCH --qos=debug
#SBATCH --output={out_dir}/{mk_catch_log}
#

source {Bin}/g5_modules
set echo

#limit stacksize unlimited
unlimit

#set lakelandRestartX = ( {Bin}/mk_LakeLandiceSaltRestarts )
#set split_saltwaterX = ({Bin}/SaltIntSplitter )
#set routRestartX = ( {Bin}/mk_RouteRestarts )
set esma_mpirun_X = ( {Bin}/esma_mpirun -np 84 )
set mk_CatchRestarts_X   = ( $esma_mpirun_X {Bin}/mk_CatchRestarts )
set mk_CatchCNRestarts_X = ( $esma_mpirun_X {Bin}/mk_CatchCNRestarts )
set Scale_Catch_X   = {Bin}/Scale_Catch
set Scale_CatchCN_X = {Bin}/Scale_CatchCN

set OUT_til   = OutData/*.til
set IN_til    = InData/*.til

#if ( {saltwaterFLG}) then
#   $lakelandRestartX $OUT_til $IN_til InData/{saltwater} 0 {zoom}
#   if ({split_saltwaterFLG}) then
#     $split_saltwaterX $OUT_til OutData/{saltwater}
#   endif
#endif

#if ( {openwaterFLG}) then
#   $lakelandRestartX $OUT_til $IN_til InData/{openwater} 0 {zoom} 
#endif
#
#if ( {seaiceFLG}) then
#   $lakelandRestartX $OUT_til $IN_til InData/{seaice} 0 {zoom} 
#endif
#
#if ( {lakeFLG}) then
#   $lakelandRestartX $OUT_til $IN_til InData/{lake} 19 {zoom} 
#endif
#
#if ( {landiceFLG}) then
#   $lakelandRestartX $OUT_til $IN_til InData/{landice} 20 {zoom} 
#endif
#
#if ( {routeFLG}) then
#   $routeRestartX $OUT_til {yyyymm} 
#endif

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

        rm $catch_regrid
        #mv $catch_regrid $catch_regrid.1
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

        rm $catchcn_regrid
        #mv $catchcn_regrid $catchcn_regrid.1
        mv $catchcn_scaled $catchcn_regrid
    endif
endif
"""
     saltwaterFLG = '1' if self.saltwater else '0'
     openwaterFLG = '1' if self.openwater else '0'
     seaiceFLG    = '1' if self.seaice    else '0'
     landiceFLG   = '1' if self.landice   else '0'
     lakeFLG      = '1' if self.lake      else '0'
     routeFLG     = '1' if self.route     else '0'
     if (split_saltwaterFLG == '1' ):
       openwaterFLG = '0'
       seaiceFLG = '0'

     catchFLG = '1' if self.catch else '0'
     catchcnFLG = '1' if (self.catchcnclm40 or self.catchcnclm45) else '0'

     catch1script =  mk_catch_j_template.format(Bin = bindir, account = self.slurm_options['account'], \
                  out_dir = out_dir, mk_catch_log = 'mk_catch_log.1', surflay = self.surf_out['surflay'],  \
                  wemin = self.surf_in['wemin'], wemout = self.surf_out['wemout'] ,  \
                  fromGCM = '0', catchFLG = catchFLG, catchcnFLG = catchcnFLG, rescale = '0', rsttime = self.common_in['yyyymmddhh'], \
                  RESTART_ID = self.surf_in.get('restart_id'), RESTART_PATH = self.surf_in.get('restart_path'), \
                  RESTART_DOMAIN = self.surf_in.get('restart_domain'), CN_VERSION = self.surf_in.get('cn_version'), \
                  saltwaterFLG = saltwaterFLG,  saltwater = self.saltwater,  \
                  openwaterFLG = openwaterFLG,  openwater = self.openwater,  \
                  seaiceFLG    = seaiceFLG,     seaice    = self.seaice,  \
                  landiceFLG   = landiceFLG,    landice   = self.landice,  \
                  lakeFLG      = lakeFLG,       lake      = self.lake,   \
                  split_saltwaterFLG = split_saltwaterFLG, \
                  routeFLG     = routeFLG,      route     = self.route, yyyymm = str(self.common_in['yyyymmddhh'])[0:6],  \
                  zoom = self.surf_in['zoom'] ) 
     catch1 = open('mk_catch.j.1','wt')
     catch1.write(catch1script)
     catch1.close()
     print("step 1: sbatch -W mk_catch.j.1")
     subprocess.call(['sbatch','-W', 'mk_catch.j.1'])
     # step 2
     if (self.surf_in['rescale'] and ( self.catch or self.catchcnclm40 or self.catchcnclm45)) :
       if os.path.exists('InData.step1') : subprocess.call(['rm','-rf', 'InData.step1'])
       print("\n Move Indata to InData.step1")
       shutil.move('InData', 'InData.step1')
       os.makedirs('InData')
       for catchfile in glob.glob("OutData/*catch*"):
         print('\n Move ' + catchfile + ' to InData/')
         shutil.move(catchfile,"InData/")
       print('\n Link ' + self.out_til + ' to ' + in_til)
       os.symlink(self.out_til, in_til)
       
       if ( (not (self.catchcnclm40 or self.catchcnclm45)) and self.catch) :
          dirname = os.path.dirname(self.out_til)
          clsm = dirname+'/clsm'
          print('\n Link ' + clsm + ' to ' + 'OutData/clsm')
          os.symlink(clsm, 'OutData/clsm')

       saltwaterFLG =  '0'
       openwaterFLG =  '0'
       seaiceFLG    =  '0'
       landiceFLG   =  '0'
       lakeFLG      =  '0'
       routeFLG     =  '0'
       catch2script =  mk_catch_j_template.format(Bin = bindir, account = self.slurm_options['account'], \
                  out_dir = out_dir, mk_catch_log = 'mk_catch_log.2', surflay = self.surf_out['surflay'],  \
                  wemin = self.surf_in['wemin'], wemout = self.surf_out['wemout'] ,  \
                  fromGCM = '0', catchFLG = catchFLG, catchcnFLG = catchcnFLG, rescale = '1', rsttime = self.common_in['yyyymmddhh'], \
                  RESTART_ID = self.surf_in.get('restart_id'), RESTART_PATH = self.surf_in.get('restart_path'), \
                  RESTART_DOMAIN = self.surf_in.get('restart_domain'), CN_VERSION = self.surf_in.get('cn_version'), \
                  saltwaterFLG = saltwaterFLG,  saltwater = self.saltwater,  \
                  openwaterFLG = openwaterFLG,  openwater = self.openwater,  \
                  seaiceFLG    = seaiceFLG,     seaice    = self.seaice,  \
                  landiceFLG   = landiceFLG,    landice   = self.landice,  \
                  lakeFLG      = lakeFLG,       lake      = self.lake,   \
                  split_saltwaterFLG = split_saltwaterFLG, \
                  routeFLG     = routeFLG,      route     = self.route, yyyymm = str(self.common_in['yyyymmddhh'])[0:6], \
                  zoom = self.surf_in['zoom'] ) 
     catch2 = open('mk_catch.j.2','wt')
     catch2.write(catch2script)
     catch2.close()
     print("step 2: sbatch -W mk_catch.j.2")
     subprocess.call(['sbatch','-W', 'mk_catch.j.2'])
     cwd = os.getcwd()
     for out_rst in glob.glob("OutData/*_rst*"):
       filename = os.path.basename(out_rst)
       print( "\nmove " + out_rst + " to " + cwd+"/"+filename)
       shutil.move(out_rst, cwd+"/"+filename)     
     os.chdir(bindir)     
