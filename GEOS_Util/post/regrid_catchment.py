#!/usr/bin/env python3
#
import os
import subprocess
import shutil
import glob
import yaml
import shlex

class catchment(object):
  def __init__(self, params_file):
     stream = open(params_file, 'r')
     self.config = yaml.full_load(stream)
     stream.close()

  def regrid(self):
     print("\nRegridding catchment.....\n")
     config = self.config
     bindir  = os.getcwd()
     in_bcsdir  = config['input']['shared']['bcs_dir']
     out_bcsdir = config['output']['shared']['bcs_dir']
     out_dir    = config['output']['shared']['out_dir']
     if not os.path.exists(out_dir) : os.makedirs(out_dir)
     print( "cd " + out_dir)
     os.chdir(out_dir)

     InData_dir = out_dir+'/InData/'
     if os.path.exists(InData_dir) : subprocess.call(['rm', '-rf',InData_dir])
     print ("mkdir " + InData_dir)
     os.makedirs(InData_dir)

     OutData_dir = out_dir+'/OutData/'
     if os.path.exists(OutData_dir) : subprocess.call(['rm', '-rf',OutData_dir])
     print ("mkdir " + OutData_dir)
     os.makedirs(OutData_dir)

     rst_dir = config['input']['surface']['rst_dir']
     if (not rst_dir): rst_dir = config['input']['shared']['rst_dir']
     restarts_in = glob.glob(rst_dir +'surface/*') 
     yyyymmddhh_ = str(config['input']['shared']['yyyymmddhh'])
     suffix = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]+'z.nc4'
     if (not restarts_in) :
         # from merra-2
         restarts_in = glob.glob(rst_dir +'/*catch*')
         if restarts_in[0].find('z.bin') != -1 : suffix = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]+'z.bin'

     def find_rst(rsts,rst_name):
        for rst in rsts:
           f = os.path.basename(rst)
           if (f.find(rst_name) != -1):
              return rst
        return None  
     catch = find_rst(restarts_in,'catch_')
     if (not catch) :
        print(" No catchment restart")
        return

     catch_f = os.path.basename(catch)
     dest    = InData_dir+'/'+ catch_f
     if os.path.exists(dest) : shutil.remove(dest)
     print('\nCopy ' + catch + ' to ' +dest)
     shutil.copy(catch , dest)

     in_tile_file  = glob.glob(in_bcsdir+ '/*-Pfafstetter.til')[0]
     out_tile_file = glob.glob(out_bcsdir+ '/*-Pfafstetter.til')[0]

     in_til = InData_dir+'/' + in_tile_file.split('/')[-1]
     out_til = OutData_dir+'/'+ out_tile_file.split('/')[-1]

     if os.path.exists(in_til)  : shutil.remove(in_til)
     if os.path.exists(out_til) : shutil.remove(out_til)
     print('\n Copy ' + in_tile_file + ' to ' + in_til)
     shutil.copy(in_tile_file, in_til)
     print('\n Copy ' + out_tile_file + ' to ' + out_til)
     shutil.copy(out_tile_file, out_til)

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

set esma_mpirun_X = ( {Bin}/esma_mpirun -np 84 )
set mk_CatchRestarts_X   = ( $esma_mpirun_X {Bin}/mk_CatchRestarts )
set Scale_Catch_X   = {Bin}/Scale_Catch

set catchIN = InData/{catch_f}
set params = ( {out_til}  {in_til} $catchIN {surflay} )
$mk_CatchRestarts_X $params

if ({rescale}) then
   
   set catch_regrid = OutData/{catch_f}
   set catch_scaled = $catch_regrid.scaled
   set params = ( $catchIN $catch_regrid $catch_scaled {surflay} )
   set params = ( $params {wemin} {wemout} )
   $Scale_Catch_X $params

   rm $catch_regrid
   #mv $catch_regrid $catch_regrid.1
   mv $catch_scaled $catch_regrid
endif
"""
     # step 1
     expid = config['output']['shared']['expid']
     wemin  = config['input']['surface']['wemin']
     wemout = config['output']['surface']['wemout']
     surflay= config['output']['surface']['surflay']
     account = config['slurm']['account']

     catch1script =  mk_catch_j_template.format(Bin = bindir, account = account, \
                  out_dir = out_dir, mk_catch_log = 'mk_catch_log.1', surflay = surflay,  \
                  wemin = wemin, wemout = wemout, out_til = out_til, in_til = in_til, catch_f = catch_f, rescale = '0')
     catch1 = open('mk_catch.j.1','wt')
     catch1.write(catch1script)
     catch1.close()
     print("step 1: sbatch -W mk_catch.j.1")
     subprocess.call(['sbatch','-W', 'mk_catch.j.1'])

     # step 2
     if os.path.exists('InData.step1') : subprocess.call(['rm','-rf', 'InData.step1'])
     print("\n Move Indata to InData.step1")
     shutil.move('InData', 'InData.step1')
     os.makedirs('InData')
     for catchfile in glob.glob("OutData/*catch*"):
       print('\n Move ' + catchfile + ' to InData/')
       shutil.move(catchfile,"InData/")
     print('\n Link ' + out_til + ' to ' + in_til)
     os.symlink(out_til, in_til)

     dirname = os.path.dirname(out_til)
     clsm = dirname+'/clsm'
     print('\n Link ' + clsm + ' to ' + 'OutData/clsm')
     os.symlink(clsm, 'OutData/clsm')

     catch2script = mk_catch_j_template.format(Bin = bindir, account = account, \
                  out_dir = out_dir, mk_catch_log = 'mk_catch_log.1', surflay = surflay,  \
                  wemin = wemin, wemout = wemout, out_til = out_til, in_til = in_til, catch_f = catch_f, rescale = '1')

     catch2 = open('mk_catch.j.2','wt')
     catch2.write(catch2script)
     catch2.close()
     print("step 2: sbatch -W mk_catch.j.2")
     subprocess.call(['sbatch','-W', 'mk_catch.j.2'])
#
# post process
#
     if (expid) :
        expid = expid + '.'
     else:
        expid = ''
     suffix = '_rst.' + suffix
     for out_rst in glob.glob("OutData/*_rst*"):
       filename = expid + os.path.basename(out_rst).split('_rst')[0].split('.')[-1]+suffix
       print('\n Move ' + out_rst + ' to ' + out_dir+"/"+filename)
       shutil.move(out_rst, out_dir+"/"+filename)
     os.chdir(bindir)

if __name__ == '__main__' :
   catch = catchment('regrid_params.yaml')
   catch.regrid()

     
