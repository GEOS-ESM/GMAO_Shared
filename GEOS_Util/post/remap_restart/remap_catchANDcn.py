#!/usr/bin/env python3
#
import os
import sys
import subprocess
import shutil
import glob
import ruamel.yaml
import shlex
from remap_base import remap_base

class catchANDcn(remap_base):
  def __init__(self, **configs):
     super().__init__(**configs)
     if self.config['input']['shared']['MERRA-2']:
       self.copy_merra2()

  def remap(self):
     if not self.config['output']['surface']['remap_catch']:
        return

     config = self.config
     rst_dir = config['input']['shared']['rst_dir']
     model = config['input']['surface']['catch_model']
     in_rstfile = ''
     yyyymmddhh_ = str(self.config['input']['shared']['yyyymmddhh'])
     time = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]
     in_rstfiles = glob.glob(rst_dir+'/*'+model+'_*'+time+'*')
     if len(in_rstfiles) == 0:
        print('\n try catchXX file without time stamp')
        in_rstfiles = glob.glob(rst_dir+'/*'+model+'_*')
     if len(in_rstfiles) == 0:
        return
     in_rstfile = in_rstfiles[0]  
             
     print("\nRemapping " + model + ".....\n")

     cwdir  = os.getcwd()
     bindir = os.path.dirname(os.path.realpath(__file__))
     in_bcsdir  = config['input']['shared']['bcs_dir']
     out_bcsdir = config['output']['shared']['bcs_dir']
     out_dir    = config['output']['shared']['out_dir']
     expid      = config['output']['shared']['expid']
     in_wemin   = config['input']['surface']['wemin']
     out_wemin  = config['output']['surface']['wemin']
     surflay    = config['output']['surface']['surflay']
     in_tilefile  = glob.glob(in_bcsdir+ '/*-Pfafstetter.til')[0]
     out_tilefile = glob.glob(out_bcsdir+ '/*-Pfafstetter.til')[0]
     account    = config['slurm']['account']
     # even the input is binary, the output si nc4
     suffix     = time+'z.nc4'

     if (expid) :
        expid = expid + '.'
     else:
        expid = ''
     suffix = '_rst.' + suffix
     out_rstfile = expid + os.path.basename(in_rstfile).split('_rst')[0].split('.')[-1]+suffix

     if not os.path.exists(out_dir) : os.makedirs(out_dir)
     print( "cd " + out_dir)
     os.chdir(out_dir)

     InData_dir = out_dir+'/InData/'
     print ("mkdir -p " + InData_dir)
     os.makedirs(InData_dir, exist_ok = True)
    
     f = os.path.basename(in_rstfile)
     dest = InData_dir+'/'+f
     # file got copy because the computing node cannot access archive
     print('\nCopy ' + in_rstfile + ' to ' +dest)
     shutil.copyfile(in_rstfile,dest)
     in_rstfile = dest

     log_name = out_dir+'/'+'mk_catchANDcn_log'
     mk_catch_j_template = """#!/bin/csh -f
#SBATCH --account={account}
#SBATCH --ntasks=56
#SBATCH --time=1:00:00
#SBATCH --job-name=mk_catchANDcn
#SBATCH --qos=debug
#SBATCH --output={log_name}
#

source {Bin}/g5_modules

limit stacksize unlimited

set esma_mpirun_X = ( {Bin}/esma_mpirun -np 56 )
set mk_catchANDcnRestarts_X   = ( {Bin}/mk_catchANDcnRestarts.x )

set params = ( -model {model}  -time {time} -in_tilefile {in_tilefile} )
set params = ( $params -out_bcs {out_bcs} -out_tilefile {out_tilefile} -out_dir {out_dir} )
set params = ( $params -surflay {surflay} -in_wemin {in_wemin} -out_wemin {out_wemin} ) 
set params = ( $params -in_rst {in_rstfile} -out_rst {out_rstfile} ) 
$esma_mpirun_X $mk_catchANDcnRestarts_X $params

"""
     catch1script =  mk_catch_j_template.format(Bin = bindir, account = account, out_bcs = out_bcsdir, \
                  model = model, out_dir = out_dir, surflay = surflay, log_name = log_name,  \
                  in_wemin   = in_wemin, out_wemin = out_wemin, out_tilefile = out_tilefile, in_tilefile = in_tilefile, \
                  in_rstfile = in_rstfile, out_rstfile = out_rstfile, time = yyyymmddhh_ )

     script_name = './mk_catchANDcn.j'

     catch_scrpt = open(script_name,'wt')
     catch_scrpt.write(catch1script)
     catch_scrpt.close()
    
     interactive = os.getenv('SLURM_JOB_ID', default = None)
     if ( interactive ) :
       print('interactive mode\n')
       ntasks = os.getenv('SLURM_NTASKS', default = None)
       if ( not ntasks):
         nnodes = int(os.getenv('SLURM_NNODES', default = '1'))
         ncpus  = int(os.getenv('SLURM_CPUS_ON_NODE', default = '28'))
         ntasks = nnodes * ncpus
       ntasks = int(ntasks)
       NPE = 56
       if (ntasks < NPE):
         print("\nYou should have at least {NPE} cores. Now you only have {ntasks} cores ".format(NPE=NPE, ntasks=ntasks))

       subprocess.call(['chmod', '755', script_name])
       print(script_name+  '  1>' + log_name  + '  2>&1')
       os.system(script_name + ' 1>' + log_name+ ' 2>&1')

     else:
       print("sbatch -W " + script_name +"\n")
       subprocess.call(['sbatch','-W', script_name])

     print( "cd " + cwdir)
     os.chdir(cwdir)

  def copy_merra2(self):
    if not self.config['input']['shared']['MERRA-2']:
      return

    expid = self.config['input']['shared']['expid']
    yyyymmddhh_ = str(self.config['input']['shared']['yyyymmddhh'])
    yyyy_ = yyyymmddhh_[0:4]
    mm_   = yyyymmddhh_[4:6]
    dd_   = yyyymmddhh_[6:8]
    hh_   = yyyymmddhh_[8:10]

    suffix = yyyymmddhh_[0:8]+'_'+ hh_ + 'z.bin'
    merra_2_rst_dir = '/archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/'+expid +'/rs/Y'+yyyy_ +'/M'+mm_+'/'
    rst_dir = self.config['input']['shared']['rst_dir'] + '/'
    os.makedirs(rst_dir, exist_ok = True)
    print(' Copy MERRA-2 catchment Restart \n from \n    ' + merra_2_rst_dir + '\n to\n    '+ rst_dir +'\n')

    f =  merra_2_rst_dir +  expid+'.catch_internal_rst.'    + suffix

    fname = os.path.basename(f)
    dest = rst_dir + '/'+fname
    print("Copy file "+f +" to " + rst_dir)
    shutil.copy(f, dest)

if __name__ == '__main__' :
   catch = catchANDcn(params_file='remap_params.yaml')
   catch.remap()
