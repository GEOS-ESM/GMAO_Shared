#!/usr/bin/env python
#
import os
import subprocess

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

     self.lake = self.restarts_in.get('lake')
     self.catch = self.restarts_in.get('catch')
     self.catchcnclm40 = self.restarts_in.get('catchcnclm40')
     self.catchcnclm45 = self.restarts_in.get('catchcnclm45')
         
  def regrid(self):
     bindir = os.getcwd()
     if (self.lake):
       cmd = bindir + 'OutData/*.til Indata/*.til Indata/' + self.lake + '19 ' + str(self.surf_in['zoom'])
       outdir = self.common_out['outdir']
       if not os.path.exists(outdir) : os.makedirs(outdir)
       os.chdir(self.common_out['outdir'])
       print(cmd)
       #subprocess.run([cmd])

     mk_catch_j_1_template = """
#!/bin/csh -f
#SBATCH --account={account}
#SBATCH --ntasks={numtasks}
#SBATCH --time=1:00:00
#SBATCH --job-name=catchj
#SBATCH --output={outdir}/{mk_catch_log}
#SBATCH --qos={qos}

source {Bin}/g5_modules
set echo

#limit stacksize unlimited
unlimit

set esma_mpirun_X = ( {Bin}/esma_mpirun -np {numtasks} )
set mk_CatchRestarts_X   = ( $esma_mpirun_X {Bin}/mk_CatchRestarts )
set Scale_Catch_X   = {Bin}/Scale_Catch
set OUT_til   = OutData/*.til
set IN_til    = InData/*.til
set catchIN = InData/*catch_internal_rst*
set params = ( $OUT_til $IN_til $catchIN {surflay} )
$mk_CatchRestarts_X $params

if ({rescale}) then
  set catch_regrid = OutData/$catchIN:t
  set catch_scaled = $catch_regrid.scaled
  set params = ( $catchIN $catch_regrid $catch_scaled {surflay} )
  set params = ( $params {weminIN} {weminOUT} )
  $Scale_Catch_X $params

  mv $catch_regrid $catch_regrid.1
  mv $catch_scaled $catch_regrid
endif   
exit
"""
     mk_catch_j_1 = """
#!/bin/csh -f
#SBATCH --account={account}
#SBATCH --ntasks={numtasks}
#SBATCH --time=1:00:00
#SBATCH --job-name=catchj
#SBATCH --output={outdir}/{mk_catch_log}
#SBATCH --qos={qos}

source {Bin}/g5_modules
set echo
set esma_mpirun_X = ( {Bin}/esma_mpirun -np {numtasks} )
"""
     catch1script =  mk_catch_j_1_template.format(Bin = bindir, account = self.slurm_options['account'], \
                  outdir = outdir, mk_catch_log = 'mk_catch_log.1', surflay = self.surf_in['surflay'],  \
                  weminIN = self.surf_in['wemin'], weminOUT = self.surf_out['wemout'] , qos = self.slurm_options['qos'], \
                  numtasks = '84', rescale = '1' ) 
     catch1 = open('mk_catch.j.1','wt')
     catch1.write(catch1script)
     catch1.close()
