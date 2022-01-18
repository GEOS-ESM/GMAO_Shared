#!/usr/bin/env python3
#
import os
import subprocess
import shutil
import glob
from regrider_base import *

class upperair(regrider):
  def __init__(self, config):
     # v3.0
     #super().__init__(config)
     # v2.7
     super(upperair, self).__init__(config)
     self.upper_out = config['output']['parameters']['UPPERAIR']
     self.restarts_in = self.restarts_in['UPPERAIR']
      
     # verify agrid
     agrid = self.common_in['agrid']
     ogrid = self.common_out['ogrid']
     if not self.common_in.get('MERRA-2'):
       fvcore = '' 
       for f in self.restarts_in:
           if 'fvcore' in f:
             fvcore = f
             break
       path_ = os.path.dirname(os.path.realpath(__file__)) 
       cmd = path_ + '/fvrst.x -h ' + fvcore
       print(cmd +'\n')
       p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
       (output, err) = p.communicate()
       p_status = p.wait()
       ss = output.decode().split()
       if (agrid):
         if agrid[0].upper() == "C":
            n=int(agrid[1:])
            o=int(ss[0])
            assert n==o, "input agrid is not consistent with fvcore restart"
       else:
         self.common_in['agrid'] = "C"+ss[0]

     tagout = self.common_out['tag']
     bctag  = self.get_bcTag(tagout, ogrid) 
     tagrank = self.tagsRank[bctag]

     self.upper_out['drymassFLG'] = 0
     if tagrank >=12 :
       self.upper_out['drymassFLG'] = 1

     agrid = self.common_out['agrid']
     if agrid[0].upper() == 'C':
       self.upper_out['imout']= int(agrid[1:])
     else:
       print("Only support cs grid so far")
       exit()

     imout = self.upper_out['imout']
     if (imout <=90):
       self.upper_out['NPE'] = 12;   self.upper_out['nwrit'] = 1
     elif (imout==180):
       self.upper_out['NPE'] = 24;   self.upper_out['nwrit'] = 1
     elif (imout<=500):
       self.upper_out['NPE'] = 96;   self.upper_out['nwrit'] = 1
     elif (imout==720):
       self.upper_out['NPE'] = 192;  self.upper_out['nwrit'] = 2
     elif (imout==1000):
       self.upper_out['NPE'] = 384;  self.upper_out['nwrit'] = 2
     elif (imout==1440):
       self.upper_out['NPE'] = 576;  self.upper_out['nwrit'] = 2
     elif (imout==2000):
       self.upper_out['NPE'] = 768;  self.upper_out['nwrit'] = 2
     elif (imout>=2880):
       self.upper_out['NPE'] = 5400; self.upper_out['nwrit'] = 6

     self.upper_out['QOS'] = "#"
     if self.upper_out['NPE'] <= 532: self.upper_out['QOS'] = "#SBATCH --qos=debug"

  def regrid(self):
     print( "\nRegridding upper air......\n")
     bindir = os.getcwd()
     out_dir = self.common_out['out_dir']
     if not os.path.exists(out_dir) : os.makedirs(out_dir)
     print( "cd " + self.common_out['out_dir'])
     os.chdir(self.common_out['out_dir'])
     tmpdir = out_dir+'/upper_data/'
     if os.path.exists(tmpdir) : subprocess.call('rm -rf '+ tmpdir, shell = True)
     print ("mkdir " + tmpdir)
     os.makedirs(tmpdir)

     print( "cd " + tmpdir)
     os.chdir(tmpdir)
     print('\nUpper air restart files should end with "_rst" \n')
     for rst in self.restarts_in :
       f = os.path.basename(rst)
       f = f.replace('_rst','_restart_in')
       cmd = '/bin/cp  '+rst+' '+f
       print('\n'+cmd)
       subprocess.call(cmd, shell = True)
 
     # link topo file
     topoin = glob.glob(self.in_bcsdir+'/topo_DYN_ave*')[0]
     cmd = '/bin/cp ' + topoin + ' .'
     print('\n'+cmd)
     subprocess.call(cmd, shell = True)

     topoout = glob.glob(self.out_bcsdir+'/topo_DYN_ave*')[0]
     cmd = '/bin/cp ' + topoout + ' .'
     print('\n'+cmd)
     subprocess.call(cmd, shell = True)
     fname = os.path.basename(topoout)
     cmd = '/bin/ln -s ' + fname + ' topo_dynave.data'
     print('\n'+cmd)
     subprocess.call(cmd, shell = True)

     regrid_template="""#!/bin/csh -xf
#!/bin/csh -xf
#SBATCH --account={account}
#SBATCH --time=1:00:00
#SBATCH --ntasks={NPE}
#SBATCH --job-name=regrid_upper
#SBATCH --output={out_dir}/{out_log}
{QOS}

unlimit

cd {out_dir}/upper_data
source {Bin}/g5_modules
/bin/touch input.nml

# The MERRA fvcore_internal_restarts don't include W or DZ, but we can add them by setting 
# HYDROSTATIC = 0 which means HYDROSTATIC = FALSE
set HYDROSTATIC = 0

if ($?I_MPI_ROOT) then
  # intel scaling suggestions
  #--------------------------
  
  setenv I_MPI_DAPL_UD on

  setenv DAPL_UCM_CQ_SIZE 4096
  setenv DAPL_UCM_QP_SIZE 4096

  setenv I_MPI_DAPL_UD_SEND_BUFFER_NUM 4096
  setenv I_MPI_DAPL_UD_RECV_BUFFER_NUM 4096
  setenv I_MPI_DAPL_UD_ACK_SEND_POOL_SIZE 4096
  setenv I_MPI_DAPL_UD_ACK_RECV_POOL_SIZE 4096
  setenv I_MPI_DAPL_UD_RNDV_EP_NUM 2
  setenv I_MPI_DAPL_UD_REQ_EVD_SIZE 2000

  setenv DAPL_UCM_REP_TIME 2000
  setenv DAPL_UCM_RTU_TIME 2000
  setenv DAPL_UCM_RETRY 7
  setenv DAPL_ACK_RETRY 7
  setenv DAPL_ACK_TIMER 20
  setenv DAPL_UCM_RETRY 10
  setenv DAPL_ACK_RETRY 10

else if ($?MVAPICH2) then
  setenv MV2_ENABLE_AFFINITY 0

endif
set infiles = ()
set outfils = ()
foreach infile ( *_restart_in )
   if ( $infile == fvcore_internal_restart_in ) continue
   if ( $infile == moist_internal_restart_in  ) continue

   set infiles = ( $infiles $infile )
   set outfil = `echo $infile | sed "s/restart_in/rst_out/"`
   set outfils = ($outfils $outfil)
end

if ( $#infiles ) then
    set ioflag = "-input_files $infiles -output_files $outfils"
else
    set ioflag = ""
endif

set drymassFLG = {drymassFLG}
if ($drymassFLG) then
    set dmflag = ""
else
    set dmflag = "-scalers F"
endif

set interp_restartsX = {Bin}/interp_restarts.x

{Bin}/esma_mpirun -np {NPE} $interp_restartsX -im {imout} -lm {nlevel} \\
   -do_hydro $HYDROSTATIC $ioflag $dmflag -nwriter {nwrit}

"""
     regrid_upper_script = regrid_template.format(Bin=bindir, account = self.slurm_options['account'], \
             out_dir = out_dir, out_log = 'regrid_upper_log', drymassFLG = self.upper_out['drymassFLG'], \
             imout = self.upper_out['imout'], nwrit = self.upper_out['nwrit'], NPE = self.upper_out['NPE'], \
             QOS = self.upper_out['QOS'], nlevel = self.upper_out['nlevel'])
     upper = open('regrider_upper.j','wt')
     upper.write(regrid_upper_script)
     upper.close()
     print('sbatch -W regrider_upper.j\n')
     subprocess.call('sbatch -W regrider_upper.j', shell= True)
     for out_rst in glob.glob("*_rst*"):
       filename = os.path.basename(out_rst)
       print('\n Move ' + out_rst + ' to ' + out_dir+"/"+filename)
       shutil.move(out_rst, out_dir+"/"+filename)
     os.chdir(bindir)
