#!/usr/bin/env python3
#
import os
import ruamel.yaml
import subprocess
import shlex
import shutil
import glob

class upperair(object):
  def __init__(self, params_file):
     yaml = ruamel.yaml.YAML()
     stream ='' 
     with  open(params_file, 'r') as f:
        stream = f.read()
     self.config = yaml.load(stream)

  def regrid(self):
     restarts_in = self.find_rst()
     if len(restarts_in) == 0:
       return

     print( "\nRegridding upper air......\n")
     config = self.config
     bindir  = os.getcwd()
     in_bcsdir  = config['input']['shared']['bcs_dir'] 
     out_bcsdir = config['output']['shared']['bcs_dir'] 
     out_dir    = config['output']['shared']['out_dir']

     if not os.path.exists(out_dir) : os.makedirs(out_dir)
     print( "cd " + out_dir)
     os.chdir(out_dir)

     tmpdir = out_dir+'/upper_data/'
     if os.path.exists(tmpdir) : subprocess.call(['rm', '-rf',tmpdir])
     print ("mkdir " + tmpdir)
     os.makedirs(tmpdir)

     print( "cd " + tmpdir)
     os.chdir(tmpdir)

     print('\nUpper air restart file names changed from "_rst" to "_restart_in" \n')

     types = 'z.bin'
     type_str = subprocess.check_output(['file','-b', restarts_in[0]])
     type_str = str(type_str)
     if type_str.find('Hierarchical') >=0:
        types = 'z.nc4'
     yyyymmddhh_ = str(config['input']['shared']['yyyymmddhh'])
     suffix = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]+ types
     
     for rst in restarts_in :
       fs = os.path.basename(rst).split('.')
       f = fs[0]
       if (len(fs) >=2): f = fs[1] 
       f = f.replace('_rst','_restart_in')
       cmd = '/bin/cp  ' + rst + ' ' + f
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))
 
     # link topo file
     topoin = glob.glob(in_bcsdir+'/topo_DYN_ave*')[0]
     cmd = '/bin/cp ' + topoin + ' .'
     print('\n'+cmd)
     subprocess.call(shlex.split(cmd))

     topoout = glob.glob(out_bcsdir+'/topo_DYN_ave*')[0]
     cmd = '/bin/cp ' + topoout + ' .'
     print('\n'+cmd)
     subprocess.call(shlex.split(cmd))
     fname = os.path.basename(topoout)
     cmd = '/bin/ln -s ' + fname + ' topo_dynave.data'
     print('\n'+cmd)
     subprocess.call(shlex.split(cmd))

     agrid  = config['output']['shared']['agrid']
     if agrid[0].upper() == 'C':
       imout = int(agrid[1:])
     else:
       exit("Only support cs grid so far")

     if (imout <=90):
       NPE = 12; nwrit = 1
     elif (imout==180):
       NPE = 24; nwrit = 1
     elif (imout<=500):
       NPE = 96; nwrit = 1
     elif (imout==720):
       NPE = 192; nwrit = 2
     elif (imout==1000):
       NPE = 384; nwrit = 2
     elif (imout==1440):
       NPE = 576; nwrit = 2
     elif (imout==2000):
       NPE = 768; nwrit = 2
     elif (imout>=2880):
       NPE = 5400; nwrit= 6

     QOS = "#"
     if NPE <= 532: QOS = "#SBATCH --qos=debug"

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

set interp_restartsX = {Bin}/interp_restarts.x
if ( $#infiles ) then
    set ioflag = "-input_files $infiles -output_files $outfils"
    set ftype = `file -b --mime-type fvcore_internal_restart_in`
    if ($ftype =~ *stream*) then
      set interp_restartsX = {Bin}/interp_restarts_bin.x
    endif
else
    set ioflag = ""
endif

set drymassFLG = {drymassFLG}
if ($drymassFLG) then
    set dmflag = ""
else
    set dmflag = "-scalers F"
endif

{Bin}/esma_mpirun -np {NPE} $interp_restartsX -im {imout} -lm {nlevel} \\
   -do_hydro {hydrostatic} $ioflag $dmflag -nwriter {nwrit}

"""
     account = config['slurm']['account']
     drymassFLG  = config['input']['air']['drymass']
     hydrostatic = config['input']['air']['hydrostatic']
     nlevel = config['output']['air']['nlevel']
     
     regrid_upper_script = regrid_template.format(Bin=bindir, account = account, \
             out_dir = out_dir, out_log = 'regrid_upper_log', drymassFLG = drymassFLG, \
             imout = imout, nwrit = nwrit, NPE = NPE, \
             QOS = QOS, nlevel = nlevel, hydrostatic = hydrostatic)
     upper = open('regrid_upper.j','wt')
     upper.write(regrid_upper_script)
     upper.close()
     print('sbatch -W regrid_upper.j\n')
     subprocess.call(['sbatch', '-W', 'regrid_upper.j'])

#
#    post process
#
     expid = config['output']['shared']['expid']
     if (expid) :
        expid = expid + '.'
     else:
        expid = ''
     suffix = '_rst.' + suffix
     for out_rst in glob.glob("*_rst*"):
       filename = expid + os.path.basename(out_rst).split('_rst')[0].split('.')[-1]+suffix
       print('\n Move ' + out_rst + ' to ' + out_dir+"/"+filename)
       shutil.move(out_rst, out_dir+"/"+filename)

     print('\n Move regrid_upper.j to ' + out_dir)
     shutil.move('regrid_upper.j', out_dir+"/regrid_upper.j")
     print('cd ' + bindir)
     os.chdir(bindir)

  def find_rst(self):
     air_restarts =["fvcore_internal_rst"      , 
                    "moist_internal_rst"       ,
                    "agcm_import_rst"          ,
                    "agcm_internal_rst"        ,
                    "carma_internal_rst"       ,
                    "achem_internal_rst"   ,
                    "geoschemchem_internal_rst",
                    "gmichem_internal_rst"     ,
                    "gocart_internal_rst"      ,
                    "hemco_internal_rst"       ,
                    "mam_internal_rst"         ,
                    "matrix_internal_rst"      ,
                    "pchem_internal_rst"       ,
                    "stratchem_internal_rst"   ,
                    "ss_internal_rst"          ,
                    "du_internal_rst"          ,
                    "cabr_internal_rst"        ,
                    "cabc_internal_rst"        ,
                    "caoc_internal_rst"        ,
                    "ni_internal_rst"          ,
                    "su_internal_rst"          ,
                    "tr_internal_rst"]

     rst_dir = self.config['input']['shared']['rst_dir']
     restarts_in=[]
     for f in air_restarts :
        files = glob.glob(rst_dir+ '/*'+f+'*')
        if len(files) >0:
          restarts_in.append(files[0])
     return restarts_in

if __name__ == '__main__' :
   air = upperair('regrid_params.yaml')
   air.regrid()