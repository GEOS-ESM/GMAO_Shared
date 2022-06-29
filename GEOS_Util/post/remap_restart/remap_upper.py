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
     f = os.path.basename(params_file)
     out_dir    = self.config['output']['shared']['out_dir']
     if not os.path.exists(out_dir) : os.makedirs(out_dir)
     dest = out_dir+'/'+f
     try:
       shutil.copy(params_file, dest)  
     except shutil.SameFileError:
       pass
  def remap(self):
     if not self.config['output']['air']['remap'] :
        return
     self.air_restarts =["fvcore_internal_rst"      , 
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
     restarts_in = self.find_rst()
     if len(restarts_in) == 0:
       return

     print( "\nRemapping upper air......\n")
     config = self.config
     cwdir  = os.getcwd()
     bindir  = os.path.dirname(os.path.realpath(__file__))
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

     print('\nUpper air restart file names link from "_rst" to "_restart_in" \n')

     types = 'z.bin'
     type_str = subprocess.check_output(['file','-b', restarts_in[0]])
     type_str = str(type_str)
     if type_str.find('Hierarchical') >=0:
        types = 'z.nc4'
     yyyymmddhh_ = str(config['input']['shared']['yyyymmddhh'])
     suffix = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]+ types
     
     for rst in restarts_in :
       f = os.path.basename(rst).split('_rst')[0].split('.')[-1]+'_restart_in'
       cmd = '/bin/ln -s  ' + rst + ' ' + f
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

     if (imout <90):
       NPE = 12; nwrit = 1
     elif (imout<=180):
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

     QOS = "#SBATCH --qos="+config['slurm']['qos']
     if NPE > 532: QOS = "###" + QOS
     CONSTR = "#SBATCH --constraint=" + config['slurm']['constraint']    

     log_name = out_dir+'/remap_upper_log'

     remap_template="""#!/bin/csh -xf
#!/bin/csh -xf
#SBATCH --account={account}
#SBATCH --time=1:00:00
#SBATCH --ntasks={NPE}
#SBATCH --job-name=remap_upper
#SBATCH --output={log_name}
{QOS}
{CONSTR}

unlimit

cd {out_dir}/upper_data
source {Bin}/g5_modules
/bin/touch input.nml

# The MERRA fvcore_internal_restarts don't include W or DZ, but we can add them by setting 
# HYDROSTATIC = 0 which means HYDROSTATIC = FALSE

if ($?I_MPI_ROOT) then
  # intel scaling suggestions
  #--------------------------
  setenv I_MPI_ADJUST_ALLREDUCE 12
  setenv I_MPI_ADJUST_GATHERV 3

  setenv I_MPI_SHM_HEAP_VSIZE 512
  setenv PSM2_MEMORY large
  setenv I_MPI_EXTRA_FILESYSTEM 1
  setenv I_MPI_EXTRA_FILESYSTEM_FORCE gpfs
  setenv ROMIO_FSTYPE_FORCE "gpfs:"

else if ($?MVAPICH2) then

  setenv MV2_ENABLE_AFFINITY 0
  setenv MV2_ENABLE_AFFINITY     0
  setenv SLURM_DISTRIBUTION block
  setenv MV2_MPIRUN_TIMEOUT 100
  setenv MV2_GATHERV_SSEND_THRESHOLD 256

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
    set ftype = `file -Lb --mime-type fvcore_internal_restart_in`
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
     
     remap_upper_script = remap_template.format(Bin=bindir, account = account, \
             out_dir = out_dir, log_name = log_name, drymassFLG = drymassFLG, \
             imout = imout, nwrit = nwrit, NPE = NPE, \
             QOS = QOS,CONSTR = CONSTR, nlevel = nlevel, hydrostatic = hydrostatic)

     script_name = './remap_upper.j'

     upper = open(script_name,'wt')
     upper.write(remap_upper_script)
     upper.close()

     interactive = os.getenv('SLURM_JOB_ID', default = None)

     if (interactive) :
       print('interactive mode\n')
       ntasks = os.getenv('SLURM_NTASKS', default = None)
       if ( not ntasks):
         nnodes = int(os.getenv('SLURM_NNODES', default = '1'))
         ncpus  = int(os.getenv('SLURM_CPUS_ON_NODE', default = '28'))
         ntasks = nnodes * ncpus
       ntasks = int(ntasks)
       if (ntasks < NPE ):
         print("\nYou should have at least {NPE} cores. Now you only have {ntasks} cores ".format(NPE=NPE, ntasks=ntasks))
         
       subprocess.call(['chmod', '755', script_name])
       print(script_name+  '  1>' + log_name  + '  2>&1')
       os.system(script_name + ' 1>' + log_name+ ' 2>&1')
     else : 
       print('sbatch -W '+ script_name +'\n')
       subprocess.call(['sbatch', '-W', script_name])

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

     print('\n Move remap_upper.j to ' + out_dir)
     shutil.move('remap_upper.j', out_dir+"/remap_upper.j")
     print('cd ' + cwdir)
     os.chdir(cwdir)

  def find_rst(self):
     rst_dir = self.config['input']['shared']['rst_dir']
     yyyymmddhh_ = str(self.config['input']['shared']['yyyymmddhh'])
     time = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]
     restarts_in=[]
     for f in self.air_restarts :
        files = glob.glob(rst_dir+ '/*'+f+'*'+time+'*')
        if len(files) >0:
          restarts_in.append(files[0])
     if (len(restarts_in) == 0) :
        print("\n try restart file names without time stamp\n")
        for f in self.air_restarts :
           fname = rst_dir+ '/'+f
           if os.path.exists(fname):
             restarts_in.append(fname)
        
     return restarts_in

if __name__ == '__main__' :
   air = upperair('remap_params.yaml')
   air.remap()
