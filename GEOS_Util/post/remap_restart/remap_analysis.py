#!/usr/bin/env python3
#
# source install/bin/g5_modules
#
# Newer GEOS code should load a module with GEOSpyD Python3 if not run:
#   module load python/GEOSpyD/Min4.10.3_py3.9
#

import os
from datetime import datetime, timedelta
import subprocess
import shlex
import shutil
import glob
import fileinput
import ruamel.yaml
from remap_base import remap_base

class analysis(remap_base):
  def __init__(self, **configs):
     super().__init__(**configs)
     self.copy_merra2()

  def remap(self):
     config = self.config
     bkg = config['output']['analysis']['bkg']
     if ( not bkg ): return

     analysis_in = self.find_analysis()
     if len(analysis_in) ==0 :
       print("\n There are no analysis files. \n")
       return

     print("\n Remapping or copying analysis files...\n")

     cwdir  = os.getcwd()
     bindir = os.path.dirname(os.path.realpath(__file__))
     in_bcsdir  = config['input']['shared']['bcs_dir']
     out_bcsdir = config['output']['shared']['bcs_dir']
     out_dir    = config['output']['shared']['out_dir']
     if not os.path.exists(out_dir) : os.makedirs(out_dir)
     print( "cd " + out_dir)
     os.chdir(out_dir)

     tmpdir = out_dir+'/ana_data/'
     if os.path.exists(tmpdir) : subprocess.call(['rm', '-rf',tmpdir])
     print ("mkdir " + tmpdir)
     os.makedirs(tmpdir)

     print( "cd " + tmpdir)
     os.chdir(tmpdir)

     yyyymmddhh_ = str(config['input']['shared']['yyyymmddhh'])
     yyyy_ = yyyymmddhh_[0:4]
     mm_   = yyyymmddhh_[4:6]
     dd_   = yyyymmddhh_[6:8]
     hh_   = yyyymmddhh_[8:10]
     rst_time = datetime(year=int(yyyy_), month=int(mm_), day=int(dd_), hour = int(hh_))
     expid_in  = config['input']['shared']['expid']
     expid_out = config['output']['shared']['expid']
     if (expid_out) :
        expid_out = expid_out + '.'
     else:
        expid_out = ''

     aqua   = config['output']['analysis']['aqua']
     local_fs=[]
     for f in analysis_in:
       print(f)
       fname    = os.path.basename(f)
       out_name = fname.replace(expid_in + '.', expid_out)
       f_tmp = tmpdir+'/'+out_name
       local_fs.append(f_tmp) 
       shutil.copy(f,f_tmp)
       if out_name.find('satbias') != -1 :
         if (aqua):
            f_ = open(f_tmp, 'w')
            for line in fileinput.input(f):
               f_.write(line.replace('airs281SUBSET_aqua', 'airs281_aqua      '))
            f_.close() 

     nlevel    = config['output']['air']['nlevel']
     agrid_out = config['output']['shared']['agrid']
     flags = "-g5 -res " + self.get_grid_kind(agrid_out.upper()) + " -nlevs " + str(nlevel)
     bkg_files = glob.glob(tmpdir+'/*.bkg??_eta_rst*')
     for f in bkg_files:
        f_orig = f + ".orig"
        shutil.move(f,f_orig)
        cmd = bindir + '/dyn2dyn.x ' + flags + ' -o ' + f + ' ' + f_orig
        print(cmd)
        subprocess.call(shlex.split(cmd))

     for f in local_fs:
       fname = os.path.basename(f)
       shutil.move(f, out_dir+'/'+fname)
    # write lcv 
     lcv = config['output']['analysis']['lcv']
     if lcv :
       ymd_ = yyyymmddhh_[0:8]
       hh_  = yyyymmddhh_[8:10]
       hms_ = hh_+'0000'
       rstlcvOut = out_dir+'/'+expid_out+'rst.lcv.'+ymd_+'_'+hh_+'z.bin'
       cmd = bindir+'/mkdrstdate.x ' + ymd_ + ' ' + hms_ +' ' + rstlcvOut
       print(cmd)
       subprocess.call(shlex.split(cmd))
     print( "cd " + cwdir)
     os.chdir(cwdir)

     self.remove_merra2()

  def get_grid_kind(this, grid):
     hgrd = {}
     hgrd['C12']   = 'a'
     hgrd['C24']   = 'a'
     hgrd['C48']   = 'b'
     hgrd['C90']   = 'c'
     hgrd['C180']  = 'd'
     hgrd['C360']  = 'd'
     hgrd['C500']  = 'd'
     hgrd['C720']  = 'e'
     hgrd['C1000'] = 'e'
     hgrd['C1440'] = 'e'
     hgrd['C2000'] = 'e'
     hgrd['C2880'] = 'e'
     hgrd['C5760'] = 'e'
     return hgrd[grid]

  def find_analysis(self):
     analysis_in = []
     rst_dir = self.config['input']['shared']['rst_dir']
     bkgs = glob.glob(rst_dir + '/*_eta_rst*')
     sfcs = glob.glob(rst_dir + '/*_sfc_rst*')
     anasat = glob.glob(rst_dir + '/*ana_satb*')
     traks= glob.glob(rst_dir + '/*.trak.GDA.rst*')
     analysis_in = bkgs + sfcs + traks + anasat
     return list(dict.fromkeys(analysis_in))

  def copy_merra2(self):
    if not self.config['input']['shared']['MERRA-2']:
      return
    bkg = self.config['output']['analysis']['bkg']
    if ( not bkg ): return

    expid = self.config['input']['shared']['expid']
    yyyymmddhh_ = str(self.config['input']['shared']['yyyymmddhh'])
    yyyy_ = yyyymmddhh_[0:4]
    mm_   = yyyymmddhh_[4:6]
    dd_   = yyyymmddhh_[6:8]
    hh_   = yyyymmddhh_[8:10]

    merra_2_rst_dir = '/archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/'+expid +'/rs/Y'+yyyy_ +'/M'+mm_+'/'
    rst_dir = self.config['input']['shared']['rst_dir'] + '/'
    os.makedirs(rst_dir, exist_ok = True)
    print(' Copy MERRA-2 analysis files \n from \n    ' + merra_2_rst_dir + '\n to\n    '+ rst_dir +'\n')

    rst_time = datetime(year=int(yyyy_), month=int(mm_), day=int(dd_), hour = int(hh_))

    anafiles=[]
    for h in [3,4,5,6,7,8,9]:
       delt = timedelta(hours = h-3)
       new_time = rst_time + delt
       yyyy = "Y"+str(new_time.year)
       mm   = 'M%02d'%new_time.month
       ymd  = '%04d%02d%02d'%(new_time.year,new_time.month, new_time.day)
       hh   = '%02d'%h
       newhh= '%02d'%new_time.hour
       m2_rst_dir = merra_2_rst_dir.replace('Y'+yyyy_,yyyy).replace('M'+mm_,mm)
       # bkg files
       for ftype in ['sfc', 'eta']:
          fname = expid + '.bkg'+hh+'_'+ftype+'_rst.'+ymd+'_'+newhh+'z.nc4'
          f = m2_rst_dir+'/'+fname
          if(os.path.isfile(f)):
             anafiles.append(f)
          else:
             print('Warning: Cannot find '+f)

       # gaas_bkg_sfc files
       if (h==6 or h==9):
          fname = expid+'.gaas_bkg_sfc_rst.'+ymd+'_'+newhh+'z.nc4'
          f = m2_rst_dir+'/'+fname
          if (os.path.isfile(f)):
            anafiles.append(f)
          else:
            print('Warning: Cannot find '+f)
    # trak.GDA.rst file
    delt = timedelta(hours = 3)
    new_time = rst_time - delt
    yyyy = "Y"+str(new_time.year)
    mm   = 'M%02d'%new_time.month
    ymdh = '%04d%02d%02d%02d'%(new_time.year, new_time.month, new_time.day, new_time.hour)
    m2_rst_dir = merra_2_rst_dir.replace('Y'+yyyy_,yyyy).replace('M'+mm_,mm)
    fname = expid+'.trak.GDA.rst.'+ymdh+'z.txt'
    f = m2_rst_dir+'/'+fname
    if (os.path.isfile(f)): anafiles.append(f)

    for f in anafiles:
      fname    = os.path.basename(f)
      f_tmp = rst_dir+'/'+fname
      print("Copy file "+f +" to " + rst_dir)
      shutil.copy(f,f_tmp)

if __name__ == '__main__' :
   ana = analysis(params_file='remap_params.yaml')
   ana.remap()
