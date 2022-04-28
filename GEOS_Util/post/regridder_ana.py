#!/usr/bin/env python3
#
import os
import subprocess
import shlex
import shutil
import glob
import fileinput
from regridder_base import *

class analysis(regridder):
  def __init__(self, config):
     # v3.0
     #super().__init__(config)
     # v2.7
     super(analysis, self).__init__(config)
     self.ana_in = config['input']['parameters']['ANALYSIS']
     self.bkg    = self.ana_in.get('bkg')
     if not self.bkg:
       return
     #ana bkg_eta grid mapping from input atmosphere grid
     self.hgrd = {}
     self.hgrd['C12']   = 'a'
     self.hgrd['C24']   = 'a'
     self.hgrd['C48']   = 'b'
     self.hgrd['C90']   = 'c'
     self.hgrd['C180']  = 'd'
     self.hgrd['C360']  = 'd'
     self.hgrd['C500']  = 'd'
     self.hgrd['C720']  = 'e'
     self.hgrd['C1000'] = 'e'
     self.hgrd['C1440'] = 'e'
     self.hgrd['C2000'] = 'e'
     self.hgrd['C2880'] = 'e'
     self.hgrd['C5760'] = 'e'

     rst_time = datetime(year=int(self.yyyy), month=int(self.mm), day=int(self.dd), hour = int(self.hh))
     rst_dir_orig = self.common_in['rst_dir']
     expid = self.common_in['expid']
     anafiles = self.restarts_in['ANALYSIS']
     if len(anafiles) == 0 and self.common_in.get('MERRA-2'): 
       self.anafiles=[]
       for h in [3,4,5,6,7,8,9]:
          delt = timedelta(hours = h-3)
          new_time = rst_time + delt
          yyyy = "Y"+str(new_time.year)
          mm   = 'M%02d'%new_time.month
          ymd  = '%04d%02d%02d'%(new_time.year,new_time.month, new_time.day)
          hh   = '%02d'%h
          newhh= '%02d'%new_time.hour
          rst_dir = rst_dir_orig.replace('Y'+self.yyyy,yyyy).replace('M'+self.mm,mm)
          # bkg files
          for ftype in ['sfc', 'eta']:
             fname = expid+'.bkg'+hh+'_'+ftype+'_rst.'+ymd+'_'+newhh+'z.nc4'
             f = rst_dir+'/'+fname
             if(os.path.isfile(f)):
               self.anafiles.append(f)
             else:
               print('Warning: Cannot find '+f)

          # cbkg file
          fname = expid + '.cbkg' + hh + '_eta_rst.' + ymd + '_' + newhh + 'z.nc4'
          f = rst_dir+'/'+fname
          if(os.path.isfile(f)):
            self.anafiles.append(f)
          else:
            print('Warning: Cannot find '+f)
          # gaas_bkg_sfc files
          if (h==6 or h==9):
             fname = expid+'.gaas_bkg_sfc_rst.'+ymd+'_'+newhh+'z.nc4'
             f = rst_dir+'/'+fname
             if (os.path.isfile(f)):
               self.anafiles.append(f)
             else:
               print('Warning: Cannot find '+f)
       # satbang and satbias
       ymd  = '%04d%02d%02d'%(rst_time.year,rst_time.month, rst_time.day)
       hr ='%02d'%rst_time.hour
       for ftype in ["ana_satbang_rst", "ana_satbias_rst", "ana_satbiaspc_rst"]:
          fname = expid+'.'+ftype+'.'+ymd+'_'+hr+'z.txt'
          f = rst_dir_orig+'/'+fname
          if(os.path.isfile(f)):
            self.anafiles.append(f)
          else:
            print('Warning: Cannot find '+f)

       # trak.GDA.rst file
       delt = timedelta(hours = 3)
       new_time = rst_time - delt
       yyyy = "Y"+str(new_time.year)
       mm   = 'M%02d'%new_time.month
       ymdh = '%04d%02d%02d%02d'%(new_time.year, new_time.month, new_time.day, new_time.hour)
       rst_dir = rst_dir_orig.replace('Y'+self.yyyy,yyyy).replace('M'+self.mm,mm)
       fname = expid+'.trak.GDA.rst.'+ymdh+'z.txt'
       f = rst_dir+'/'+fname
       if (os.path.isfile(f)): self.anafiles.append(f)
     else:
       self.anafiles = anafiles

     agrid_in  = self.common_in['agrid']
     agrid_out = self.common_out['agrid']

     self.ana_regrid = True
     if (self.hgrd[agrid_in] == self.hgrd[agrid_out]):
       self.ana_regrid = False

  def regrid(self):
    if not self.bkg:
      return
    expid_in = self.common_in.get('expid')+'.'
    expid_out = self.common_out.get('expid')
    if not expid_out:
      expid = ''
    else:
      expid_out = expid_out+'.'
    bindir = os.getcwd()
    out_dir = self.common_out.get('out_dir')

    if not os.path.exists(out_dir) : os.makedirs(out_dir)
    print( "cd " + self.common_out['out_dir'])
    os.chdir(self.common_out['out_dir'])
    tmp_dir = out_dir+'/ana_data/'
    if os.path.exists(tmp_dir) : subprocess.call(['rm', '-rf', tmp_dir])
    print ("mkdir " + tmp_dir)
    os.makedirs(tmp_dir)
    print( "cd " + tmp_dir)
    os.chdir(tmp_dir)

    ogrid = self.common_out['ogrid']
    tagout = self.common_out['tag']
    bctag  = self.get_bcTag(tagout, ogrid)
    tagrank = self.tagsRank[bctag]
    local_fs=[]
    for f in self.anafiles:
      fname    = os.path.basename(f)
      out_name = fname.replace(expid_in, expid_out)
      f_tmp = tmp_dir+'/'+out_name
      local_fs.append(f_tmp) 
      shutil.copy(f,f_tmp)
      if (out_name.find('satbias') != -1):
         if (tagrank >= self.tagsRank["Ganymed-4_0_Reynolds"]):
            f_orig = f_tmp+'.orig'
            shutil.move(f_tmp, f_orig)
            f_ = open(f_tmp, 'w')
            for line in fileinput.input(f_orig):
               f_.write(line.replace('airs281SUBSET_aqua', 'airs281_aqua      '))
            f_.close() 

    if self.ana_regrid:
       agrid_out = self.common_out['agrid']
       nlevel = self.upper_out['nlevel']
       flags = "-g5 -res " + self.hgrd[agrid_out] + " -nlevs " + str(nlevel)
       bkg_files = glob.glob(tmp_dir+'/*.bkg??_eta_rst*')
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
    lcv = self.ana_in.get('lcv')
    if lcv :
       ymd = self.ymd
       hh  = self.hh
       hms = hh+'0000'
       rstlcvOut = out_dir+'/'+expid_out+'rst.lcv.'+ymd+'_'+hh+'z.bin'
       cmd = bindir+'/mkdrstdate.x ' + ymd + ' ' + hms +' ' + rstlcvOut
       print(cmd)
       subprocess.call(shlex.split(cmd))
    os.chdir(bindir)
