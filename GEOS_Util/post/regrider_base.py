#!/usr/bin/env python
#
import os
import subprocess
import shutil
import glob
import time

class regrider:
  def __init__(self, config):
     self.common_in   = config['input']['parameters']['COMMON']
     self.restarts_in = config['input']['restarts']
     self.common_out  = config['output']['parameters']['COMMON']
     self.slurm_options = config['slurm_options']

     self.init_tags()

     # get bc directory and tile file
     self.in_bcsdir  = self.get_bcdir("IN")
     self.in_til     = glob.glob(self.in_bcsdir+ '/*-Pfafstetter.til')[0] 
     self.out_bcsdir = self.get_bcdir("OUT")
     self.out_til     = glob.glob(self.out_bcsdir+ '/*-Pfafstetter.til')[0] 

     #ToDo
     # 1) get rst source

  def init_tags(self):
     # copy and paste from regrid.pl
     # minor change. Add "D" to the number for each group
     # BCS Tag: Fortuna-1_4
     F14  = ( 'F14',              'Fortuna-1_4',            'Fortuna-1_4_p1' )
     D214 = ( 'D214',              'GEOSdas-2_1_4',          'GEOSdas-2_1_4-m1', 
              'GEOSdas-2_1_4-m2', 'GEOSdas-2_1_4-m3',       'GEOSdas-2_1_4-m4' )
     D540 = ( 'D540',              'GEOSadas-5_4_0',         'GEOSadas-5_4_0_p1',
              'GEOSadas-5_4_0_p2',  'GEOSadas-5_4_0_p3',    'GEOSadas-5_4_0_p4',
              'GEOSadas-5_4_1',     'GEOSadas-5_4_1_p1',    'GEOSadas-5_4_2',
              'GEOSadas-5_4_3',     'GEOSadas-5_4_4',       'GEOSadas-5_5_0',
              'GEOSadas-5_5_1',     'GEOSadas-5_5_2',       'GEOSadas-5_5_3' )

    # BCS Tag: Fortuna-2_0
    #---------------------
     F20  = ( 'F20',                    'Fortuna-2_0')

    # BCS Tag: Fortuna-2_1
    #---------------------
     F21  = ( 'F21',                  'Fortuna-2_1',         'Fortuna-2_1_p1',
              'Fortuna-2_1_p2',       'Fortuna-2_1_p3',      'Fortuna-2_2',
              'Fortuna-2_2_p1',       'Fortuna-2_2_p2',      'Fortuna-2_3',
              'Fortuna-2_3_p1',       'Fortuna-2_4',         'Fortuna-2_4_p1',
              'Fortuna-2_4_p2',       'Fortuna-2_5',         'Fortuna-2_5_BETA0',
              'Fortuna-2_5_p1',       'Fortuna-2_5_p2',      'Fortuna-2_5_p3',
              'Fortuna-2_5_p4',       'Fortuna-2_5_p5',      'Fortuna-2_5_p6',
              'Fortuna-2_5_pp2' )
     D561 = ( 'D561',                  'GEOSadas-5_6_1',      'GEOSadas-5_6_1_p1',
              'GEOSadas-5_6_1_p2',    'GEOSadas-5_6_1_p3',   'GEOSadas-5_6_1_p4',
              'GEOSadas-5_6_2',       'GEOSadas-5_6_2_p1',   'GEOSadas-5_6_2_p2',
              'GEOSadas-5_6_2_p3',    'GEOSadas-5_6_2_p4',   'GEOSadas-5_6_2_p5',
              'GEOSadas-5_6_2_p6',    'GEOSadas-5_7_1',      'GEOSadas-5_7_1_p1',
              'GEOSadas-5_7_1_p2',    'GEOSadas-5_7_2',      'GEOSadas-5_7_2_p1',
              'GEOSadas-5_7_2_p2',    'GEOSadas-5_7_2_p2_m1','GEOSadas-5_7_2_p3',
              'GEOSadas-5_7_2_p3_m1', 'GEOSadas-5_7_2_p3_m2','GEOSadas-5_7_2_p4',
              'GEOSadas-5_7_2_p5',    'GEOSadas-5_7_2_p5_m1','GEOSadas-5_7_3',
              'GEOSadas-5_7_3_p1',    'GEOSadas-5_7_3_p2',   'GEOSadas-5_7_3_p2' )

   # BCS Tag: Ganymed-1_0
    #---------------------
     G10 =  ( 'G10',                  'Ganymed-1_0',          'Ganymed-1_0_BETA',
              'Ganymed-1_0_BETA1',    'Ganymed-1_0_BETA2',    'Ganymed-1_0_BETA3',
              'Ganymed-1_0_BETA4' )

     D580 = ( 'D580',                  'GEOSadas-5_8_0',       'GEOSadas-5_9_0',
              'GEOSadas-5_9_1' )

    # BCS Tags: Ganymed-1_0_M and Ganymed-1_0_D
    #------------------------------------------
     G10p = ( 'G10p',                 'Ganymed-1_0_p1',       'Ganymed-1_0_p2',
              'Ganymed-1_0_p3',       'Ganymed-1_0_p4',       'Ganymed-1_0_p5',
              'Ganymed-1_0_p6' )

     D591p= ( 'D591p',                 'GEOSadas-5_9_1_p1',    'GEOSadas-5_9_1_p2',
              'GEOSadas-5_9_1_p3',    'GEOSadas-5_9_1_p4',    'GEOSadas-5_9_1_p5',
              'GEOSadas-5_9_1_p6',    'GEOSadas-5_9_1_p7',    'GEOSadas-5_9_1_p8',
              'GEOSadas-5_9_1_p9' )

    # BCS Tags: Ganymed-1_0_M and Ganymed-1_0_D w/ new landice rst
    #------------------------------------------------------------------------
     G20  = ( 'G20',                  'Ganymed-2_0',          'Ganymed-2_1',
              'Ganymed-2_1_p1',       'Ganymed-2_1_p2',       'Ganymed-2_1_p3',
              'Ganymed-2_1_p4',       'Ganymed-2_1_p5',       'Ganymed-2_1_p6' )
     D5A0 = ( 'D5A0',                  'GEOSadas-5_10_0',      'GEOSadas-5_10_0_p1' )


    # BCS Tags: Ganymed-1_0_Reynolds and Ganymed-1_0_Ostia
    #-----------------------------------------------------
     G30  = ( 'G30',                  'Ganymed-3_0',         'Ganymed-3_0_p1' )
     D5B0 = ( 'D5B0',                  'GEOSadas-5_10_0_p2',  'GEOSadas-5_11_0' )

    # BCS Tags: Ganymed-4_0_Reynolds, Ganymed-4_0_MERRA-2, and Ganymed-4_0_Ostia
    #---------------------------------------------------------------------------
     G40  = ( 'G40',                  'Ganymed-4_0',         'Ganymed-4_0_p1',
              'Ganymed-4_1',          'Heracles-1_0',        'Heracles-1_1',
              'Heracles-2_0',         'Heracles-2_1',        'Heracles-3_0',
              'Heracles-4_0',         'Heracles-5_4_p3' )
     D512 = ( 'D512',                  'GEOSadas-5_12_2',     'GEOSadas-5_12_4',
              'GEOSadas-5_12_4_p1',   'GEOSadas-5_12_4_p2',  'GEOSadas-5_12_4_p3',
              'GEOSadas-5_12_5',      'GEOSadas-5_13_0_p1',  'GEOSadas-5_13_0_p2',
              'GEOSadas-5_13_1',      'GEOSadas-5_16_5' )

    # BCS Tags: Icarus (New Land Parameters, New Topography)
    #---------------------------------------------------------------------------
     ICA  = ( 'ICA',                  'Icarus',              'Jason' )
     D517 = ( 'D517', 'GEOSadas-5_17_0',      'GEOSadas-5_17_1',     'GEOSadas-5_18_0',
              'GEOSadas-5_18_1',      'GEOSadas-5_18_2',     'GEOSadas-5_18_3',
              'GEOSadas-5_18_3_p1',   'GEOSadas-5_19_0',     'GEOSadas-5_20_0',
              'GEOSadas-5_20_0_p1',   'GEOSadas-5_20_0_p2',  'GEOSadas-5_21_0',
              'GEOSadas-5_21_2',      'GEOSadas-5_21_3_p1',  'GEOSadas-5_22_0',
              'GEOSadas-5_22_0_p1',   'GEOSadas-5_22_0_p2',  'GEOSadas-5_23_0',
              'GEOSadas-5_23_0_p1',   'GEOSadas-5_24_0',     'GEOSadas-5_24_0_p1' )

    # BCS Tags: Icarus-NLv3 (New Land Parameters)
    #---------------------------------------------------------------------------
     INL  = ( 'INL', 'Icarus-NL', 'Icarus-NLv3', 'Jason-NL' )

     self.bcsTag={}
     for tag in F14:  self.bcsTag[tag]= "Fortuna-1_4"
     for tag in F20:  self.bcsTag[tag]= "Fortuna-2_0"
     for tag in F21:  self.bcsTag[tag]= "Fortuna-2_1"
     for tag in G10:  self.bcsTag[tag]= "Ganymed-1_0"
     for tag in G10p: self.bcsTag[tag]= "Ganymed-1_0_M"
     for tag in G20:  self.bcsTag[tag]= "Ganymed-1_0_M"
     for tag in G30:  self.bcsTag[tag]= "Ganymed-1_0_Reynolds"
     for tag in G40:  self.bcsTag[tag]= "Ganymed-4_0_Reynolds"
     for tag in ICA:  self.bcsTag[tag]= "Icarus_Reynolds"
     for tag in INL:  self.bcsTag[tag]= "Icarus-NLv3_Reynolds"

     
     for tag in D214:  self.bcsTag[tag]= "Fortuna-1_4"
     for tag in D540:  self.bcsTag[tag]= "Fortuna-1_4"
     for tag in D561:  self.bcsTag[tag]= "Fortuna-2_1"
     for tag in D580:  self.bcsTag[tag]= "Ganymed-1_0"
     for tag in D591p: self.bcsTag[tag]= "Ganymed-1_0_M"
     for tag in D5A0:  self.bcsTag[tag]= "Ganymed-1_0_M"
     for tag in D5B0:  self.bcsTag[tag]= "Ganymed-1_0_Reynolds"
     for tag in D512:  self.bcsTag[tag]= "Ganymed-4_0_Reynolds"
     for tag in D517:  self.bcsTag[tag]= "Icarus_Reynolds"
    
     self.tagsRank ={}
     self.tagsRank['Fortuna-1_4'] = 1 
     self.tagsRank['Fortuna-1_5'] = 2 
     self.tagsRank['Fortuna-2_0'] = 3 
     self.tagsRank['Fortuna-2_1'] = 4 
     self.tagsRank['Ganymed-1_0'] = 5 
     self.tagsRank['Ganymed-1_0_m1'] = 6 
     self.tagsRank['Ganymed-1_0_m2'] = 7 
     self.tagsRank['Ganymed-1_0_M']  = 8 
     self.tagsRank['Ganymed-1_0_D']  = 9 
     self.tagsRank['Ganymed-1_0_Reynolds'] = 10
     self.tagsRank['Ganymed-1_0_Ostia']    = 11
     self.tagsRank['Ganymed-4_0_Reynolds'] = 12
     self.tagsRank['Ganymed-4_0_Ostia']    = 13
     self.tagsRank['Ganymed-4_0_MERRA-2']  = 14
     self.tagsRank['Icarus_Reynolds']      = 15
     self.tagsRank['Icarus_MERRA-2']       = 16
     self.tagsRank['Icarus_Ostia']         = 17
     self.tagsRank['Icarus-NLv3_Reynolds'] = 18
     self.tagsRank['Icarus-NLv3_MERRA-2']  = 19
     self.tagsRank['Icarus-NLv3_Ostia']    = 20

     self.bcbase={}
     self.bcbase['discover_ops'] = "/discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/bcs"
     self.bcbase['discover_lt']  = "/discover/nobackup/ltakacs/bcs"

  def get_bcbase(self, opt):
     base = ''
     if opt.upper() == 'IN':
        base = self.common_in.get('bc_base')
     if opt.upper() == 'OUT':
        base = self.common_out.get('bc_base')
     assert base, 'please specify bc_base: discover_ops, discover_lt or an absolute path'
     if base == 'discover_ops' or base == 'discover_lt':
        return self.bcbase[base]
     else:
        return base 

  def get_bcdir(self, opt):
    bc_base = self.get_bcbase(opt)
    tag = self.common_in['tag']
    ogrid = self.common_in['ogrid']
    if opt.upper() == "OUT":
       tag = self.common_out['tag']
       ogrid = self.common_out['ogrid']
    bctag = self.get_bcTag(tag,ogrid)

    tagrank = self.tagsRank[bctag]
    if (tagrank >= self.tagsRank['Icarus-NLv3_Reynolds']) :
       bcdir = bc_base+'/Icarus-NLv3/'+bctag+'/'
    elif (tagrank >= self.tagsRank['Icarus_Reynolds']):
       if bc_base == self.bcbase['discover_ops']:
          bcdir = bc_base+'/Icarus_Updated/'+bctag+'/'
       else:
          bcdir = bc_base+'/Icarus/'+bctag+'/'
    elif(tagrank >= self.tagsRank["Ganymed-4_0_Reynolds"]):
       bcdir = bc_base + '/Ganymed-4_0/'+bctag+'/'
    else:
       bcdir = bc_base + '/' + bctag

    if not os.path.exists(bcdir):
       print( "Cannot find bc dir " +  bcdir)
       exit()

    gridStr = self.get_grid_subdir(bcdir,opt)
    bcdir =  bcdir + '/' + gridStr

    return bcdir
        
  def get_grid_subdir(self, bcdir, opt):

     def get_name_with_grid( grid, names):
       if not grid :
         return names
       namex = []
       if (grid[0].upper() == 'C'):
         n = int(grid[1:])
         s1 =f'{n}x6C'
         j=n*6
         s2 =str(n)
         s3 =str(j)
         # first try
         for name in names:
           if (name.find(s1) != -1):
              namex.append(name)
         if len(namex) ==0:
           for name in names:
             if (name.find(s2) != -1 and name.find(s3) != -1): namex.append(name)
       else:
         xy = grid.upper().split('X')
         s2 = xy[0]
         s3 = xy[1]
         for name in names:
           if (name.find(s2) != -1 and name.find(s3) != -1): namex.append(name)
       return namex

     dirnames = [ f.name for f in os.scandir(bcdir) if f.is_dir()]
     agrid_ = self.common_in['agrid']   
     ogrid_ = self.common_in['ogrid'] 
     if opt.upper() == "OUT" :
       agrid_ = self.common_out['agrid']   
       ogrid_ = self.common_out['ogrid'] 
          
     anames = get_name_with_grid(agrid_, dirnames)
     gridID = get_name_with_grid(ogrid_, anames)
     if len(gridID) ==0 or len(gridID) >=2 :
       print(" cannot find the grid string: " + bcdir)
     return gridID[0]

  def get_bcTag(self, tag, ogrid):
    bctag = self.bcsTag[tag]
    if ogrid[0].upper() == "C": 
       bctag=bctag.replace('_Reynolds','_Ostia')
    else:
       xy = ogrid.upper().split('X')
       x = int(xy[0])
       if x == 1440:   bctag=bctag.replace('_Reynolds','_MERRA-2')
       if x == 2800:  
          bctag=bctag.replace('_Reynolds','_Ostia')
          bctag=bctag.replace('_M','_D')
    return bctag
