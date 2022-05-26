#
# This template file can be filled with questionary or manually
# 
# By default, the values in child dictionary overwrites the values
# in parent dictionary. For example, rst_dir in input is the location with 
# restart files to be regrrided. The surface restart files can be in
# surfface:rst_dir if it has value 
# 

input:
  air:
    drymass: 1
    hydrostatic: 0
  shared:
    agrid: 
    bcs_dir:
    expid:
    ogrid:
    rst_dir:
    yyyymmddhh:
  surface:
    wemin:
    catch_model: null

output:
  shared:
    agrid:
    bcs_dir: 
    expid: 
    ogrid: 
    out_dir:
  air:
    nlevel:
  surface:
    zoom:
    split_saltwater: false
    surflay: 20.
    wemin:
  analysis:
    bkg: true
    aqua: False
    lcv: false

slurm:
  account: 
  debug: 
  partition: 
