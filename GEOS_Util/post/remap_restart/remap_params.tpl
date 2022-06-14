#
# This template file can be filled with questionary or manually
# 
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
    zoom:
    wemin:
    # it supports three models: catch, catchcnclm40, catchcnclm45
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
    split_saltwater: false
    surflay: 20.
    wemin:
  analysis:
    bkg: true
    aqua: False
    lcv: false

slurm:
  account: 
  qos: 
  constraint: 
