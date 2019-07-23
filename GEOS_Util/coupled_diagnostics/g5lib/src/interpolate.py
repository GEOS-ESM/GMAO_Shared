#!/usr/bin/env python
import netCDF4 as nc
import scipy as sp
import glob,sys
from optparse import OptionParser

modulepath='/home/yvikhlia/python'
if not(sys.path.count(modulepath)): sys.path.append(modulepath)
import my_lib.utils as utl

usage = '''
usage: %prog option arg

Reads input collection on tripolar grid and interpolates it to latlon grid. 
Appends variable VAR_interpolated to the output collection. Run %prog -h for help.
'''

parser=OptionParser(usage=usage)
parser.add_option("--idir",help="input collection")
parser.add_option("--odir",help="output collection")
parser.add_option("--tfile",help="tile file")
parser.add_option("--var",help="Var name to interpolate in the input collection")
parser.add_option("--iext",help="File extension to search in the idir and idir", default='nc')
parser.add_option("--oext",help="File extension to search in the idir and odir", default='nc')

(opts, args) = parser.parse_args()

for arg in ['idir','odir','tfile','var']:
    if getattr(opts,arg)==None:
        parser.error('set '+arg+' correctly')


iflist=glob.glob(opts.idir+'/*'+opts.iext); iflist.sort() # Input file list
oflist=glob.glob(opts.odir+'/*'+opts.oext); oflist.sort() # Output file list

if len(iflist)!=len(oflist):
    parser.error('Files in input and output collections do not match.')

# Read tile information
x=sp.loadtxt(opts.tfile, skiprows=8,dtype=sp.float32)
otdata=x[sp.where(x[:,0]==0)[0]]
t=sp.float32
rec=[('iin',t),('jin',t),('iout',t),('jout',t),('frac',t)]
tdata=sp.zeros(otdata.shape[0],dtype=sp.dtype(rec))
tdata['iin']=otdata[:,8]; tdata['jin']=otdata[:,9]
tdata['iout']=otdata[:,4]; tdata['jout']=otdata[:,5]; tdata['frac']=otdata[:,6]
del x,otdata,t

def interpolate(varname,ifilename,ofilename,tiledata):
    
    # Read input field
    fin=nc.Dataset(ifilename)
    varin=fin.variables[varname]; varin.set_auto_maskandscale(True)
    typein=varin.dtype
    dimsin=varin.dimensions
    datain=varin[:]
    shin=datain.shape
    tin=fin.variables['time'].units

    # Read output dimensions
    fout=nc.Dataset(ofilename,'a')
    shout=shin[:-2]+fout.variables['lat'].shape+fout.variables['lon'].shape
    tout=fout.variables['time'].units
    
    if tin!=tout: 
        raise Exception('Input date does not match output date.')
    
    dataout=utl.g2g(datain,tdata,shout[-2:])

    # Put field into output netCDF
    ovarname=varname+'_interpolated'
    varout=fout.createVariable(ovarname,typein,dimsin)
    for attr in varin.ncattrs():
        varout.setncattr(attr,getattr(varin,attr))
    varout[:]=dataout
    fin.close(); fout.close()
    return dataout

for num,fname in enumerate(iflist):
    print 'Input file: ',fname
    print 'Output file: ',oflist[num]
    var=interpolate(opts.var,fname,oflist[num],tdata)

