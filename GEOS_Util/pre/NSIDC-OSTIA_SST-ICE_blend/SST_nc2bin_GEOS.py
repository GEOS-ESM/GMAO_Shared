from numpy import *
import xarray as xr
import pandas as pd
from pandas.tseries.offsets import Day, BDay
from scipy.io import FortranFile

# this script reads SST/sea ice fraction netcdf data, interpolates,
# writes in fortran bianry to make compaitable with GEOS
# this example uses CMCC-CM2 cmip6 12 monthly mean climatological data

# reading netcdf data
sst=xr.open_dataset('CMCC_CM2_clim_360x180_180.nc').ts

# assumes SST is in -180 to 180 fo lon. If not, needs to regrid. 

lon=arange(-179.5,179.5+1)
lat=arange(-89.5,89.5+1,1)
sst=sst.interp(lon=lon,lat=lat)

# linearly interpolating nan grids created by xarray interp at 180 longitudes
sst[:,:,-1]=(sst[:,:,-2]+sst[:,:,0])/2


# data setup

hlen=14 # header length
nlon=360 # size of lon
nlat=180 # size of lat
nt=12 # size of time dimension

# creating date stamps to write for header file. 
start_date=pd.date_range('1978-12-31', '1980-01-01', freq='M')+ Day(1)

# writing output
f = FortranFile('dataoceanfile_CMIP6_his_.1950-1979.clim.360x180.LE', 'w')
for m in tqdm(range(nt)):
	header=zeros(14,dtype='int')
	header[0]=int(str(start_date[m])[0:4])    # record start year
	header[1]=int(str(start_date[m])[5:7])    # record start mon
	header[2]=int(str(start_date[m])[8:10])   # record start day
	header[3]=0 							  # record start hr

	header[6]=int(str(start_date[m+1])[0:4])  # record end year
	header[7]=int(str(start_date[m+1])[5:7])  # record end mon
	header[8]=int(str(start_date[m+1])[8:10]) # record end day	
	header[9]=0 							  # record end hr

	header[12]= nlon
	header[13]= nlat

	clim_data= (sst[m,:,:].data)
	clim_data=reshape(clim_data,nlon*nlat)

	f.write_record(np.array(header, dtype = np.float32))
	f.write_record(np.array(clim_data, dtype = np.float32))

	print(header)
f.close()

# write fraci

# same method to write fraci (range values from 0 to 1)

#afahad Nov 29, 2021
