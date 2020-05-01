#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 10:33:12 2019

@author: smahanam
"""

from bs4 import BeautifulSoup
from netCDF4 import Dataset
import numpy as np
import re, os
import matplotlib.pyplot as plt
from pyhdf.SD import SD, SDC
import numpy.ma as ma
from datetime import datetime, timedelta
import requests
import shutil

#----------------------------------------------------------#
#                 BEGIN USER DEFINED VARIABLES             #

"""
MODIS LAI Product Name (M6_NAME):
    MOD15A2H (Terra   ) 20000218 (DOY : 049) - to date   
    MYD15A2H (Auqa    ) 20020704 (DOY : 185) - to date
    MCD15A2H (combined) 20020704 (DOY : 185) - to date
"""

M6_NAME     = 'MCD15A2H'
process_lai = True
process_alb = False
# Output file dimensions
IM = 720
JM = 360

#                  END USER DEFINED VARIABLES              # 
#----------------------------------------------------------#

# ---- Global Parameters
 
MAPL_UNDEF = np.float(1.e15)
DY = 180. / JM 
DX = 360. / IM
CWD= os.getcwd()
M6_DIR       = {'MCD15A2H':'MOTA', 'MOD15A2H':'MOLT','MYD15A2H':'MOLA'}
MODIS_DOWN   = "https://e4ftl01.cr.usgs.gov/" + M6_DIR[M6_NAME] + '/' + M6_NAME + '.006/'
MODIS_PATH   = "https://ladsweb.modaps.eosdis.nasa.gov/opendap/hyrax/allData/6/" + M6_NAME + '/'
EXTDATA_PATH = "/discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001/"

class DriverFunctions (object):

    def create_netcdf (FILE_NAME,VAR_NAMES):
        import datetime
        ncFidOut = Dataset(FILE_NAME,'w',format='NETCDF4')
        LatDim  = ncFidOut.createDimension('lat', JM)
        LonDim  = ncFidOut.createDimension('lon', IM)
        timeDim = ncFidOut.createDimension('time', None)
        
        ncFidOut.description = "MODIS MCD43GF v006 VISDF and NIRDF @ 1km aggrregated to the 0.5 degree grid"
        ncFidOut.history     = "Created on " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " by Sarith Mahanama (sarith.p.mahanama@nasa.gov)"
        
        # Variables
        lonout  = ncFidOut.createVariable('lon','f4',('lon',))
        latout  = ncFidOut.createVariable('lat','f4',('lat',))
        timeout = ncFidOut.createVariable('time', 'i', ('time',))
        mskout  = ncFidOut.createVariable('mask', 'i', ('lat','lon'))
        for l in range (len(VAR_NAMES)):
            varout = ncFidOut.createVariable(VAR_NAMES[l], 'f4', ('time','lat','lon'), fill_value=np.float(1.e15))
            varout.units  = '1'
            setattr(ncFidOut.variables[VAR_NAMES[l]],'missing_value',np.float32(1.e15))
            setattr(ncFidOut.variables[VAR_NAMES[l]],'fmissing_value',np.float32(1.e15))

        datout  = ncFidOut.createVariable('REFERENCE_DATE','i', ('time',))
        
        # Attributes
        timeout.units = 'days since 2000-01-01 00:00:00'
        latout.units  = 'degrees north'
        lonout.units  = 'degrees east'
            
        varr    = np.full (IM, 0.)
        for i in range (IM):
            varr [i] = -180. + DX/2. + DX*i
        lonout [:] = varr
        
        varr    = np.full (JM,0.)
        for i in range (JM):
            varr [i] = DY*i + -90. + DY/2.  
        latout [:] = varr 
        
        # Mask
        IX = np.array(np.loadtxt (EXTDATA_PATH + 'DE_00720x00360_PE_0720x0360.til',skiprows=8, usecols=4,dtype='int'))
        JX = np.array(np.loadtxt (EXTDATA_PATH + 'DE_00720x00360_PE_0720x0360.til',skiprows=8, usecols=5,dtype='int'))
        geos5_mask    = np.full ((JM, IM), 0)
        
        for i in range (IX.size):
            geos5_mask [JX[i]-1,IX[i]-1] = 1
        mskout[:] = geos5_mask
        ncFidOut.close()  
    
    def nearest_cell (array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def fill_gaps (data, fill_value=None, ocean=None):
        NX = min (data.shape)
        data = ma.masked_array(data,data==fill_value)
        odata= ma.masked_array(data,data==ocean)

        for zoom in range (1,NX//2):
            for direction in (-1,1): 
                shift = direction * zoom
                if not np.any(data.mask): break
                for axis in (0,1):                
                    a_shifted = np.roll(data ,shift=shift,axis=axis)
                    o_shifted = np.roll(odata,shift=shift,axis=axis)
                    idx=~a_shifted.mask * data.mask*~o_shifted.mask
                    data[idx]=a_shifted[idx]
    
    def regrid_to_coarse (data):
        global NC, NR, IM, JM
        JX = NR // JM
        IX = NC // IM
        temp = data.reshape((data.shape[0] // JX, JX, data.shape[1] // IX, IX))
        return np.nanmean(temp, axis=(1,3))

    def get_tag_list(url,label):
        page = requests.get(url)
        soup = BeautifulSoup(page.content, 'html.parser')
        allfiles = soup.find_all(label)
        thislist = [tag.text for tag in allfiles]
        return thislist
        
if process_lai:   

    if not os.path.isfile(EXTDATA_PATH +'MCD15A2H.006_LAI_ExtData.nc4'):
        DriverFunctions.create_netcdf(EXTDATA_PATH +'MCD15A2H.006_LAI_ExtData.nc4',['MODIS_LAI'])
    
    LAI_PATH = '/l_data/model_parameters/LAI/MODIS6/MCD15A2H.006/'
    NC = 86400
    NR = 43200
    DXY = 360./NC
    N_MODIS = 2400  
    
    class MCD15A2H (object):    
    
        #- read MCD15A2H granules    
        
        def __init__ (self, FILE_NAME):
    
            hdf = SD(FILE_NAME, SDC.READ)
    
            # Read datasets
            #--------------
            #- Query global attributes
    
            sattr = str(hdf.attributes())
    
            xdim = int(re.findall(r"XDim=(\d+)", sattr)[0])
            ydim = int(re.findall(r"YDim=(\d+)", sattr)[0])
            
            earth_radius = np.double( re.findall(r"ProjParams=\((-?\d+\.\d*|\d*\.\d+)", sattr)[0] )
            upleft = np.double( re.findall(r"UpperLeftPointMtrs=\((-?\d+\.\d*|\d*\.\d+),(-?\d+\.\d*|\d*\.\d+)", sattr)[0] )
            lowright = np.double( re.findall(r"LowerRightMtrs=\((-?\d+\.\d*|\d*\.\d+),(-?\d+\.\d*|\d*\.\d+)", sattr)[0] )
            # upleft = np.double( re.findall(r"UpperLeftPointMtrs=\((-?\d+\.\d*|\d*\.\d+),-?(\d+\.\d*|\d*\.\d+)", sattr)[0] )
            # lowright = np.double( re.findall(r"LowerRightMtrs=\((-?\d+\.\d*|\d*\.\d+),-?(\d+\.\d*|\d*\.\d+)", sattr)[0] )
            
            #- MODIS pixel geolocation
            #- x,y grid
            dx = (lowright[0] - upleft[0]) / xdim
            dy = (upleft[1] - lowright[1]) / ydim
            
            xm, ym = np.meshgrid(np.arange(1, xdim+1), np.arange(1, ydim+1))
            x_sin = (xm-xm+1)*(upleft[0] + dx/2.) + (xm-1)*dx
            y_sin = (ym-ym+1)*(upleft[1] - dy/2.) - (ym-1)*dy
            
            #- Create Longitude and Latitude maps from Sinusoidal projection
            lat = (y_sin / earth_radius * 180.0/np.pi)
            lon = (x_sin / (earth_radius * np.cos(y_sin/earth_radius)) * 180.0/np.pi)
            eastern_hemispehere = np.where(lon < -180.)
            western_hemispehere = np.where(lon >  180.)
            lon [eastern_hemispehere] = lon [eastern_hemispehere] + 360.
            lon [western_hemispehere] = lon [western_hemispehere] - 360.
            self.x_index = np.array(np.floor ((lon + 180.)/DXY),dtype=np.int)         
            self.y_index = np.array(np.floor ((lat + 90.)/DXY) ,dtype=np.int)
    #        self.earth_radius = earth_radius
            
            # Read LAI
            # --------
    
            DATAFIELD_NAME = 'Lai_500m'
            data2D = hdf.select(DATAFIELD_NAME)
            data = data2D[:,:].astype(np.double)
    
            # Read attributes.
            # 255 = _Fillvalue, assigned when:
            #   * the MOD09GA suf. reflectance for channel VIS, NIR was assigned its _Fillvalue, or
            #   * land cover pixel itself was assigned _Fillvalus 255 or 254.
            # 254 = land cover assigned as perennial salt or inland fresh water.
            # 253 = land cover assigned as barren, sparse vegetation (rock, tundra, desert.)
            # 252 = land cover assigned as perennial snow, ice.
            # 251 = land cover assigned as "permanent" wetlands/inundated marshlands.
            # 250 = land cover assigned as urban/built-up.
            # 249 = land cover assigned as "unclassified" or not able to determine.
            # 248 = no standard deviation available, pixel produced using backup method.
            # 0 <= modis_data <= 100 MODIS good LAI data        
            # -----------------------------------------------------------------------------
    
            attrs = data2D.attributes(full=1)
            #lna=attrs["long_name"]
            #self.name = lna[0]
            vra=attrs["valid_range"]
            valid_range = vra[0]
            fva=attrs["_FillValue"]
            _FillValue = fva[0]
            sfa=attrs["scale_factor"]
            scale_factor = sfa[0]
            
            data [data == 250] = 5
            data [data == 253] = 0.1
            data [np.where (data == 248)] = _FillValue
            data [np.where (data == 249)] = _FillValue
            data [np.where (data == 252)] = _FillValue
            data [np.where (data == 254)] = _FillValue
            
            # Apply MODIS QC 
            DATAFIELD_NAME = 'FparLai_QC'
            data2D = hdf.select(DATAFIELD_NAME)
            dataQC = data2D[:,:].astype(np.int)
            # 00XXXXX0 
            dataQC [np.where (dataQC % 2 ==1)] = -9999
            dataQC [np.where (dataQC > 62)] = -9999
            data [np.where (dataQC == -9999)] = _FillValue
            
    #        DATAFIELD_NAME = 'FparExtra_QC'
    #        data2D  = hdf.select(DATAFIELD_NAME)
    #        dataQCX = data2D[:,:].astype(np.int)
                            
            # Fill marshland (251) with nearest land neighbor
            DriverFunctions.fill_gaps (data, fill_value=251, ocean=_FillValue)
            invalid = data == _FillValue
            invalid = np.logical_or(invalid, data < valid_range[0])
            invalid = np.logical_or(invalid, data > valid_range[1])
            data[invalid] = np.nan
    
            self.lai_500 =  data * scale_factor
            hdf.end()    
    os.chdir(CWD)
    ncFidOut = Dataset(EXTDATA_PATH +'MCD15A2H.006_LAI_ExtData.nc4',mode='a')
    IM = np.array (ncFidOut.variables['lon'][:]).size
    JM = np.array (ncFidOut.variables['lat'][:]).size    
    datestamp = ncFidOut.variables['REFERENCE_DATE']
    LAIOUT    = ncFidOut.variables['MODIS_LAI']
    mask      = np.array (ncFidOut.variables['mask'][:])
    timeout   = ncFidOut.variables['time']
    noland    = mask == 0
    firstdate = np.str(np.array(datestamp)[ 0])
    lastdate  = np.array(datestamp)[-1 -3] -4
    date0     = datetime(int(firstdate[0:4]), int(firstdate[4:6]), int(firstdate[6:8]))
    next_day  = np.array(datestamp)[:].shape[0] - 3

    # Processomh 8-day composites

    years = DriverFunctions.get_tag_list(MODIS_PATH,"a")[1:-6]

    for year in years:
        doys  = DriverFunctions.get_tag_list(MODIS_PATH + year,"a")[1:-6]
        for doy in doys:
            this_date = datetime(int(year[0:4]),1,1) + timedelta (days=int(doy[0:3])-1)
            yyyymmdd  =int(this_date.strftime(('%Y%m%d')))
            if yyyymmdd > lastdate:
                files    = DriverFunctions.get_tag_list(MODIS_PATH + year + doy,"span")[1:-1]
                date1 = datetime(int(year[0:4]),1,1) + timedelta (days=int(doy[0:3])-1) 
                date2 = date1 + timedelta (days=8)
                mday  = date1 + (date2 - date1)/2
                os.chdir(CWD)
                DATADIR = CWD + '/download/' + year + doy + '/'
                if not os.path.exists(DATADIR):
                    os.makedirs(DATADIR)
                    os.chdir(DATADIR)
                    for f in range(len(files)):
                        FILE_NAME = MODIS_DOWN + date1.strftime('%Y.%m.%d') + '/' + files[f]
                        os.system('wget --user=sarith --password=GusBoy2017 --max-redirect=20 ' + FILE_NAME)
                    os.chdir(CWD)

                # ---- Stitch MODIS granules
            
                print ('Processing Date :', mday.strftime("%Y-%m-%d"))
                                
                lai_high = np.full((NR,NC),np.nan)
            
                for f in range(len(files)):   
                    FILE_NAME = DATADIR + files[f]
                    #LAI_PATH + laidirs[next_day] + '/' + files[f]
                    print(FILE_NAME)
                    thistile = MCD15A2H (FILE_NAME)
                    lai_mask = np.where((thistile.lai_500 >= 0.) & (thistile.lai_500 <= 10.))
                    lai_array= thistile.lai_500.reshape (N_MODIS*N_MODIS)
                    xin_array= thistile.x_index.reshape (N_MODIS*N_MODIS)
                    yin_array= thistile.y_index.reshape (N_MODIS*N_MODIS)
                    lai_mask = np.where((lai_array >= 0.) & (lai_array <= 10.))
                    lai_high [yin_array[lai_mask],xin_array[lai_mask]] = lai_array[lai_mask]
        
                # regrid to IMxJM
                lai_low   = DriverFunctions.regrid_to_coarse(lai_high)
                invalid   = np.ma.masked_invalid(lai_low)    
                #    lai_low [invalid.mask] = -9999.
                lai_low [invalid.mask] = 0.
                lai_low [noland] = MAPL_UNDEF
                #    DriverFunctions.fill_gaps (lai_low, fill_value=-9999., ocean=MAPL_UNDEF)
                LAIOUT   [next_day] = lai_low
                datestamp[next_day] = int(mday.strftime("%Y%m%d")) 
                timeout  [next_day] = 1440 * (mday - date0).days
                next_day = next_day + 1
                shutil.rmtree(DATADIR)

    # Add 3 MODIS dates (24 calendar days)
    for ext_day in range(3):
        mday = mday + timedelta (days=8)
        datestamp[next_day] = int(mday.strftime("%Y%m%d")) 
        timeout  [next_day] = 1440 * (mday - date0).days
        LAIOUT[next_day] = lai_low
        next_day = next_day + 1
    ncFidOut.close()
    
if process_alb:   

    if not os.path.isfile(EXTDATA_PATH + 'MCD43GF.006_ALBEDO_ExtData.nc4'):
        DriverFunctions.create_netcdf(EXTDATA_PATH +'MCD43GF.006_ALBEDO_ExtData.nc4',['MODIS_VISDF', 'MODIS_NIRDF'])

    from datetime import datetime
    
    ALB_PATH = '/l_data/model_parameters/MODIS-Albedo/MCD43GF.006/'
    NC = 43200
    NR = 21600
    DXY = 360./NC
    
    class MCD43GF (object):    
    
        #- read MCD43GF global files    
        def read_alb (self, alb_file, var_name):

            #hdf_qa  = SD(qa_file,  SDC.READ) 
            hdf_alb = SD(alb_file, SDC.READ) 
            
            # Read MODIS VISDF and NIRDF
            # --------------------------
            # QC data
            # -------
            #DATAFIELD_NAME = 'QA_' + var_name 
            #data2D = hdf_qa.select(DATAFIELD_NAME)
            #qa = data2D[:,:].astype(np.double)
            
            # Albedo data
            # -----------            
            DATAFIELD_NAME = var_name
            data2D = hdf_alb.select(DATAFIELD_NAME)
            alb    = data2D[:,:].astype(np.double)            
            attrs = data2D.attributes(full=1)
            vra=attrs["valid_range"]
            valid_range = vra[0]
            fva=attrs["_FillValue"]
            _FillValue = fva[0]                        
            sfa=attrs["scale_factor"]
            scale_factor = sfa[0]

            # alb [np.where (qa > 5)] = _FillValue
            invalid = alb == _FillValue
            invalid = np.logical_or(invalid, alb < valid_range[1])
            invalid = np.logical_or(invalid, alb > valid_range[0])
            alb[invalid] = np.nan 
            alb = alb * scale_factor 
            alb = np.flip (alb,0)
            
            # regrid to IMxJM
            alb_low = DriverFunctions.regrid_to_coarse(alb)
            # hdf_qa.end() 
            hdf_alb.end()            
            return alb_low
        
        def __init__ (self, ALB_DIR):

            albfiles = sorted(os.listdir(ALB_DIR))                        
            self.nir = self.read_alb (ALB_DIR + '/' + albfiles[0], 'Albedo_Map_0.7-5.0') 
            self.vis = self.read_alb (ALB_DIR + '/' + albfiles[1], 'Albedo_Map_0.3-0.7') 
    
    ncFidOut = Dataset('ExtData/MCD43GF.006_ALBEDO_ExtData.nc4',mode='a')
    IM = np.array (ncFidOut.variables['lon'][:]).size
    JM = np.array (ncFidOut.variables['lat'][:]).size    
    datestamp = ncFidOut.variables['REFERENCE_DATE']
    VISOUT    = ncFidOut.variables['MODIS_VISDF']
    NIROUT    = ncFidOut.variables['MODIS_NIRDF']
    mask      = np.array (ncFidOut.variables['mask'][:])
    timeout   = ncFidOut.variables['time']
    noland    = mask == 0
    
    # Process dailies
                    
    albdirs = sorted(os.listdir(ALB_PATH))    
    albdirs.remove('wget_mod.csh')
    albdirs.remove('wget_qa.csh')
    
    for d in range(len(albdirs)-1):

        # Processing MODIS Date
        date1 = datetime(int(albdirs[d][0:4]), int(albdirs[d][5:7]), int(albdirs[d][8:10]))
        date2 = datetime(int(albdirs[d+1][0:4]), int(albdirs[d+1][5:7]), int(albdirs[d+1][8:10]))
        mday  = date1 + (date2 - date1)/2
        if d == 0:
            date0 = mday
            timeout.units = 'days since ' + mday.strftime("%Y-%m-%d %H:%M:%S")
        datestamp[d] = int(mday.strftime("%Y%m%d")) 
        timeout[d]   = (mday - date0).days
        print ('Processing Date :', mday.strftime("%Y-%m-%d"))
        thisdate = MCD43GF (ALB_PATH + albdirs[d])                
        nir_low   = thisdate.nir
        vis_low   = thisdate.vis

        invalid   = np.ma.masked_invalid(nir_low)    
    #    nir_low [invalid.mask] = -9999.
        nir_low [invalid.mask] = 0.
        nir_low [noland] = MAPL_UNDEF
    #    DriverFunctions.fill_gaps (nir_low, fill_value=-9999., ocean=MAPL_UNDEF)
        NIROUT[d] = nir_low

        invalid   = np.ma.masked_invalid(vis_low)    
    #    vis_low [invalid.mask] = -9999.
        vis_low [invalid.mask] = 0.
        vis_low [noland] = MAPL_UNDEF
    #    DriverFunctions.fill_gaps (vis_low, fill_value=-9999., ocean=MAPL_UNDEF)
        VISOUT[d] = vis_low     
    ncFidOut.close()
    
