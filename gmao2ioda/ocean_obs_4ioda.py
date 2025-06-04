#! /usr/bin/env python
# /////////////////////////////////////////////////////////////////////////
# /**
# * Title: ocean_obs.py
# *
# * Description: Extract the ocean observations for a specific date/time window.
# *
# * Args:
# *              yyyy     : 4 digit year
# *              mm       : 2 digit month
# *              dd       : 2 digit day
# *              hh       : 2 digit hour
# *              obs_type : Entered as a list
# *                         Argo
# *                         CTD
# *                         XBT
# *                         TAO
# *                         PIRATA
# *                         RAMA
# *                         SMOS (Level 2)
# *                         SMOSSUB (SUbset Level 2)
# *                         SMOSL3 (Level 3)
# *                         AQUARIUS (Level 2)
# *                         SMAP (Level 2)
# *                         SMAPV4.1 (Level 2)
# *                         SMAPV4.2 (Level 2)
# *                         SMAPV4.3 (Level 2)
# *                         AQUARIUSL3
# *                         SMAPL3
# *                         CryoSat-2
# *                         Jason-1
# *                         Jason-2
# *                         Jason-3
# *                         Saral
# *                         Sentinel-3a
# *                         ERS-1
# *                         ERS-2
# *                         TOPEX
# *                         GEOSAT-2
# *                         Envisat
# *                         HY-2A
# *                         Reynolds
# *                         OSTIA
# *                         AVHRR-18
# *                         NOAA-16
# *                         METOP-A
# *                         NASA-TEAM-2
# *                         NOAA-AICE
# *                         OIB-HICE
# *
# * Example:
# *
# *             ocean_obs.py 2012 03 15 03 Argo CTD Jason-2 CryoSat-2
# *
# *             Will create a netcdf file and a series of png figures containing/showing the observations
# *             that are within the assimilation window.
# *
# * @author: Guillaume Vernieres
# */
#
# /////////////////////////////////////////////////////////////////////////
# Date: Dec 2015

import matplotlib
matplotlib.use('agg')
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys
import datetime
import os
import ocean_obs_utils
from ocean_obs_utils import cs2_reader

#   test of 50 layers Argo data for 2012 only
OBSDIR='/discover/nobackup/projects/gmao/ssd/g5odas/production/GEOS5odas-5.00/RC/OBS/'
LEV50INSITUOBSDIR='/discover/nobackup/projects/gmao/m2oasf/aogcm/g5odas/obs/assim/'
ERICLEV50INSITUOBSDIR='/discover/nobackup/projects/gmao/ssd/decadal/ehackert/obs/assim/'
ADTOBSDIR='/discover/nobackup/projects/gmao/m2oasf/aogcm/g5odas/obs/assim/AVISO/V3/'
SatSSSOBSDIR='/discover/nobackup/projects/gmao/m2oasf/aogcm/g5odas/obs/assim/'
SYNOBSDIR='/discover/nobackup/projects/gmao/m2oasf/aogcm/g5odas/obs/assim/SYN_7.0/'


#   artificially ramp up smos error
smosrampuperr=1.

# OIB data
OIBDIR='/discover/nobackup/projects/gmao/ssd/g5odas/gvernier/SAND_BOXES/sea_ice_da/ICE_BRIDGE/'

# NESDIS L2-SST
NESDISDIR='/discover/nobackup/projects/gmao/m2oasf/aogcm/g5odas/obs/raw/L2_SST/NESDIS/'
prefix_dict = {
    'NOAA16'  : '0000-STAR-L2P_GHRSST-SSTskin-AVHRR16_G-ACSPO_V2.40-v02.0-fv01.0.nc',
    'METOPA'  : '0000-STAR-L2P_GHRSST-SSTskin-AVHRRMTA_G-ACSPO_V2.40-v02.0-fv01.0.nc'
    }

# Note on sigos from ORAS4 from 2012 tech-memo "The NEMOVAR ocean data assimilation system ..."
# Mostly errors of representation, so same for XBT, CTD and Argo (while measurement errors is quite different)
# T sigo: 0.18 deg C, single global surface value but depth dependent, lower bound is 0.07 deg C
# S sigo: 0.18 psu, single global surface value but depth dependent, lower bound is 0.02 psu
# SSH: 5 cm

T_prof_sigo = 0.25   # deg C
#S_prof_sigo = 0.18  # psu
S_prof_sigo = 0.025  # psu
ADT_sigo = 0.2
SSS_sigo = 1.0         #psu


xsigo_s    = 1.0 #1.0
xsigo_t    = 1.0 #1.0
xsigo_sst  = 2.0
xsigo_sss  = 1.0
xsigo_ssh  = 0.25
xsigo_aice = 1.0
xsigo_hice = 0.01

try:
    EXP_NDAYS = float(os.environ['ODAS_NDAYS'])
except:
    EXP_NDAYS = 0.125 # Default value for the observation window [yyyymmdd_hh-EXP_NDAY
    #EXP_NDAYS = 20.125 # Default value for the observation window [yyyymmdd_hh-EXP_NDAY

logit_transform = False
use_obs_old = False

try:
    T_prof_sigo = float(os.environ['ODAS_T_prof_sigo'])
    S_prof_sigo = float(os.environ['ODAS_S_prof_sigo'])
    ADT_sigo = float(os.environ['ODAS_ADT_sigo'])
    SSS_sigo = float(os.environ['ODAS_SSS_sigo'])

    xsigo_t    = T_prof_sigo
    xsigo_s    = S_prof_sigo
    xsigo_ssh  = ADT_sigo
    xsigo_sss  = SSS_sigo

    SCRDIR = os.environ['SCRDIR']
    EXPDIR = os.environ['EXPDIR']
    EXPID = os.environ['EXPID']
    #logit_transform = os.environ['ODAS_logit_transform']
    #if (logit_transform=='True'):
    #    logit_transform = True
    #else:
    #    logit_transform = False
except:
    print('Environement variables not set, reverting to default:')


print('NDAYS=',EXP_NDAYS)
print('T_prof_sigo=',T_prof_sigo)
print('S_prof_sigo=',S_prof_sigo)
print('ADT_sigo=',ADT_sigo)
print('SSS_sigo=',SSS_sigo)

# Obs id as defined in the UMD_oletkf
obsid_dict = {
    'id_s_obs'   : 5521, # Salinity
    'id_t_obs'   : 3073, # Temperature
    'id_sst_obs' : 5525, # SST
    'id_sss_obs' : 5522, # Sea surface Salinity
    'id_ssh_obs' : 5526, # SSH (Not used ...)
    'id_eta_obs' : 5351, # SSH
    'id_aice_obs': 6000, # AICE
    'id_hice_obs': 6001  # HICE
    }

def inv_logit(p):
    return np.exp(p) / (1 + np.exp(p))

def da_window(yyyy, mm, dd, hh, NDAYS, OBSDATE, QCPRF):
    """
    This function exctract the observations that are in [yyyymmdd_hh-NDAYS, yyyymmdd_hh+NDAYS]
    Args:
        yyyy (int)     : 4 digit year
          mm (int)     : 2 digit month
          dd (int)     : 2 digit day
          hh (int)     : 2 digit hour
          NDAYS (float): Used to define the size of the assimilation window centered at yyyymmddhh [yyyymmddhh-NDAYS, yyyymmddhh+NDAYS].
          OBSDATE (int): Observation date in a YYYYMMDDHH format
          QCPRF (int)  : Quality control flag 1=good, 0=bad

    Returns:
        Array of indices that corresponds to the observations that are within the window and a QCPRF flag of 1
    """

    center = datetime.datetime(int(yyyy), int(mm), int(dd), int(hh))

    #Lower bound for obs time
#    obs_date_min=datetime.datetime(int(yyyy), int(mm), int(dd), int(hh))-datetime.timedelta(days=NDAYS)
    obs_date_min=(center-datetime.timedelta(days=NDAYS)).strftime('%Y%m%d%H')

#    yyyyo=str(obs_date_min.year)
#    mmo=str(obs_date_min.month).zfill(2)
#    ddo=str(obs_date_min.day).zfill(2)
#    hho=str(obs_date_min.hour).zfill(2)

    #Upper bound for obs time
#    obs_date_max=datetime.datetime(int(yyyy), int(mm), int(dd), int(hh))+datetime.timedelta(days=NDAYS)
    obs_date_max=(center+datetime.timedelta(days=NDAYS)).strftime('%Y%m%d%H')

#    yyyye=str(obs_date_max.year)
#    mme=str(obs_date_max.month).zfill(2)
#    dde=str(obs_date_max.day).zfill(2)
#    hhe=str(obs_date_max.hour).zfill(2)

    return list(np.squeeze(np.where( (OBSDATE>=int(obs_date_min)) & (OBSDATE<int(obs_date_max)) & (QCPRF==1))))


def readnc(fname, varname):
    ncfile=Dataset(fname)
    VAR=np.squeeze(ncfile.variables[varname][:])
    print(varname, np.shape(VAR))
    ncfile.close()
    return VAR

def standard_obs_reader(fname, vartype):
    print('IN standard_obs_reader',fname)
    ncfile=Dataset(fname)
    N_LEVS    = len(ncfile.dimensions['N_LEVS'])
    N_LEVS    = min(N_LEVS, 50)                               # Assumes profiles have been superobed to 50 levels
    DEPTH     = ncfile.variables['DEPTH'][:]
    VAR       = ncfile.variables[vartype][:]
    QC_LEV    = ncfile.variables['QC_LEV'][:]
    QC_PRF    = ncfile.variables['QC_PRF'][:]
    LON       = ncfile.variables['LON'][:]
    LAT       = ncfile.variables['LAT'][:]
    DATE_TIME = ncfile.variables['DATE_TIME'][:]
    OBS_ERROR = ncfile.variables['OBS_ERROR'][:]
    INST_ID   = ncfile.variables['INST_ID'][:]
    ncfile.close()

    print(np.shape(LON), np.shape(VAR))

    return N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID

def l2_smossub_reader(fname, vartype,platform='SMOSSUB'):
    print('IN l2_smossub_reader')
    ncfile=Dataset(fname)
    N_LEVS    = len(ncfile.dimensions['N_LEVS'])
    N_LEVS    = min(N_LEVS, 50)                               # Assumes profiles have been superobed to 50 levels
    DEPTH     = ncfile.variables['DEPTH'][:]
    VAR       = ncfile.variables[vartype][:]
    QC_LEV    = ncfile.variables['QC_LEV'][:]
    QC_PRF    = ncfile.variables['QC_PRF'][:]
    LON       = ncfile.variables['LON'][:]
    LAT       = ncfile.variables['LAT'][:]
    DATE_TIME = ncfile.variables['DATE_TIME'][:]
    OBS_ERROR = ncfile.variables['OBS_ERROR'][:]
    INST_ID   = ncfile.variables['INST_ID'][:]
    ncfile.close()

    print(vartype, np.shape(LON),np.shape(VAR))

    return N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID

def M2_sst_reader(yyyy, mm, dd, hh, path2scratch, path2expdir, path2expid):

    # Reads MERRA-2 sst from path_to_scratch
    LON = []
    LAT = []
    VAR = []
    DATE_TIME = []
    INST_ID = []
    QC_FLAG = []
    QC_PRF = []
    NPTS = []
    DEPTH = []
    QC_LEV = []
    OBS_ERROR = []
    DATA_ID = []
    N_PROF = []


    fname=path2scratch+'/sst_'+yyyy+mm+dd+'_'+hh+'00z.nc'
    fname2=path2scratch+'/AICE_'+yyyy+mm+dd+'_'+hh+'00z.nc'
#   this change is for ERIC only 12/13/22
    fname3=path2scratch+'/'+path2expid+'.ice_inst_6hr_glo_T1440x1080_slv.'+yyyy+mm+dd+'_1200z.nc4'
#   fname3=path2scratch+'/'+path2expid+'.geosgcm_seaice2.'+yyyy+mm+dd+'_1200z.nc4'
#   ddm1 = str(int(dd)-1).zfill(2)
#   fname3=path2scratch+'/'+path2expid+'.geosgcm_seaice.'+yyyy+mm+ddm1+'_1200z.nc4'
#   fname3=path2scratch+'/'+path2expid+'.geosgcm_seaice2.'+yyyy+mm+dd+'_1200z.nc4'
    print (fname)
    print (fname2)
    print (fname3)
    if (os.path.exists(fname)):
        ncf = Dataset(fname, 'r')
        sst = (np.squeeze(ncf.variables['sst'][:]))
        omask = (np.squeeze(ncf.variables['mask'][:]))
        ncf.close()
        print('in M2_sst_reader past data read')

        if (os.path.exists(fname2)):
            ncf = Dataset(fname2, 'r')
            sic = (np.squeeze(ncf.variables['AICE'][:]))
            ncf.close()
            print('in M2_sst_reader past AICE data read')
            print('length of sic',sic.shape)
        else:
            print('MISSING AICE FILE',fname2)
            sys.exit(1)

        if (os.path.exists(fname3)):
            ncf = Dataset(fname3, 'r')
            sicmod = (np.squeeze(ncf.variables['AICE'][:]))
            ncf.close()
            print('in M2_sst_reader past AICE model data read')
            print('length of sicmod',sicmod.shape)
        else:
            print('MISSING AICE MODEL FILE',fname3)
            sys.exit(1)

        #read mom's grid
        fname=path2scratch+'/grid_spec.nc'
        print('grid file:',fname)
        ncf = Dataset(fname, 'r')
        xt = ncf.variables['x_T'][:]
        yt = ncf.variables['y_T'][:]
        ncf.close()

        #Inline vars
        xt=xt[:]
        yt=yt[:]
        sst=sst[:]
        sic=sic[:]
        sicmod=sicmod[:]
        omask=omask[:]
        print('length of sic, sicmod, sst:',sic.shape,sicmod.shape,sst.shape)


        #Get rid of land values
        I=np.where(omask==1.0)
        xt=xt[I]
        yt=yt[I]
        sst=sst[I]
        sic=sic[I]
        sicmod=sicmod[I]
        print('length of sic,sicmod,sst after omask',sicmod.shape,sic.shape,sst.shape)

        #Get rid of SST under ice
        I=np.where(sic==0.0)
        xt=xt[I]
        yt=yt[I]
        sst=sst[I]
        sic=sic[I]
        sicmod=sicmod[I]

        #Get rid of SST under model ice
        I=np.where(sicmod==0.0)
        xt=xt[I]
        yt=yt[I]
        sst=sst[I]
        print ('length of sst is after ICE is',len(sst))

        xt[xt<-180.0]=xt[xt<-180.0]+360.0

        #subsample
        skip=4
#       skip=8
#       skip=16
        xt=xt[0::skip]
        yt=yt[0::skip]
        sst=sst[0::skip]

        N=len(sst)
        LON = xt
        LAT = yt
        VAR = sst
        yyyymmddhh = yyyy+mm+dd+hh
        DATE_TIME = int(yyyymmddhh)*np.ones(N)
        INST_ID = 516*np.ones(N)
        QC_FLAG = np.ones(N)
        QC_PRF = np.ones(N)
        DATA_ID = np.ones(N)
        N_PROF = np.ones(N)

        NPTS = np.ones(N)
        DEPTH = np.ones( (N,1) )
        QC_LEV = np.ones( (N,1) )
#       OBS_ERROR = 0.5*np.ones( (N,1) )
#       OBS_ERROR = 0.05*np.ones( (N,1) )
        OBS_ERROR = 0.20*np.ones( (N,1) )

        VAR = np.reshape(VAR, (len(VAR),1))

    N_LEVS = 1

    return N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID

def l2_gmi_reader(yyyy, mm, dd, hh, NDAYS=0.125):

    import read_l2_gmi_rss

    LON = []
    LAT = []
    VAR = []
    DATE_TIME = []
    INST_ID = []
    QC_FLAG = []
    QC_PRF = []
    NPTS = []
    DEPTH = []
    QC_LEV = []
    OBS_ERROR = []
    DATA_ID = []
    N_PROF = []
    obs_cnt = 0
    cnt = 0

    lat, lon, time_, sst, qc, sigo_ = read_l2_gmi_rss.gather_l2_gmi_rss(int(yyyy), int(mm), int(dd), int(hh), NDAYS = NDAYS, plot = False )

    obs_cnt = obs_cnt + np.shape(sst)[0]

    N=len(sst)
    LON = np.append( LON, lon[:] )
    LAT = np.append( LAT, lat[:] )
    VAR = np.append( VAR, sst[:] )
    yyyymmddhh = yyyy+mm+dd+hh
    DATE_TIME = np.append( DATE_TIME, int(yyyymmddhh)*np.ones(N) )
    INST_ID = np.append( INST_ID, 516*np.ones(N) )
    QC_FLAG = np.append( QC_FLAG, np.ones(N) )
    QC_PRF = np.append( QC_PRF, np.ones(N) )
    DATA_ID = np.append( DATA_ID, np.ones(N) )
    N_PROF = np.append( N_PROF, np.ones(N) )

    NPTS = np.append( NPTS, np.ones(N) )
    DEPTH = np.append( DEPTH, np.ones( (N,1) ) )
    QC_LEV = np.append( QC_LEV, np.ones( (N,1) ) )
    OBS_ERROR = np.append( OBS_ERROR, 0.5*np.ones( (N,1) ) )

    cnt+=N
    N_LEVS = 1

    VAR = np.reshape(VAR, (len(VAR),1))

    return N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID

def l2_sst_reader(yyyy, mm, dd, hh, platform='NOAA16', NDAYS=0.125):

    # Select L2 files
    center = datetime.datetime(int(yyyy), int(mm), int(dd), int(hh))
    date_end = center + datetime.timedelta(days=NDAYS)
    doy_end = (date_end - datetime.datetime(int(date_end.year), 1, 1)).days + 1

    date_start = center - datetime.timedelta(days=NDAYS)
    doy_start = (date_start - datetime.datetime(int(date_start.year), 1, 1)).days + 1

    prefix = prefix_dict[platform] #='0000-STAR-L2P_GHRSST-SSTskin-AVHRR16_G-ACSPO_V2.40-v02.0-fv01.0.nc'

    Nhrs=int(NDAYS*24*2.0)
    print('Nhrs:',Nhrs)
    list_of_dates = list(range(Nhrs))
    list_of_files = list(range(Nhrs))
    for n in range(Nhrs):
        if n == 0:
            list_of_dates[0] = date_start
        else:
            list_of_dates[n] = list_of_dates[n-1] + datetime.timedelta(hours=1)

        doy   = (list_of_dates[n] - datetime.datetime(int(list_of_dates[n].year), 1, 1)).days + 1
        year  = str(list_of_dates[n].year)
        month = str(list_of_dates[n].month).zfill(2)
        day   = str(list_of_dates[n].day).zfill(2)
        hour  = str(list_of_dates[n].hour).zfill(2)
        list_of_files[n] = NESDISDIR+'/'+platform+'/'+str(list_of_dates[n].year)+'/'+str(doy).zfill(3)+'/'+year+month+day+hour+prefix

    LON = []
    LAT = []
    VAR = []
    DATE_TIME = []
    INST_ID = []
    QC_FLAG = []
    QC_PRF = []
    NPTS = []
    DEPTH = []
    QC_LEV = []
    OBS_ERROR = []
    DATA_ID = []
    N_PROF = []
    obs_cnt = 0
    cnt = 0
    for fname in list_of_files:
        ncf = Dataset(fname, 'r')
        #!!!!!!!!!!!!!!!!!!!!!!!!!
        #Need to also read bias and std
        #!!!!!!!!!!!!!!!!!!!!!!!!!
        sst = (np.squeeze(ncf.variables['sea_surface_temperature'][:]))-273.15
        qc  = (np.squeeze(ncf.variables['quality_level'][:]))
        lon = np.squeeze(ncf.variables['lon'][:])
        lat = np.squeeze(ncf.variables['lat'][:])
        yyyymmddhh = ncf.time_coverage_start[0:8]+ncf.time_coverage_start[9:11]
        ncf.close()

        sst[qc<5]=np.nan
        sst=sst.flatten()
        I=np.where(np.isfinite(sst))
        sst=sst[I]

        lon = lon.flatten()
        lon=lon[I]

        lat = lat.flatten()
        lat=lat[I]

        obs_cnt = obs_cnt + np.shape(sst)[0]

        N=len(sst)
        LON = np.append( LON, lon[:] )
        LAT = np.append( LAT, lat[:] )
        VAR = np.append( VAR, sst[:] )
        DATE_TIME = np.append( DATE_TIME, int(yyyymmddhh)*np.ones(N) )
        INST_ID = np.append( INST_ID, 516*np.ones(N) )
        QC_FLAG = np.append( QC_FLAG, np.ones(N) )
        QC_PRF = np.append( QC_PRF, np.ones(N) )
        DATA_ID = np.append( DATA_ID, np.ones(N) )
        N_PROF = np.append( N_PROF, np.ones(N) )

        NPTS = np.append( NPTS, np.ones(N) )
        DEPTH = np.append( DEPTH, np.ones( (N,1) ) )
        QC_LEV = np.append( QC_LEV, np.ones( (N,1) ) )
        OBS_ERROR = np.append( OBS_ERROR, 0.5*np.ones( (N,1) ) )

        cnt+=N
    N_LEVS = 1

    VAR = np.reshape(VAR, (len(VAR),1))

    return N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID

def oib_reader(yyyy, mm, dd, hh, platform='AIR-BORN', NDAYS=3.0):

    fname = OIBDIR+'IDCSI4_'+yyyy+mm+dd+'.txt'

    LON = []
    LAT = []
    VAR = []
    DATE_TIME = []
    INST_ID = []
    QC_FLAG = []
    QC_PRF = []
    NPTS = []
    DEPTH = []
    QC_LEV = []
    OBS_ERROR = []
    DATA_ID = []
    N_PROF = []
    obs_cnt = 0
    cnt = 0

    if (hh!='12'):
        return N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID

    f = open(fname, 'r')
    for line in f:
        try:
            columns = line.split()
            lat=float(columns[0][0:-2])
            lon=float(columns[1][0:-2])
            sit=float(columns[2][0:-2])    # Thickness
            sigo=float(columns[3][0:-2])   # Thickness error
            LON=np.append(LON, lon)
            LAT=np.append(LAT, lat)
            VAR=np.append(VAR, sit)
            OBS_ERROR=np.append(OBS_ERROR, sigo)

        except:
            pass

    #plt.plot(LON, ocean_obs_utils.smooth(VAR,window_len=200),'-r')
    #plt.plot(LON, VAR,'--k',alpha=0.2)
    #plt.savefig('OIB.png')

    VAR = ocean_obs_utils.smooth(VAR,window_len=200)

    f.close()
    VAR = np.reshape(VAR, (len(VAR),1))
    yyyymmddhh=yyyy+mm+dd+hh
    DATE_TIME = int(yyyymmddhh)*np.ones(np.shape(LON))
    DEPTH = 15.0*np.ones(np.shape(VAR))
    N_LEVS=1
    QC_LEV = 1.0*np.ones(np.shape(VAR))
    QC_PRF = 1.0*np.ones(np.shape(LON))
    INST_ID = 1.0*np.ones(np.shape(INST_ID))

    print(np.shape(LON), np.shape(VAR))

    return N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID

class Obs:
    """
    """
    def __init__(self, yyyy, mm, dd, hh, fname, id_obs, vartype='TEMP', xsigo=1.0, descriptor='', platform = '', color='r', markersize=5, NDAYS=EXP_NDAYS):
        """ This function exctract the observations that are within the assimilation window
        Args:
            yyyy (int)     : 4 digit year
              mm (int)     : 2 digit month
              dd (int)     : 2 digit day
              hh (int)     : 2 digit hour
              fname        : Name of file containing the obs
              id_obs       : Obs identifier
              vartype      : Variable type, one of ['TEMP','SALT','ADT','ICE'] Note that ICE is ice fraction
              xsigo        : Used to rescale the obserror. new_sigo = xsigo x old_sigo
              descriptor   : String describing the instance of the object (ex: T-Argo, Jason-1-ADT, ...)
              NDAYS (float): Used to define the size of the assimilation window centered at yyyymmddhh [yyyymmddhh-NDAYS, yyyymmddhh+NDAYS].
    """

        print('platform is', platform)

        if ( (platform == 'NOAA16') | (platform == 'METOPA') ):
            print('L2-SST')
            try:
                N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR = l2_sst_reader(yyyy, mm, dd, hh, platform=platform, NDAYS=EXP_NDAYS)
            except:
                self.no_obs(descriptor = descriptor, platform = platform, color = color, size = markersize, present=False)
                print('Failed looking up ',self.descriptor)
                return
        elif (platform == 'GMI'):
            #try:
            N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR = l2_gmi_reader(yyyy, mm, dd, hh, NDAYS=EXP_NDAYS)
            #except:
            #    self.no_obs(descriptor = descriptor, platform = platform, color = color, size = markersize, present=False)
            #    print 'Failed looking up ',self.descriptor
            #    return
        elif (platform == 'AIR-BORN'):
            try:
                N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR = oib_reader(yyyy, mm, dd, hh, NDAYS=EXP_NDAYS)
            except:
                self.no_obs(descriptor = descriptor, platform = platform, color = color, size = markersize, present=False)
                print('Failed looking up ',self.descriptor)
                return
        elif (platform == 'SMOSSUB'):
            try:
                print('GOT HERE SMOSSUB')
                N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID = l2_smossub_reader(fname, vartype, platform=platform) # Reads smos observation format
            except:
                self.no_obs(descriptor = descriptor, platform = platform, color = color, size = markersize, present=False)
                print('Failed looking up ',self.descriptor)
                return
        elif (platform == 'CS2-HICE'):
            try:
                N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR = cs2_reader(yyyy, mm, dd, hh, NDAYS=EXP_NDAYS)
            except:
                self.no_obs(descriptor = descriptor, platform = platform, color = color, size = markersize, present=False)
                print('Failed looking up ',self.descriptor)
                return
        elif (platform == 'M2-SST'):
            print ('in M2-SST')
            N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID = M2_sst_reader(yyyy, mm, dd, hh, path2scratch=SCRDIR, path2expdir=EXPDIR, path2expid=EXPID)
        else:
            if not(os.path.exists(fname)):
                self.no_obs(descriptor = descriptor, platform = platform, color = color, size = markersize, present=False)
                #self.typ   = []
                #self.lon   = []
                #self.lat   = []
                #self.depth = []
                #self.value = []
                #self.oerr  = []
                #self.descriptor = descriptor
                #self.platform = platform
                #self.color=color
                #self.size=markersize
                #self.present=False
                return
            N_LEVS, DEPTH, VAR, QC_LEV, QC_PRF, LON, LAT, DATE_TIME, OBS_ERROR, INST_ID = standard_obs_reader(fname, vartype) # Reads standard iodas observation format

        print(':::::::::::::::::::',descriptor)
        I=da_window(yyyy, mm, dd, hh, NDAYS, DATE_TIME, QC_PRF)
        if (descriptor=='AVHRR18 L2 SST'):
            print('L2 SST')
            I=da_window(yyyy, mm, dd, hh, NDAYS, np.floor_divide(DATE_TIME,10000), QC_PRF)

        if (len(I)==0): # File present but no obs within window
            self.no_obs(descriptor = descriptor, platform = platform, color = color, size = markersize, present=False)
            return

#   here for all insitu data
        DATE_TIME=DATE_TIME[I]
        LON=LON[I]
        LAT=LAT[I]
        VAR=VAR[I,0:N_LEVS]
        DEPTH=DEPTH[I,0:N_LEVS]
        QC_LEV=QC_LEV[I,0:N_LEVS]
        OBS_ERROR=OBS_ERROR[I,0:N_LEVS]
        INST_ID=INST_ID[I]

        lon=np.zeros(np.shape(VAR))
        lat=np.zeros(np.shape(VAR))
        instid=np.zeros(np.shape(VAR))
        date_time=np.zeros(np.shape(VAR))
        for prof in range(np.shape(VAR)[0]):
            lon[prof,:]=LON[prof]
            lat[prof,:]=LAT[prof]
            instid[prof,:]=INST_ID[prof]
            date_time[prof,:]=DATE_TIME[prof]
        lon=lon.flatten()
        lat=lat.flatten()
        depth=DEPTH.flatten()
        value=VAR.flatten()
        instid=instid.flatten()
        date_time=date_time.flatten()

        # Compute sigo's for insitu profiles
        print('id_obs is',id_obs)
        if ( (id_obs==obsid_dict['id_t_obs']) | (id_obs==obsid_dict['id_s_obs']) ): # Profiles
#           print ('compute profile sigos ',np.shape(VAR),range(np.shape(VAR)[0]),range(N_LEVS-1))
            dTdZ=np.zeros(np.shape(VAR))
            for indexp in range(np.shape(VAR)[0]):
#     first get the real depth of this profile
                nlevprof=0
                for indexz in range(np.shape(VAR)[1]):
                    if(DEPTH[indexp,indexz]<=100000.):   #Assume obs depth < 100000m
                        nlevprof=nlevprof+1

#                   print('Real depth of profile ',nlevprof)
                for indexz in range(0,nlevprof-1):
                    dTdZ[indexp, indexz] = np.abs(( VAR[indexp, indexz] - VAR[indexp, indexz+1] )/( DEPTH[indexp, indexz] - DEPTH[indexp, indexz+1] ))

                dTdZ[indexp, nlevprof-1]=dTdZ[indexp, nlevprof-2] # make the bottom gradient equal to one up

#                if (descriptor == 'CTD-S' or descriptor == 'CTD-T' ):
#   now eliminate CTD under the ice
#                       print 'SCRDIR IS', SCRDIR
#                    fname=SCRDIR+'/rawM2sic_'+yyyy+mm+dd+'.nc'
#                    AICE = readnc(fname,'AICE')
#                    AICE.x = readnc(fname,'lonout')
#                    AICE.y = readnc(fname,'latout')
#                       print AICE.y
#                       print np.shape(AICE)
#                find the index of the obs
#                    iax = np.shape(AICE.x)
#                    iay = np.shape(AICE.y)
#                       print iax,iay
#                    lonminstore=1000.
#                    latminstore=1000.
#                    for ilon in range(1440):
#                        print AICE.x[ilon], LON[indexp]
#                        lonmin = np.abs(AICE.x[ilon] - LON[indexp])
                        # print 'lonmin',lonmin,lonminstore,AICE.x[ilon],LON[indexp]
#                        if (lonmin < lonminstore):
#                            ilon1 = ilon
#                            lonminstore = lonmin
#                    for ilat in range(720):
#                        latmin = np.abs(AICE.y[ilat] - LAT[indexp])
#                        if (latmin < latminstore):
#                            ilat1 = ilat
#                            latminstore = latmin
##                       print 'found',ilon1,ilat1,LON[indexp],LAT[indexp],indexp
#                    if (AICE[ilat1,ilon1]>0.1):   # it's ice
#                            print 'eliminate ice CTD',ilon1,ilat1,LON[indexp],LAT[indexp],indexp,AICE[ilat1,ilon1]
#                        print(('eliminate ice CTD',LON[indexp],LAT[indexp],indexp,AICE[ilat1,ilon1]), descriptor)
#                            for ilev in range(N_LEVS):
#                                if(depth[ilev] < 60.):
#                                   QC_LEV[indexp,ilev]=0
#                                   dTdZ[indexp,ilev]=9999.
#                                   print 'removing data at',depth[ilev],ilev,indexp
#                        QC_LEV[indexp,0:N_LEVS]=0
#                        dTdZ[indexp,0:N_LEVS]=9999.
#   12/15/22 ERIC
#                if (descriptor == 'XBT-T'):
#   now eliminate XBT under the ice
#                       print 'SCRDIR IS', SCRDIR
#                    fname=SCRDIR+'/rawM2sic_'+yyyy+mm+dd+'.nc'
#                    AICE = readnc(fname,'AICE')
#                    AICE.x = readnc(fname,'lonout')
#                    AICE.y = readnc(fname,'latout')
#                       print AICE.y
#                       print np.shape(AICE)
#                find the index of the obs
#                    iax = np.shape(AICE.x)
#                    iay = np.shape(AICE.y)
#                       print iax,iay
#                    lonminstore=1000.
#                    latminstore=1000.
#                    for ilon in range(1440):
#                        print AICE.x[ilon], LON[indexp]
#                        lonmin = np.abs(AICE.x[ilon] - LON[indexp])
                        # print 'lonmin',lonmin,lonminstore,AICE.x[ilon],LON[indexp]
#                        if (lonmin < lonminstore):
#                            ilon1 = ilon
#                            lonminstore = lonmin
#                    for ilat in range(720):
#                        latmin = np.abs(AICE.y[ilat] - LAT[indexp])
#                        if (latmin < latminstore):
#                            ilat1 = ilat
#                            latminstore = latmin
#                       print 'found',ilon1,ilat1,LON[indexp],LAT[indexp],indexp
#                    if (AICE[ilat1,ilon1]>0.1):   # it's ice
#                        print(('eliminate ice XBT',LON[indexp],LAT[indexp],indexp,AICE[ilat1,ilon1]))
#                        QC_LEV[indexp,0:N_LEVS]=0
#                        dTdZ[indexp,0:N_LEVS]=9999.
#                if (descriptor == 'Argo-S' or descriptor == 'Argo-T' ):
#                    print ('IN ARGO CHECK FOR ICE')
#   now eliminate Argo under the ice
#                       print 'SCRDIR IS', SCRDIR
#                    fname=SCRDIR+'/rawM2sic_'+yyyy+mm+dd+'.nc'
#                    AICE = readnc(fname,'AICE')
#                    AICE.x = readnc(fname,'lonout')
#                    AICE.y = readnc(fname,'latout')
#                       print AICE.y
#                       print np.shape(AICE)
#                find the index of the obs
#                    iax = np.shape(AICE.x)
#                    iay = np.shape(AICE.y)
#                       print iax,iay
#                    lonminstore=1000.
#                    latminstore=1000.
#                    for ilon in range(1440):
#                        print AICE.x[ilon], LON[indexp]
#                        lonmin = np.abs(AICE.x[ilon] - LON[indexp])
                        # print 'lonmin',lonmin,lonminstore,AICE.x[ilon],LON[indexp]
#                        if (lonmin < lonminstore):
#                            ilon1 = ilon
#                            lonminstore = lonmin
#                    for ilat in range(720):
#                        latmin = np.abs(AICE.y[ilat] - LAT[indexp])
#                        if (latmin < latminstore):
#                            ilat1 = ilat
#                            latminstore = latmin
#                       print 'found',ilon1,ilat1,LON[indexp],LAT[indexp],indexp
#                    if (AICE[ilat1,ilon1]>0.1):   # it's ice
#                            print 'eliminate ice Argo',ilon1,ilat1,LON[indexp],LAT[indexp],indexp,AICE[ilat1,ilon1]
#                        print(('eliminate ice Argo',LON[indexp],LAT[indexp],indexp,AICE[ilat1,ilon1]))
#                            for ilev in range(N_LEVS):
#                                if(depth[ilev] < 60.):
#                                   QC_LEV[indexp,ilev]=0
#                                   dTdZ[indexp,ilev]=9999.
#                                   print 'removing data at',depth[ilev],ilev,indexp
#                if (descriptor == 'WOA18_atArgo-T' or descriptor == 'WOA18_atArgo-S' or descriptor == 'WOA18_atNA-T' or descriptor == 'WOA18_atNA-S' ):
#                   print ('IN WOA18 CHECK FOR ICE')
#   now eliminate Argo under the ice
#                       print 'SCRDIR IS', SCRDIR
#                    fname=SCRDIR+'/rawM2sic_'+yyyy+mm+dd+'.nc'
#                    AICE = readnc(fname,'AICE')
#                    AICE.x = readnc(fname,'lonout')
#                    AICE.y = readnc(fname,'latout')
#                       print AICE.y
#                       print np.shape(AICE)
#                find the index of the obs
#                    iax = np.shape(AICE.x)
#                    iay = np.shape(AICE.y)
#                       print iax,iay
#                    lonminstore=1000.
#                    latminstore=1000.
#                    ilon1=999
#                    for ilon in range(1440):
#                        print AICE.x[ilon], LON[indexp]
#                        lonmin = np.abs(AICE.x[ilon] - LON[indexp])
#                         print 'lonmin',lonmin,lonminstore,AICE.x[ilon],LON[indexp],ilon1,ilon
#                        if (lonmin < lonminstore):
#                            ilon1 = ilon
#                            lonminstore = lonmin
#                   print(ilon1)
#                    ilat1=999
#                    for ilat in range(720):
#                        latmin = np.abs(AICE.y[ilat] - LAT[indexp])
#                       print 'latmin',latmin,latminstore,AICE.y[ilat],LAT[indexp],ilat1,ilat
#                        if (latmin < latminstore):
#                            ilat1 = ilat
#                            latminstore = latmin
#                         print 'found',ilon1,ilat1,LON[indexp],LAT[indexp],indexp,AICE[ilat1,ilon1]
#                   exit()
#                   if(LAT[indexp]>70.):
#                       print ('xxxxxxxxxxxxxx',ilon1,ilat1,LON[indexp],LAT[indexp],indexp,AICE[ilat1,ilon1])
             
#                    if (AICE[ilat1,ilon1]>0.1):   # it's ice
#                        print('eliminate ice WOA',ilon1,ilat1,LON[indexp],LAT[indexp],indexp,AICE[ilat1,ilon1])
#                       print('eliminate ice WOA',LON[indexp],LAT[indexp],indexp,AICE[ilat1,ilon1])
#                            for ilev in range(N_LEVS):
#                                if(depth[ilev] < 60.):
#                                   QC_LEV[indexp,ilev]=0
#                                   dTdZ[indexp,ilev]=9999.
#                                   print 'removing data at',depth[ilev],ilev,indexp
#                        QC_LEV[indexp,0:N_LEVS]=0
#                        dTdZ[indexp,0:N_LEVS]=9999.
#                        QC_LEV[indexp,0:N_LEVS]=0
#                        dTdZ[indexp,0:N_LEVS]=9999.
            dTdZ[dTdZ<1e-3]=1e-3
            OBS_ERROR[:]=2.0*dTdZ[:]

        # Make sure bad obs are out
        I=np.where(QC_LEV.flatten()==1)

        date_time=date_time[I]
        lon=lon[I]
        lat=lat[I]
        depth=depth[I]
        value=value[I]
        instid=instid[I]
        typ=id_obs*np.ones(np.shape(value))
        sigo=np.squeeze(OBS_ERROR.flatten()[I])#[I])

#   ERIC31021 now get rid of obviously bad data
        I=np.where(np.abs(value)<100.)
        date_time=date_time[I]
        lon=lon[I]
        lat=lat[I]
        depth=depth[I]
        value=value[I]
 #      obs_error=obs_error[I]
        typ=id_obs*np.ones(np.shape(value))
        sigo=np.squeeze(sigo.flatten()[I])#[I])
        instid=np.squeeze(instid.flatten()[I])#[I])
#   ERIC102722 now get rid of obviously bad data
        I=np.where(np.abs(sigo)<100.)
        date_time=date_time[I]
        lon=lon[I]
        lat=lat[I]
        depth=depth[I]
        value=value[I]
 #      obs_error=obs_error[I]
        typ=id_obs*np.ones(np.shape(value))
        sigo=np.squeeze(sigo.flatten()[I])#[I])
        instid=np.squeeze(instid.flatten()[I])

#       print 'here',np.shape(sigo), np.shape(value)
        #raw_input('<>?')

#   7/28/22 added to eliminate Tz reporting as Argo Sz 
        if (id_obs==obsid_dict['id_s_obs']):
                I=(np.where((value >= 29) & (value <= 40)))
                date_time=date_time[I]
                lon=lon[I]
                lat=lat[I]
                depth=depth[I]
                value=value[I]
                typ=id_obs*np.ones(np.shape(value))
                sigo=np.squeeze(sigo.flatten()[I])#[I])
                instid=np.squeeze(instid.flatten()[I])#[I])

#   12/14/22 added to eliminate Tz reporting < -2C 
#       if (id_obs==obsid_dict['id_s_obs']):
#               I=(np.where(value < -2.))
#               lon=lon[I]
#               lat=lat[I]
#               depth=depth[I]
#               value=value[I]
#               typ=id_obs*np.ones(np.shape(value))
#               sigo=np.squeeze(sigo.flatten()[I])#[I])
#               instid=np.squeeze(instid.flatten()[I])#[I])

#   Yehui: 01/10/2023 
#   12/14/22 added to eliminate Tz reporting < -2C
        if (id_obs==obsid_dict['id_t_obs']):
                I=(np.where(value > -2.))
                date_time=date_time[I]
                lon=lon[I]
                lat=lat[I]
                depth=depth[I]
                value=value[I]
                typ=id_obs*np.ones(np.shape(value))
                sigo=np.squeeze(sigo.flatten()[I])#[I])
                instid=np.squeeze(instid.flatten()[I])#[I])
 
        if (id_obs==obsid_dict['id_sst_obs']):    #'SST'
            oerr_min=1e-2
            sigo[sigo<oerr_min]=oerr_min
            sigo[np.isnan(sigo)]=999.9

        if (id_obs==obsid_dict['id_sss_obs']):    #'SSS'
#           print ('hi e', descriptor)

#    double the SMOS error to see if that tames the increment some 2/1/19
#           smosrampuperr=2.
            if (descriptor=='SMOS_L2_SSS'):
                OBS_ERROR=OBS_ERROR*smosrampuperr

#   mask out where sss > 40
                I=(np.where((value >= 29) & (value <= 40)))
#              I=np.where(value<40)
                date_time=date_time[I]
                lon=lon[I]
                lat=lat[I]
                depth=depth[I]
                value=value[I]
#              obs_error=obs_error[I]
                typ=id_obs*np.ones(np.shape(value))
                sigo=np.squeeze(sigo.flatten()[I])#[I])
                instid=np.squeeze(instid.flatten()[I])#[I])

            if (descriptor=='SMAP_L2_SSS' or descriptor=='SMOS_L2_SSS' or descriptor=='Aquarius_L2_SSS'):
                oerr_min=0.1   # 0.1 psu
                sigo=np.squeeze(OBS_ERROR.flatten()[I])#[I])

            if (descriptor=='SMAP_L3_SSS' or descriptor=='SMOS_L3_SSS' or descriptor=='Aquarius_L3_SSS'):
                oerr_min=0.1   # 0.1 psu
                sigo[sigo<oerr_min]=oerr_min
                sigo[np.isnan(sigo)]=999.9
#              print ('here eric',sigo[0:30])

#              Impose large sigos in the high latitudes from https://aquarius.umaine.edu/docs/aqsci2015_meissner_2.pdf page 17 error map (from v4.0 data) linear from 30-60 max=0.5
                sigo_scale=np.zeros(np.shape(sigo))
                sigo_scale=0.4*((abs(lat)-30.)/30.0)
                sigo_scale[sigo_scale>0.4]=0.4
                sigo_scale[sigo_scale<0.0]=0.0
                sigo=(sigo_scale+sigo) #*sigo
                sigo[sigo>999.9]=999.9  #*sigo
#  hack to stop assim into arctic
                sigo[abs(lat)>50.]=999.9

            # Get rid of SSS obs that have sigo>999.9
                I=np.where(sigo<999.9)
                date_time=date_time[I]
                lon=lon[I]
                lat=lat[I]
                depth=depth[I]
                value=value[I]
#              obs_error=obs_error[I]
                typ=id_obs*np.ones(np.shape(value))
                sigo=np.squeeze(sigo.flatten()[I])#[I])
                instid=np.squeeze(instid.flatten()[I])#[I])

#   ADDED BY ERIC TO REMOVE ADT WITH NAN VALUES
        if id_obs==obsid_dict['id_eta_obs']: #'ADT'
#           I=np.where(np.isnan(value))
            I=np.where(np.abs(value)<10.)
            date_time=date_time[I]
            lon=lon[I]
            lat=lat[I]
            depth=depth[I]
            value=value[I]
#           obs_error=obs_error[I]
            typ=id_obs*np.ones(np.shape(value))
            sigo=np.squeeze(sigo.flatten()[I])#[I])
            instid=np.squeeze(instid.flatten()[I])#[I])

        if id_obs==obsid_dict['id_eta_obs']: #'ADT'
            sigo[np.isnan(sigo)]=999.9
            oerr_min=0.1 # 10 cm
            sigo[sigo<oerr_min]=oerr_min

            #Impose large sigos in the high latitudes
            sigo_scale=np.zeros(np.shape(sigo))
#           sigo_scale=0.01*(LAT[I]/90.0)**4
            sigo_scale=0.1*(lat/90.0)**4
            sigo=(sigo_scale+sigo) #*sigo
            sigo[sigo>999.9]=999.9  #*sigo
            #oerr[LAT<-60]=99.9

#   added by eric to test crop south of 55S
        if id_obs==obsid_dict['id_eta_obs']: #'ADT'
#           I=np.where(lat>-45.)
#           I=np.where(lat>-55.)
            I=np.where(lat>-50.)
            date_time=date_time[I]
            lon=lon[I]
            lat=lat[I]
            depth=depth[I]
            value=value[I]
            typ=typ[I]
            sigo=sigo[I]
            instid=instid[I]

#   eric test standard deviation of ADT observation error
        if id_obs==obsid_dict['id_eta_obs']: #'ADT'
            I=np.where(sigo!=999.9)
            date_time=date_time[I]
            lon=lon[I]
            lat=lat[I]
            depth=depth[I]
            value=value[I]
            typ=typ[I]
            sigo=sigo[I]
            instid=instid[I]
            sigo_sigma=np.std(sigo)
            if (sigo_sigma < .05):
                print ("BAD ADT SIGO FULL STOP ANALYSIS CYCLE")
                print ('PARENT ID IS', os.getppid(), os.getpid())
#      create a file to tell parents that obs are bad
                command='touch '+SCRDIR+'/BADOBS'
                os.system(command)
                sys.exit(1)

        '''
        if ( (id_obs==obsid_dict['id_t_obs']) | (id_obs==obsid_dict['id_s_obs']) ): # Profiles
            # Assumes a exponential decay of obs error
            if (id_obs==obsid_dict['id_t_obs']):
                prof_sigo = T_prof_sigo
            else:
                prof_sigo = S_prof_sigo
            #D0=200.0 #e-folding scale
            D0=2000000.0 #e-folding scale

            #print(np.shape(value))
            #raw_input('<>?')

            sigo=prof_sigo*np.ones(np.shape(value))*np.exp(-depth/D0)
            oerr_max=np.max(sigo)
            sigo[sigo<0.01*oerr_max]=0.01*oerr_max
        '''

        if ( (id_obs==obsid_dict['id_sss_obs']) ): # satellite SSS
            sigo[sigo<=0.0]=0.01  # don't allow any 0 error

        if ( (id_obs==obsid_dict['id_t_obs']) | (id_obs==obsid_dict['id_s_obs']) ): # Profiles
            sigo[sigo<=0.0]=0.01  # don't allow any 0 error for profile data

        if (id_obs==obsid_dict['id_aice_obs']):
            sigo = 0.05*np.ones(np.shape(value)) #value*0.05

        if (id_obs==obsid_dict['id_hice_obs']):
            sigo = 0.3*np.ones(np.shape(value)) #value*0.05
            #sigo[sigo<0.01]=0.01

        self.date_time = date_time
        self.typ   = typ
        self.lon   = lon
        self.lat   = lat
        self.depth = depth
        self.value = value
        self.instid = instid
        self.oerr  = xsigo*sigo       # Rescaled sigos
        self.descriptor = descriptor
        self.platform = platform
        self.color=color
        self.size=markersize
        if (len(self.typ)>0):
            self.present=True
        else:
            self.present=False



    def no_obs(self, descriptor='', platform='', color='r',size=10,present=False):
        self.date_time=[]
        self.typ   = []
        self.lon   = []
        self.lat   = []
        self.depth = []
        self.value = []
        self.oerr  = []
        self.instid= []
        self.descriptor = descriptor
        self.platform = platform
        self.color=color
        self.size=size
        self.present=False

    def transform(self, transtyp='logit'):
        '''
        Transform the obs and sigo, currently only supports the logit transform ...
        '''
        if self.present:
            if transtyp=='logit':
            #Deal with 0% and 100%
                tiny=1e-3
                self.value[self.value<tiny]=tiny
                self.value[self.value>1.0-tiny]=1.0-tiny
                self.value=np.log(self.value) - np.log(1 - self.value)
            #print(np.min(self.value[np.isfinite(self.value)].flatten()))
            #print(np.max(self.value[np.isfinite(self.value)].flatten()))
            #self.oerr= ... sigos need to be transormed as well

            if transtyp=='invlogit':
                self.value = np.exp(self.value) / (1 + np.exp(self.value))

    def plot(self, pngname):
        if self.present:
#            fig = plt.figure(num=1, figsize=(10,8), facecolor='w')
#            fig.add_subplot(111)
            x, y = self.lon, self.lat
            fig = plt.figure(figsize=(10,8))

#            map = Basemap(projection='moll', llcrnrlat=-90, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c', lon_0=-80)
#            x, y = list(map(self.lon, self.lat))
            if ( (self.typ[0]==5351) | (self.typ[0]==5525) | (self.typ[0]==5522) ):
                if (self.typ[0]==5351):
                    valmin=-1.5
                    valmax=1.5
                    errmax=0.25
                if (self.typ[0]==5525):
                    print('plotting SST ....')
                    print('sst max:',np.min(self.lat))
                    valmin=-2.0
                    valmax=31.0
                    errmax=1.0
                if (self.typ[0]==5522):
                    print('plotting SSS ....')
                    print('sss max:',np.min(self.lat))
                    valmin=30.0
                    valmax=38.0
                    if (self.descriptor=='SMAP_L2_SSS'):
                        efact=2.
                    if (self.descriptor=='Aquarius_L2_SSS'):
                        efact=0.5
                    if (self.descriptor=='SMOS_L2_SSS'):
                        efact=2.0*smosrampuperr
                    if (self.descriptor=='SMOS_L3_SSS' or self.descriptor=='SMAP_L3_SSS' or self.descriptor=='Aquarius_L3_SSS'):
                        efact=0.5
                    #errmax=0.5*xsigo_sss
                    errmax=efact*xsigo_sss

                ax = plt.subplot(2,1,1, projection=ccrs.Mollweide(central_longitude=-80))
                ax.set_global()
                ax.coastlines()
                ax.add_feature(cfeature.BORDERS, lw=.5)
                ax.add_feature(cfeature.RIVERS)
                c= ax.scatter(x, y, s=1, c=self.value,transform=ccrs.PlateCarree(),
                                  cmap=cm.jet,vmin=valmin,vmax=valmax,edgecolor=None,lw=0)
                fig.colorbar(c, shrink=0.5, ax = ax)
 
                ax2 = plt.subplot(2,1,2, projection=ccrs.Mollweide(central_longitude=-80))
                ax2.set_global()
                ax2.coastlines()
                ax2.add_feature(cfeature.BORDERS, lw=.5)
                ax2.add_feature(cfeature.RIVERS)
                ax2.add_feature(cfeature.LAND, facecolor=("coral"))
                c= ax2.scatter(x, y, s=1, c=self.oerr,transform=ccrs.PlateCarree(),
                                   cmap=cm.jet,vmin=0,vmax=errmax,edgecolor=None,lw=0)
                fig.colorbar(c, shrink=0.5, ax=ax2)


            elif ( (self.typ[0]==6000) | (self.typ[0]==6001)):
                valmin=0.
                if (self.typ[0]==6000):
                    valmax=1.0
                if (self.typ[0]==6001):
                    valmax=4.0
                errmax=0.5

                ax = plt.subplot(2,1,1, projection=ccrs.NorthPolarStereo())
                ax.coastlines(resolution='110m')
                ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())
                ax.add_feature(cfeature.BORDERS, lw=.5)
                ax.add_feature(cfeature.RIVERS)
                ax.add_feature(cfeature.LAND, facecolor=("coral"))
                if logit_transform:
                    c = ax.scatter(x, y, s=1, c=inv_logit(self.value), transform=ccrs.PlateCarree(),
                                cmap=cm.jet,vmin=valmin,vmax=valmax,edgecolor=None,lw=0)
                else:
                    c = ax.scatter(x, y, 1, c=self.value, transform=ccrs.PlateCarree(),
                                cmap=cm.jet,vmin=valmin,vmax=valmax,edgecolor=None,lw=0)

                fig.colorbar(c, shrink=0.5, ax = ax)

                ax2 = plt.subplot(2,1,2, projection=ccrs.SouthPolarStereo())
                ax2.coastlines(resolution='110m')
                ax2.set_extent([-180, 180, -90, -55], crs=ccrs.PlateCarree())
                ax2.add_feature(cfeature.BORDERS, lw=.5)
                ax2.add_feature(cfeature.RIVERS)
                ax2.add_feature(cfeature.LAND, facecolor=("coral"))
                if logit_transform:
                    c = ax2.scatter(x, y, 1, c=inv_logit(self.value), transform=ccrs.PlateCarree(),
                                    cmap=cm.jet,vmin=valmin,vmax=valmax,edgecolor=None,lw=0)
                else:
                    c = ax2.scatter(x, y, 1, c=self.value, transform=ccrs.PlateCarree(),
                                    cmap=cm.jet,vmin=valmin,vmax=valmax,edgecolor=None,lw=0)
                fig.colorbar(c, shrink=0.5, ax = ax2)

            else:
                ax = fig.add_subplot(111, projection=ccrs.Mollweide(central_longitude=-80))
                ax.set_global()
                ax.coastlines()
                ax.add_feature(cfeature.BORDERS, lw=.5)
                ax.add_feature(cfeature.RIVERS)
                ax.add_feature(cfeature.LAND, facecolor=("coral"))
                ax.plot(x, y, color=self.color, transform=ccrs.PlateCarree(), marker='.', markersize= self.size, linestyle='None',alpha=0.2)

            titlestr=str(len(self.value))+' Obs'
            ax.set_title(titlestr)
            fig.savefig(self.descriptor+pngname)
            plt.clf()



def update_list_of_obs(list_of_obs, obs):
    if obs.present:
        list_of_obs.append(obs)

# Profiling drifters, moorings, ...
#==================================
def argo(list_of_obs):
#    argoName = '/gpfsm/dnb33/lren1/pre_proc/NRT'
    argo_t  = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/ARGO/V3/FINAL/T_ARGO_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', markersize=5, descriptor='Argo-T', NDAYS=EXP_NDAYS)
    argo_s  = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/ARGO/V3/FINAL/S_ARGO_'+yyyy+'.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='Argo-S', NDAYS=EXP_NDAYS)
#    argo_t  = Obs(yyyy, mm, dd, hh, fname=argoName+'/ARGO/V3/FINAL/T_ARGO_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', markersize=5, descriptor='Argo-T', NDAYS=EXP_NDAYS)
#    argo_s  = Obs(yyyy, mm, dd, hh, fname=argoName+'/ARGO/V3/FINAL/S_ARGO_'+yyyy+'.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='Argo-S', NDAYS=EXP_NDAYS)

    update_list_of_obs(list_of_obs, argo_t)
    update_list_of_obs(list_of_obs, argo_s)
    return list_of_obs

#   5 Randomly generated subsets of Argo 
def argo_1(list_of_obs):  
    argo_t_1  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/T_ARGO_'+yyyy+'_1.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', markersize=5, descriptor='Argo_1-T', NDAYS=EXP_NDAYS)
    argo_s_1  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/S_ARGO_'+yyyy+'_1.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='Argo_1-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, argo_t_1)
    update_list_of_obs(list_of_obs, argo_s_1)
    return list_of_obs

def argo_2(list_of_obs):  
    argo_t_2  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/T_ARGO_'+yyyy+'_2.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', markersize=5, descriptor='Argo_2-T', NDAYS=EXP_NDAYS)
    argo_s_2  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/S_ARGO_'+yyyy+'_2.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='Argo_2-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, argo_t_2)
    update_list_of_obs(list_of_obs, argo_s_2)
    return list_of_obs

def argo_3(list_of_obs):  
    argo_t_3  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/T_ARGO_'+yyyy+'_3.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', markersize=5, descriptor='Argo_3-T', NDAYS=EXP_NDAYS)
    argo_s_3  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/S_ARGO_'+yyyy+'_3.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='Argo_3-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, argo_t_3)
    update_list_of_obs(list_of_obs, argo_s_3)
    return list_of_obs

def argo_4(list_of_obs):  
    argo_t_4  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/T_ARGO_'+yyyy+'_4.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', markersize=5, descriptor='Argo_4-T', NDAYS=EXP_NDAYS)
    argo_s_4  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/S_ARGO_'+yyyy+'_4.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='Argo_4-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, argo_t_4)
    update_list_of_obs(list_of_obs, argo_s_4)
    return list_of_obs

def argo_5(list_of_obs):  
    argo_t_5  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/T_ARGO_'+yyyy+'_5.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', markersize=5, descriptor='Argo_5-T', NDAYS=EXP_NDAYS)
    argo_s_5  = Obs(yyyy, mm, dd, hh, fname=ERICLEV50INSITUOBSDIR+'/ARGO/V3/FINAL/S_ARGO_'+yyyy+'_5.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='Argo_5-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, argo_t_5)
    update_list_of_obs(list_of_obs, argo_s_5)
    return list_of_obs

def ctd(list_of_obs):
#    ctd_t   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/CTD/V3.1/T_CTD_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', descriptor='CTD-T', NDAYS=EXP_NDAYS)
    ctd_t   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/CTD/V3/FINAL/T_CTD_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', descriptor='CTD-T', NDAYS=EXP_NDAYS)
#    ctd_s   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/CTD/V3.1/S_CTD_'+yyyy+'.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='CTD-S', NDAYS=EXP_NDAYS)
    ctd_s   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/CTD/V3/FINAL/S_CTD_'+yyyy+'.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='CTD-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, ctd_t)
    update_list_of_obs(list_of_obs, ctd_s)
    return list_of_obs

def xbt(list_of_obs):
#    xbt_t   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/XBT/V3.1/T_XBT_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', descriptor='XBT-T', NDAYS=EXP_NDAYS)
    xbt_t   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/XBT/V3/FINAL/T_XBT_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', descriptor='XBT-T', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, xbt_t)
    return list_of_obs

def xbt_synS(list_of_obs):
    print(SYNOBSDIR+'/SYN_XBT_')
    xbt_s   = Obs(yyyy, mm, dd, hh, fname=SYNOBSDIR+'/XBT/SYN_XBT_'+yyyy+'.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='XBT-SYN-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, xbt_s)
    return list_of_obs

def woa18_atargo(list_of_obs):
    woa18_atargo_t  = Obs(yyyy, mm, dd, hh, fname=WOA18+'/WOA_ARGO/T_WOAatARGO_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', markersize=5, descriptor='WOA18_atArgo-T', NDAYS=EXP_NDAYS)
    woa18_atargo_s  = Obs(yyyy, mm, dd, hh, fname=WOA18+'/WOA_ARGO/S_WOAatARGO_'+yyyy+'.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='WOA18_atArgo-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, woa18_atargo_t)
    update_list_of_obs(list_of_obs, woa18_atargo_s)
    return list_of_obs

def woa18_atna(list_of_obs):
    woa18_atna_t= Obs(yyyy, mm, dd, hh, fname=WOA18+'/WOA_NA/T_WOAatNA_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='c', markersize=5, descriptor='WOA18_atNA-T', NDAYS=EXP_NDAYS)
    woa18_atna_s= Obs(yyyy, mm, dd, hh, fname=WOA18+'/WOA_NA/S_WOAatNA_'+yyyy+'.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='c', descriptor='WOA18_atNA-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, woa18_atna_t)
    update_list_of_obs(list_of_obs, woa18_atna_s)
    return list_of_obs

def tao(list_of_obs):
#    tao_t   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/TAO/V3.1/T_TAO_'+yyyy+'.nc',   id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='r', descriptor='TAO-T', NDAYS=EXP_NDAYS)
    tao_t   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/TAO/V3/FINAL/T_TAO_'+yyyy+'.nc',   id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='r', descriptor='TAO-T', NDAYS=EXP_NDAYS)
#    tao_s   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/TAO/V3.1/S_TAO_'+yyyy+'.nc',   id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='r', descriptor='TAO-S', NDAYS=EXP_NDAYS)
    tao_s   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/TAO/V3/FINAL/S_TAO_'+yyyy+'.nc',   id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=xsigo_s, color='r', descriptor='TAO-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, tao_t)
    update_list_of_obs(list_of_obs, tao_s)
    return list_of_obs

def pirata(list_of_obs):
#    pir_t   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/PIRATA/V3.1/T_PIR_'+yyyy+'.nc',   id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='r', descriptor='PIRATA-T', NDAYS=EXP_NDAYS)
    pir_t   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/PIRATA/V3/FINAL/T_PIR_'+yyyy+'.nc',   id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='r', descriptor='PIRATA-T', NDAYS=EXP_NDAYS)
#    pir_s   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/PIRATA/V3.1/S_PIR_'+yyyy+'.nc',   id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=10.0*xsigo_s, color='r', descriptor='PIRATA-S', NDAYS=EXP_NDAYS)
    pir_s   = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/PIRATA/V3/FINAL/S_PIR_'+yyyy+'.nc',   id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=10.0*xsigo_s, color='r', descriptor='PIRATA-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, pir_t)
    update_list_of_obs(list_of_obs, pir_s)
    return list_of_obs

def rama(list_of_obs):
#   rama_t  = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/RAMA/V3.1/T_RAMA_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='r', descriptor='RAMA-T', NDAYS=EXP_NDAYS)
    rama_t  = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/RAMA/V3/FINAL/T_RAMA_'+yyyy+'.nc', id_obs=obsid_dict['id_t_obs'], vartype='TEMP', xsigo=xsigo_t, color='r', descriptor='RAMA-T', NDAYS=EXP_NDAYS)
#   rama_s  = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/RAMA/V3.1/S_RAMA_'+yyyy+'.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=10.0*xsigo_s, color='r', descriptor='RAMA-S', NDAYS=EXP_NDAYS)
    rama_s  = Obs(yyyy, mm, dd, hh, fname=LEV50INSITUOBSDIR+'/RAMA/V3/FINAL/S_RAMA_'+yyyy+'.nc', id_obs=obsid_dict['id_s_obs'], vartype='SALT', xsigo=10.0*xsigo_s, color='r', descriptor='RAMA-S', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, rama_t)
    update_list_of_obs(list_of_obs, rama_s)
    return list_of_obs

# SSS satellite
#===================
def aq_L2_sss(list_of_obs):     #Aquarius V5.0
    fname=SatSSSOBSDIR+'L2_AQ_SSS_7.0/V5/'+'SSS_TRK_AQ_V5_'+yyyy+'.nc'
    print('fname is',fname)
#   aq_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='Aquarius_L3_SSS', color='y')
    aq_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='Aquarius_L2_SSS', color='y')
    update_list_of_obs(list_of_obs, aq_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

def smap_L2_sss(list_of_obs):     #SMAP V4.0
    fname=SatSSSOBSDIR+'L2_SMAP_SSS_7.0/V4/'+'SSS_TRK_SMAP_V4_'+yyyy+'.nc'
    print('fname is',fname)
    smap_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMAP_L2_SSS', color='y')
    update_list_of_obs(list_of_obs, smap_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

def smap_v4_1_L2_sss(list_of_obs):     #SMAP V4.1  currently only for 6/18-9/18
    fname=SatSSSOBSDIR+'L2_SMAP_SSS_7.0/V4.1/'+'SSS_TRK_SMAP_V4.1_'+yyyy+'.nc'
    print('fname is',fname)
    smap_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMAP_L2_SSS', color='y')
    update_list_of_obs(list_of_obs, smap_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

def smap_v4_2_L2_sss(list_of_obs):     #SMAP V4.2  currently only for 2019
    fname=SatSSSOBSDIR+'L2_SMAP_SSS_7.0/V4.2/'+'SSS_TRK_SMAP_V4.2_'+yyyy+'.nc'
    print('fname is correct',fname)
    smap_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMAP_L2_SSS', color='y')
    update_list_of_obs(list_of_obs, smap_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

def smap_v4_3_L2_sss(list_of_obs):     #SMAP V4.3
    fname=SatSSSOBSDIR+'L2_SMAP_SSS_7.0/V4.3/'+'SSS_TRK_SMAP_V4.3_'+yyyy+'.nc'
    print('fname is correct',fname)
    smap_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMAP_L2_SSS', color='y')
    update_list_of_obs(list_of_obs, smap_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

def smap_v5_0_L2_sss(list_of_obs):     #SMAP V5.0
    fname=SatSSSOBSDIR+'L2_SMAP_SSS_7.0/V5.0/'+'SSS_TRK_SMAP_V5.0_'+yyyy+'.nc'
    print('fname is correct',fname)
    smap_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMAP_L2_SSS', color='y')
    update_list_of_obs(list_of_obs, smap_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

def smos_L2_sss(list_of_obs):     #SMOS V3.0
    fname=SatSSSOBSDIR+'L2_SMOS_SSS_7.0/L32Q/'+'SSS_TRK_SMOS_L32Q_'+yyyy+'.nc'
    print('fname is',fname,hh)
    smos_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMOS_L2_SSS', color='y')
#  only grab the 12z data
#   if (hh=='12'):
#     smos_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMOS_L2_SSS', color='y')

    update_list_of_obs(list_of_obs, smos_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

def smossub_L2_sss(list_of_obs):     #SMOS V3.0
    fname=SatSSSOBSDIR+'L2_SMOS_SSS_7.0/L32Q/'+'SSS_TRK_SMOS_L32Q_'+yyyy+'.nc'
    print('fname is',fname,hh)
    smos_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMOS_L2_SSS', color='y',platform='SMOSSUB')
#  only grab the 12z data
#   if (hh=='12'):
#     smos_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMOS_L2_SSS', color='y')

    update_list_of_obs(list_of_obs, smos_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

def aq_L3_sss(list_of_obs):     #Aquarius V5.0
#  ERIC AS IT NOW STANDS IT ONLY LIKES 12z
    fname=SatSSSOBSDIR+'L3_AQ_SSS_7.0/7DAYRUN/'+'SSS_GRD_AQ_'+yyyy+'.nc'
    print('fname is',fname)
    hh12=12  # pick up 12z all day
#   aq_sss  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='Aquarius_L3_SSS', color='y')
    aq_sss  = Obs(yyyy, mm, dd, hh12, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='Aquarius_L3_SSS', color='y')
    update_list_of_obs(list_of_obs, aq_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

#def smap_L3_sss(list_of_obs):     #SMAP V2.0
##  ERIC AS IT NOW STANDS IT ONLY LIKES 12z
#    fname=SatSSSOBSDIR+'L3_SMAP_SSS_7.0/8DAYRUN/'+'SSS_GRD_SMAP_'+yyyy+'.nc'
def smap_L3_sss(list_of_obs):     #SMAP V4.0
#  ERIC AS IT NOW STANDS IT ONLY LIKES 12z
    fname=SatSSSOBSDIR+'L3_SMAP_SSS_7.0/8DAYRUN/'+'SSS_GRD_SMAP_V4_'+yyyy+'.nc'
    print('SMAP L3 fname is',fname)
    hh12=12  # pick up 12z all day
    smap_sss  = Obs(yyyy, mm, dd, hh12, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMAP_L3_SSS', color='y')
    update_list_of_obs(list_of_obs, smap_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

def smos_L3_sss(list_of_obs):     #SMOS V3.0
#  ERIC AS IT NOW STANDS IT ONLY LIKES 12z
    fname=SatSSSOBSDIR+'L3_SMOS_SSS_7.0/'+'SSS_GRD_SMOS_'+yyyy+'.nc'
    print('SMOS fname is',fname)
    hh12=12  # pick up 12z all day
    smos_sss  = Obs(yyyy, mm, dd, hh12, fname=fname, id_obs=obsid_dict['id_sss_obs'], vartype='SSS',xsigo=xsigo_sss, descriptor='SMOS_L3_SSS', color='y')
    update_list_of_obs(list_of_obs, smos_sss)
#   print ('list_of_obs is ',list_of_obs)
    return list_of_obs

# Altimeters
#===================
# Li Ren added SWON
def swotn(list_of_obs):        #SWON
    fname=ADTOBSDIR+'ADT_TRK_SWON_'+yyyy+'.nc'
    sw_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='SWOTN-ADT', color='k')
    update_list_of_obs(list_of_obs, sw_adt)
    return list_of_obs

# Li Ren added S3b
def sentinel3b(list_of_obs):        #S3B-Sentinel
    fname=ADTOBSDIR+'ADT_TRK_S3B_'+yyyy+'.nc'
    sw_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='S3B-ADT', color='c')
    update_list_of_obs(list_of_obs, sw_adt)
    return list_of_obs

def cryosat2(list_of_obs):     #CryoSat-2
    fname=ADTOBSDIR+'ADT_TRK_C2_'+yyyy+'.nc'
    c2_adt  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='CryoSat-2-ADT', color='y')
    update_list_of_obs(list_of_obs, c2_adt)
    return list_of_obs

def cryosat2_N(list_of_obs):     #CryoSat-2-N
    fname=ADTOBSDIR+'ADT_TRK_C2N_'+yyyy+'.nc'
    c2n_adt  = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='CryoSat-2-N-ADT', color='y')
    update_list_of_obs(list_of_obs, c2n_adt)
    return list_of_obs

def sentinel3a(list_of_obs):        #S3A-Sentinel
    fname=ADTOBSDIR+'ADT_TRK_S3A_'+yyyy+'.nc'
    s3a_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Sentinel-3a-ADT', color='m')
    update_list_of_obs(list_of_obs, s3a_adt)
    return list_of_obs

def sentinel6a(list_of_obs):        #S6A-Sentinel
    fname=ADTOBSDIR+'ADT_TRK_S6A_'+yyyy+'.nc'
    s6a_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Sentinel-6a-ADT', color='m')
    update_list_of_obs(list_of_obs, s6a_adt)
    return list_of_obs

def jason1(list_of_obs):        #Jason-1
    fname=ADTOBSDIR+'ADT_TRK_J1_'+yyyy+'.nc'
    j1_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Jason-1-ADT', color='m')
    update_list_of_obs(list_of_obs, j1_adt)
    return list_of_obs

def jason1_N(list_of_obs):        #Jason-1
    fname=ADTOBSDIR+'ADT_TRK_J1N_'+yyyy+'.nc'
    j1n_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Jason-1-N-ADT', color='m')
    update_list_of_obs(list_of_obs, j1n_adt)
    return list_of_obs

def jason1G(list_of_obs):        #Jason-1
    fname=ADTOBSDIR+'ADT_TRK_J1G_'+yyyy+'.nc'
    j1g_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Jason-1G-ADT', color='m')
    update_list_of_obs(list_of_obs, j1g_adt)
    return list_of_obs

def jason2(list_of_obs):        #Jason-2
    fname=ADTOBSDIR+'ADT_TRK_J2_'+yyyy+'.nc'
    j2_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Jason-2-ADT', color='m')
    update_list_of_obs(list_of_obs, j2_adt)
    return list_of_obs

def jason2_N(list_of_obs):        #Jason-2-N
    fname=ADTOBSDIR+'ADT_TRK_J2N_'+yyyy+'.nc'
    j2n_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Jason-2-N-ADT', color='m')
    update_list_of_obs(list_of_obs, j2n_adt)
    return list_of_obs

def jason3(list_of_obs):        #Jason-3
    fname=ADTOBSDIR+'ADT_TRK_J3_'+yyyy+'.nc'
# ADT_TRK_J3_2017_sles12.nc
#   fname=ADTOBSDIR+'ADT_TRK_J3_'+yyyy+'_sles12_sigo.nc'
    j3_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Jason-3-ADT', color='m')
    update_list_of_obs(list_of_obs, j3_adt)
    return list_of_obs

def jason3N(list_of_obs):        #Jason-3N
    fname=ADTOBSDIR+'ADT_TRK_J3N_'+yyyy+'.nc'
    j3n_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Jason-3N-ADT', color='m')
    update_list_of_obs(list_of_obs, j3n_adt)
    return list_of_obs

def saral(list_of_obs):        #Saral/Altica
    fname=ADTOBSDIR+'ADT_TRK_AL_'+yyyy+'.nc'
#    fname=ADTOBSDIR+'ADT_TRK_AL_'+yyyy+'_sles12_sigo.nc'
    al_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Saral-Altika-ADT', color='m')
    update_list_of_obs(list_of_obs, al_adt)
    return list_of_obs

def ers1(list_of_obs):        #ERS-1
    fname=ADTOBSDIR+'ADT_TRK_E1_'+yyyy+'.nc'
    e1_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='ERS-1-ADT', color='m')
    update_list_of_obs(list_of_obs, e1_adt)
    return list_of_obs

def ers2(list_of_obs):        #ERS-2
    fname=ADTOBSDIR+'ADT_TRK_E2_'+yyyy+'.nc'
    e2_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='ERS-2-ADT', color='m')
    update_list_of_obs(list_of_obs, e2_adt)
    return list_of_obs

def topex_poseidon(list_of_obs):        #TOPEX/POSEIDON
    fname=ADTOBSDIR+'ADT_TRK_TP_'+yyyy+'.nc'
    tp_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='TOPEX-POSEIDON-ADT', color='m')
    update_list_of_obs(list_of_obs, tp_adt)
    return list_of_obs

def topex_poseidon_N(list_of_obs):        #TOPEX/POSEIDON
    fname=ADTOBSDIR+'ADT_TRK_TPN_'+yyyy+'.nc'
    tpn_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='TOPEX-POSEIDON-N-ADT', color='m')
    update_list_of_obs(list_of_obs, tpn_adt)
    return list_of_obs

def geosat_follow_on(list_of_obs):      #GEOSAT follow on
    fname=ADTOBSDIR+'ADT_TRK_G2_'+yyyy+'.nc'
    g2_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='GEOSAT-follow-on-ADT', color='m')
    update_list_of_obs(list_of_obs, g2_adt)
    return list_of_obs

def envisat(list_of_obs):       #Envisat
    fname=ADTOBSDIR+'ADT_TRK_EN_'+yyyy+'.nc'
    en_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Envisat-ADT', color='m')
    update_list_of_obs(list_of_obs, en_adt)
    return list_of_obs

def envisat_N(list_of_obs):       #Envisat-N
    fname=ADTOBSDIR+'ADT_TRK_ENN_'+yyyy+'.nc'
    enn_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='Envisat-N-ADT', color='m')
    update_list_of_obs(list_of_obs, enn_adt)
    return list_of_obs

def hy2a(list_of_obs):        #HY-2A
    fname=ADTOBSDIR+'ADT_TRK_H2_'+yyyy+'.nc'
    h2_adt = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_eta_obs'], vartype='ADT',xsigo=xsigo_ssh, descriptor='HY-2A-ADT', color='m')
    update_list_of_obs(list_of_obs, h2_adt)
    return list_of_obs

# SST retrieval
#===================
# Need to choose between Reyn, OSTIA and Hadley depending on date <<<<<<<======= Should be on the to do list
def reyn_L3_sst(list_of_obs):
    fname='/discover/nobackup/projects/gmao/ssd/g5odas/production/GEOS5odas-5.00/RC/OBS/SST_6.0/SST_REYN_'+yyyy+'.nc'
    reyn_t = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sst_obs'], vartype='TEMP',xsigo=xsigo_sst, descriptor='Reynolds-SST', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, reyn_t)
    return list_of_obs

def ostia_L3_sst(list_of_obs):
    fname   = '/discover/nobackup/projects/gmao/m2oasf/aogcm/g5odas/obs/raw/OSTIA/ASSIM_QRT_SUB/'+yyyy+'/SST_OSTIA_'+yyyy+mm+dd+'.nc'
    ostia_t = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sst_obs'], vartype='TEMP',xsigo=xsigo_sst, descriptor='OSTIA SST', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, ostia_t)
    return list_of_obs

def avhrr18_L2_sst(list_of_obs):
    fname='/discover/nobackup/projects/gmao/m2oasf/aogcm/g5odas/obs/assim/L2_SST_7.0/'+yyyy+'/AVHRR18_G_'+yyyy+mm+dd+'.nc'
    avhrr18 = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_sst_obs'], vartype='TEMP',xsigo=0.1, descriptor='AVHRR18 L2 SST', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, avhrr18)
    return list_of_obs

def noaa16_L2_sst(list_of_obs):
    noaa16 = Obs(yyyy, mm, dd, hh, fname='', id_obs=obsid_dict['id_sst_obs'], vartype='TEMP',xsigo=0.1, descriptor='NOAA-16', platform = 'NOAA16', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, noaa16)
    return list_of_obs

def metopa_L2_sst(list_of_obs):
    metopa = Obs(yyyy, mm, dd, hh, fname='', id_obs=obsid_dict['id_sst_obs'], vartype='TEMP',xsigo=0.1, descriptor='METOP-A', platform = 'METOPA', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, metopa)
    return list_of_obs

def gmi_L2_sst(list_of_obs):
    gmi = Obs(yyyy, mm, dd, hh, fname='', id_obs=obsid_dict['id_sst_obs'], vartype='TEMP',xsigo=0.1, descriptor='GMI-RSS', platform = 'GMI', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, gmi)
    return list_of_obs

def merra2_sst(list_of_obs):
    print('in merra2_sst')
    m2sst = Obs(yyyy, mm, dd, hh, fname='', id_obs=obsid_dict['id_sst_obs'], vartype='TEMP',xsigo=1.0, descriptor='M2-SST', platform = 'M2-SST', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, m2sst)
    return list_of_obs

#def pathfinder(list_of_obs):


# Ice Fraction retrieval
#=======================
def nsidc_aice(list_of_obs):
    fname=OBSDIR+'/AICE_6.0/ICE_NSIDC_'+yyyy+'.nc'
    nsidc_a = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_aice_obs'], vartype='ICE',xsigo=xsigo_aice, descriptor='NASA-TEAM-2', NDAYS=EXP_NDAYS)
    if logit_transform:
        nsidc_a.transform() #logit transform for sea ice fraction
    update_list_of_obs(list_of_obs, nsidc_a)
    return list_of_obs

def reyn_aice(list_of_obs):
    fname=OBSDIR+'/AICE_6.0/ICE_REYN_'+yyyy+'.nc'
    reyn_a = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_aice_obs'], vartype='ICE',xsigo=xsigo_aice, descriptor='NOAA-AICE', NDAYS=EXP_NDAYS)
    if logit_transform:
        reyn_a.transform() #logit transform for sea ice fraction
    update_list_of_obs(list_of_obs, reyn_a)
    return list_of_obs

def ostia_aice(list_of_obs):
    fname   = '/discover/nobackup/projects/gmao/m2oasf/aogcm/g5odas/obs/raw/OSTIA/ASSIM_QRT_SUB/'+yyyy+'/ICE_OSTIA_'+yyyy+mm+dd+'.nc'
    ostia_a = Obs(yyyy, mm, dd, hh, fname=fname, id_obs=obsid_dict['id_aice_obs'], vartype='ICE',xsigo=xsigo_aice, descriptor='OSTIA-AICE', NDAYS=EXP_NDAYS)
    if logit_transform:
        ostia_a.transform() #logit transform for sea ice fraction
    update_list_of_obs(list_of_obs, ostia_a)
    return list_of_obs

# Ice Thickness/Freeboard/Snow depth
#===================================
def oib_hice(list_of_obs):
    oib_hi = Obs(yyyy, mm, dd, hh, fname='', id_obs=obsid_dict['id_hice_obs'], platform='AIR-BORN', xsigo=xsigo_hice, descriptor='OIB', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, oib_hi)
    return list_of_obs

def cs2_hice(list_of_obs):
    cs2_hi = Obs(yyyy, mm, dd, hh, fname='', id_obs=obsid_dict['id_hice_obs'], platform='CS2-HICE', xsigo=xsigo_hice, descriptor='CS2-HICE', NDAYS=EXP_NDAYS)
    update_list_of_obs(list_of_obs, cs2_hi)
    return list_of_obs

switch_obs = {
    "Argo"        : argo,
    "Argo_1"      : argo_1,
    "Argo_2"      : argo_2,
    "Argo_3"      : argo_3,
    "Argo_4"      : argo_4,
    "Argo_5"      : argo_5,
    "CTD"         : ctd,
    "XBT"         : xbt,
    "XBT-SYN-S"   : xbt_synS,
    "WOA18_atArgo": woa18_atargo,
    "WOA18_atNA"  : woa18_atna,
    "TAO"         : tao,
    "PIRATA"      : pirata,
    "RAMA"        : rama,
    "SMOS"        : smos_L2_sss,
    "SMOSSUB"     : smossub_L2_sss,
    "SMOSL3"      : smos_L3_sss,
    "SMAP"        : smap_L2_sss,
    "SMAPV4.1"    : smap_v4_1_L2_sss,
    "SMAPV4.2"    : smap_v4_2_L2_sss,
    "SMAPV4.3"    : smap_v4_3_L2_sss,
    "SMAPV5.0"    : smap_v5_0_L2_sss,
    "SMAPL3"      : smap_L3_sss,
    "AQUARIUS"    : aq_L2_sss,
    "AQUARIUSL3"  : aq_L3_sss,
    "CryoSat-2"   : cryosat2,
    "CryoSat-2-N" : cryosat2_N,
    "Jason-1"     : jason1,
    "Jason-1-N"   : jason1_N,
    "Jason-1G"    : jason1G,
    "Jason-2"     : jason2,
    "Jason-2-N"   : jason2_N,
    "Jason-3"     : jason3,
    "Jason-3N"    : jason3N,
    "Saral"       : saral,
    "Sentinel-3a" : sentinel3a,
    "Sentinel-6a" : sentinel6a,
    "Sentinel-3b" : sentinel3b,
    "SWOTN"       : swotn,
    "ERS-1"       : ers1,
    "ERS-2"       : ers2,
    "TOPEX"       : topex_poseidon,
    "TOPEX-N"     : topex_poseidon_N,
    "GEOSAT-2"    : geosat_follow_on,
    "Envisat"     : envisat,
    "Envisat-N"   : envisat_N,
    "HY-2A"       : hy2a,
    "Reynolds"    : reyn_L3_sst,
    "OSTIA"       : ostia_L3_sst,
    "M2-SST"      : merra2_sst,
    "AVHRR-18"    : avhrr18_L2_sst,
    "NASA-TEAM-2" : nsidc_aice,
    "NOAA-AICE"   : reyn_aice,
    "OSTIA-AICE"  : ostia_aice,
    "NOAA-16"     : noaa16_L2_sst,
    "METOP-A"     : metopa_L2_sst,
    "GMI-RSS"     : gmi_L2_sst,
    "OIB-HICE"    : oib_hice,
    "CS2-HICE"    : cs2_hice
 }

switch_inst_ids = {
    "TAO"         : 501,
    "PIRATA"      : 502,
    "XBT"         : 503,
    "RAMA"        : 504,
#   "ADCP"        : 505,
#   "Curr Meter"  : 506,
    "Argo"        : 508,
#   "Levitus"     : 509,
    "CTD"         : 513,
#   "Reynolds"    : 516,
    "GEOSAT-2"    : 510,
    "ERS-1"       : 511,
    "Envisat"     : 512,
    "TOPEX"       : 514,
    "Jason-1"     : 515,
    "Jason-2"     : 517,
    "CryoSat-2"   : 532,
    "CryoSat-2-N" : 543,
    "Envisat-N"   : 533,
    "Saral/Altika": 534,
    "Jason-1G"    : 535,
    "ERS-2"       : 536,
    "HY-2A"       : 537,
    "Jason-1-N"   : 538,
    "TOPEX-N"     : 539,
    "Jason-2-N"   : 540,
    "Jason-3"     : 541,
    "Jason-3N"    : 545,
    "Sentinel-3a" : 542,
    "Sentinel-6a" : 544,
    "SWOTN"       : 546,
    "Sentinel-3b" : 547,
    "SMOS"        : 555,
#   "SMOSSUB"     : 555,
#   "SMOSL3"      : 556,
#   "SMAP"        : 556,
#   "SMAPV4.1"    : 556,
#   "SMAPV4.2"    : 556,
#   "SMAPV4.3"    : 556,
    "SMAPV5.0"    : 556,
#   "SMAPL3"      : 551,
    "AQUARIUS"    : 553,
    "AQUARIUSL3"  : 550,
    "Levitus SSS" : 521,
#   "NSIDC AICE"  : 518,
    "NASA-TEAM-2" : 518,
    "CMIP5 AICE"  : 519,
    "CMIP5 SST"   : 520,
    "Reynold AICE": 523,
    "NOAA-AICE"   : 523,
    "OSTIA-AICE"  : 526,
    "PIOMAS HICE" : 528,
    "O1B HICE"    : 527,
    "O1B VICE"    : 529,
#   "MODS Chlor"  : 524,
    "Cryosat HICE": 530,
#   "MODIS_GSFC_A": 600,
#   "MODIS_A"     : 601,
#   "AVHRR18_G"   : 602,
    "XBT-SYN-S"   : 531,
#   "OSTIA"       : ostia_L3_sst,
    "WOA18_atArgo": 557,
    "WOA18_atNA"  : 558,
    "M2-SST"      : 516,
#   "AVHRR-18"    : avhrr18_L2_sst,
#   "NOAA-16"     : noaa16_L2_sst,
#   "METOP-A"     : metopa_L2_sst,
#   "GMI-RSS"     : gmi_L2_sst,
#   "OIB-HICE"    : oib_hice,
#   "CS2-HICE"    : cs2_hice
    
 }

yyyy     = sys.argv[1]      # '2012'
mm       = sys.argv[2]      # '12'
dd       = sys.argv[3]      # '01'
hh       = sys.argv[4]      # '01'
obs_type = sys.argv[5:]     #  platform
gmao_data= sys.argv[6:]     #  output directory

print ('top',obs_type)
list_of_obs = []
for doobs in obs_type:
    print('top2',obs_type)
    print(doobs)
    print(switch_obs[doobs])
    list_of_obs = switch_obs[doobs](list_of_obs)

print('============== Extracting obs ============')
if list_of_obs:
    cnt=0
    for obs in list_of_obs:
        print((obs.descriptor))
        pngname='obs-'+mm+'-'+dd+'-'+yyyy
        obs.plot(pngname)
        if (cnt==0):
            typ=obs.typ
            lon=obs.lon
            lat=obs.lat
            depth=obs.depth
            value=obs.value
            oerr=obs.oerr
            instid=obs.instid
            date_time = obs.date_time

            print(np.shape(oerr), np.shape(obs.oerr))
        else:
            typ=np.concatenate( (typ, obs.typ) )
            lon=np.concatenate( (lon,obs.lon) )
            lat=np.concatenate( (lat,obs.lat) )
            depth=np.concatenate( (depth,obs.depth) )
            value=np.concatenate( (value,obs.value) )
            instid=np.concatenate( (instid,obs.instid) )
            date_time=np.concatenate( (date_time, obs.date_time) )
            print(np.shape(oerr), np.shape(obs.oerr))

            oerr=np.concatenate( (oerr,obs.oerr) )

        cnt+=1

#***code to manually fix 180 data ***************************
    for i in range (0,len(lon)):
       if lon[i]==180.:
         lon[i]=-179.9999
       if lon[i]==-180.:
         lon[i]=-179.9999
#************************************************************

#***code to manually fix 24:00z issue for Argo ************************
    for i in range (0,len(date_time)):
        ss=str(int(date_time[i]))
        if float(ss[-2::])==24.0:
           date_time[i]=date_time[i]-1
#********************************************************


    nobs=len(typ)
    print('nobs=',nobs)
#   with open('Nobs', 'wb') as fh:
    with open('Nobs', 'w') as fh:
        fh.write(str(nobs)+'\n')

#   check to make sure minimum observation types are in gmao- file
#   minreq=[3073,5521,5351]  #Tz, Sz, ADT
#   minreq=[3073,5351]  #Tz, ADT
    minreq=[]  #Tz, ADT
#   minreq=[]  #Tz, ADT
#   minreq=[3073,5521]  #Tz, Sz
#   minreq=[3073]  #Tz, Sz
    for i in range(len(minreq)):
      if minreq[i] in typ:
        print("Yes,found in List : ",minreq[i])
      else:
        if minreq[i] == 3073:
           print("YOU ARE MISSING",minreq[i],"Argo T")
           print ('PARENT ID IS', os.getppid(), os.getpid())
        if minreq[i] == 5521:
           print("YOU ARE MISSING",minreq[i],"Argo S")
        if minreq[i] == 5351:
           print("YOU ARE MISSING",minreq[i],"ADT")
#   create a file to tell parents that obs are bad
        command='touch '+SCRDIR+'/BADOBS'
        os.system(command)
        sys.exit(1)
    fnameout = gmao_data + '/gmao-obs-'+yyyy+mm+dd+hh+'-'+obs_type[0]+'.nc'
#    fnameout='/discover/nobackup/lren1/jedi_obs//gmao_ocean_obs/gmao-obs-'+yyyy+mm+dd+hh+'-'+obs_type[0]+'.nc'
    ncfile = Dataset(fnameout,'w')
    ncfile.createDimension('nobs',nobs)

    tmp = ncfile.createVariable('date_time',np.dtype('int32').char,('nobs'))
    tmp[:] = date_time

    tmp = ncfile.createVariable('typ',np.dtype('int32').char,('nobs'))
    tmp[:] = typ

    tmp = ncfile.createVariable('lon',np.dtype('float32').char,('nobs'))
    tmp[:] = lon

    tmp = ncfile.createVariable('lat',np.dtype('float32').char,('nobs'))
    tmp[:] = lat

    tmp = ncfile.createVariable('depth',np.dtype('float32').char,('nobs'))
    tmp[:] = depth

    tmp = ncfile.createVariable('value',np.dtype('float32').char,('nobs'))
    tmp[:] = value

    tmp = ncfile.createVariable('oerr',np.dtype('float32').char,('nobs'))
    tmp[:] = oerr

    tmp = ncfile.createVariable('instid',np.dtype('float32').char,('nobs'))
    tmp[:] = instid

    ncfile.close()
    print('Saved ',str(nobs),'obs in ',fnameout)

#   now replace with the proper gmao- file  DANGEROUS
#    assumes you will be running in the same directory tree
#    (i.e. eh0??/ocean_das/oana-YYYYMMDD_HH/ocean_observer_YYYYMMDD_HH)
#    use_old_obs = os.environ['ODAS_USE_OLD_OBS']
#    if (use_old_obs=='True'):
#        use_old_obs = True
#    else:
#        use_old_obs = False
#   if(os.environ['ODAS_USE_OLD_OBS']):
#    if( use_old_obs == True ):
#        print('TRYING TO REPLACE OBS WITH OLD VERSION')
#        ODAS_DIR_OLD_OBS = os.environ['ODAS_DIR_OLD_OBS']
#        currentpwd = os.getcwd ()
#       lastnode=os.path.split(currentpwd)[-1]
#        path_list = currentpwd.split(os.sep)
#        lastnode = path_list[-1]
#        lastnode2 = path_list[-2]
#        checkFile=ODAS_DIR_OLD_OBS+lastnode2+'/'+lastnode+'/'+fnameout
#        print(checkFile)
#        if os.path.isfile(checkFile):
#            command = 'mv '+fnameout+' '+fnameout+'.didnotuse'
#            print(command)
#            os.system(command)
#            command = 'mv Nobs Nobs.didnotuse'
#            print(command)
#            os.system(command)
#   DANGEROUS!!!!!!!!!!!!!!!!!
#            print(('***REPLACING OBS WITH',ODAS_DIR_OLD_OBS+lastnode2+'/'+lastnode+'/'+fnameout,'**********'))
#          command='cp '+ODAS_DIR_OLD_OBS+lastnode2+'/'+lastnode+'/'+fnameout+' '+fnameout  3/16/21 from Yehui
#            command='cp -f '+ODAS_DIR_OLD_OBS+lastnode2+'/'+lastnode+'/'+fnameout+' '+fnameout
#            print(command)
#            os.system(command)
#          command='cp '+ODAS_DIR_OLD_OBS+lastnode2+'/'+lastnode+'/Nobs Nobs' 3/16/21 from Yehui
#            command='cp -f '+ODAS_DIR_OLD_OBS+lastnode2+'/'+lastnode+'/Nobs Nobs'
#            print(command)
#            os.system(command)
#  now check  to make sure it has instid, if not add a variable instid with 0's
#######################################
#            infile=fnameout
#            outfile='tempout.nc'
#            cmd = "ncdump -v instid "+infile
#          print cmd
#            returned_value = os.system(cmd)
#          print 'return is',returned_value
#            if (returned_value != 0):
#                cmd = "/discover/swdev/gmao_SIteam/Baselibs/latest-mpiuni-SLES12/Linux/bin/ncap2 -s 'instid[$nobs]=-999.' "+infile+' '+outfile
#                  print cmd
#                print('ADDING INST_ID TO NETCDF FILE')
#                os.system(cmd)
#                cmd = "mv "+outfile+' '+fnameout
#                os.system(cmd)
#######################################

#        else:
#            print(('SORRY', checkFile, 'DOES NOT EXIST'))
else:
#    with open('Nobs', 'wb') as fh:
    with open('Nobs', 'w') as fh:
        fh.write('0'+'\n')
#try:
#    command = './oceanobs_nc2bin.x -y '+yyyy+' -m '+mm+' -d '+dd+' -indir1 gmao-obs- -outdir .'
#    print(command)
#    os.system(command)
#except:
#    print('Could not convert to NCEP binary format')
