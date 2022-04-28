!
PROGRAM proc_SST_FRACI_reynolds_quart
!---------------------------------------------------------------------------
        IMPLICIT NONE

        INTEGER, PARAMETER      :: iDebug                 =  0

        REAL,    PARAMETER      :: myUNDEF                = 1.0e15
        REAL,    PARAMETER      :: TempLow                = 273.15d0 ! 0 deg C or low sst (in deg K) below which there is ice
        REAL,    PARAMETER      :: Ice_thr                = 1.0e-4   ! threshold on ice concentration- related to TempLow

        INTEGER, PARAMETER      :: NLAT                   = 720
        INTEGER, PARAMETER      :: NLON                   = 1440

        CHARACTER (LEN = 100)   :: inputBuffer, inputFile
        CHARACTER (LEN = 150)   :: fileName, maskFileName
        CHARACTER (LEN = 8)     :: today, tomrw

        CHARACTER (LEN = 40)    :: fileName_reynolds_SST, fileName_reynolds_ICE

        INTEGER                 :: iERR
        INTEGER                 :: NLAT_out
        INTEGER                 :: NLON_out
        INTEGER                 :: iLon, iLat, k

        REAL                    :: reynolds_LAT (NLAT), reynolds_LON(NLON)
        REAL                    :: reynolds_SST (NLON,  NLAT)
        REAL                    :: reynolds_ICE (NLON,  NLAT)
        REAL                    :: reynolds_MASK(NLON,  NLAT)

        REAL                    :: sstave

        REAL                    :: HEADER(14)
        CHARACTER(LEN = 4)      :: today_Year, tomrw_Year
        CHARACTER(LEN = 2)      :: today_Mon,  tomrw_Mon, today_Day, tomrw_Day
        INTEGER                 :: today_iYear, tomrw_iYear
        INTEGER                 :: today_iMon,  tomrw_iMon, today_iDay, tomrw_iDay
!       ....................................................................

!---------------------------------------------------------------------------
!       Read all input data parameters (time to proc, files to proc, output resolution)
        CALL getarg(1,inputBuffer)
        READ(inputBuffer, *) inputFile

        OPEN (UNIT = 21, FILE = inputFile, STATUS = 'old')
        READ (21, '(A)')  today
        READ (21, '(A)')  tomrw
        READ (21, '(A)')  fileName       ! Data set with SST and sea ice concentration
        READ (21, '(A)')  maskFileName   ! Data set with land-sea mask 0: ocean, 1:land
        READ (21, '(I4)') NLAT_out
        READ (21, '(I4)') NLON_out
        CLOSE(21)

        iERR = 0
        IF( today == tomrw)  THEN
          iERR = 1
          PRINT *, 'Processing Start date: ', today
          PRINT *, 'is SAME as End date:   ', tomrw
          PRINT *, 'End date must be AFTER Start date'
        END IF
!---------------------------------------------------------------------------
        IF( iERR == 0) THEN
             PRINT *, 'Processing SST and ICE data @ 1/4 deg from: ', today, '...To... ', tomrw
        ELSE
             PRINT *, 'User input is not in correct format- for this program to work!'
             PRINT *, 'SEE ABOVE LOG FOR DETAILS'
             STOP
        END IF

!------------------------------Read input files----------------------------
!       Read Reynolds 
!       SST               -> reynolds_SST_native 
!       Ice concentration -> reynolds_ICE_native
        CALL read_Reynolds(fileName, maskFileName, NLAT, NLON,           &
                           reynolds_LAT, reynolds_LON,                   &
                           reynolds_SST, reynolds_ICE, reynolds_MASK, myUNDEF)

!------------------------------Process SST & ICE fields--------------------
!       reynolds ice has undef in open water as well. make that to 0.
!       WHERE( (reynolds_SST .ne. myUNDEF) .and. (reynolds_ICE == myUNDEF)) 
!            reynolds_ICE = 0.0d0
!       END WHERE 
        WHERE( (reynolds_MASK == 0.) .and. (reynolds_ICE == myUNDEF)) 
!       WHERE( reynolds_MASK == 1.0) ! mask: land == 1.
             reynolds_ICE = 0.0d0
!            reynolds_ICE = myUNDEF
        END WHERE 

!       Reynolds(SST, ICE): (1) flip
        CALL hflip              ( reynolds_SST, NLON, NLAT )
        CALL hflip              ( reynolds_ICE, NLON, NLAT )   
!---------------------------------------------------------------------------

        CALL fill_Land (reynolds_SST,    NLON_out, NLAT_out, myUNDEF) 
        CALL fill_Land (reynolds_ICE,    NLON_out, NLAT_out, myUNDEF) 

!---------------------------------------------------------------------------
! SST values over Antarctic land - for ice, it does not matter which way, since it is over *land*
!---------------------------------------------------------------------------

        sstave = 0.0d0
        DO iLon = 1, NLON_out
           iLat = 1
           DO WHILE( reynolds_SST(iLon, iLat) .EQ. myUNDEF)
              iLat = iLat + 1
           END DO
           sstave = sstave + reynolds_SST(iLon, iLat)
        END DO
        sstave = sstave/NLON_out
        DO iLon = 1, NLON_out
           iLat = 1
           DO WHILE( reynolds_SST(iLon, iLat) .EQ. myUNDEF)
              reynolds_SST(iLon, iLat) = sstave
              iLat = iLat + 1
           END DO
        END DO

       DO iLon = 1, NLON_out
           iLat = 1
           DO WHILE( reynolds_ICE(iLon, iLat) .EQ. myUNDEF)
              iLat = iLat + 1
           END DO
           DO k = 1, iLat-1
              reynolds_ICE(iLon, k) = reynolds_ICE(iLon, iLat)
           END DO
        END DO
!---------------------------------------------------------------------------

!       Header info.  Start & end dates: format: YYYYMMDDHHMMSS; Hour,min,Sec are set to zero.
        today_Year    = today(1:4);      tomrw_Year    = tomrw(1:4)
        today_Mon     = today(5:6);      tomrw_Mon     = tomrw(5:6)
        today_Day     = today(7:8);      tomrw_Day     = tomrw(7:8)

        READ( today_Year, 98) today_iYear
        READ( tomrw_Year, 98) tomrw_iYear

        READ( today_Mon,  99) today_iMon
        READ( tomrw_Mon,  99) tomrw_iMon

        READ( today_Day,  99) today_iDay
        READ( tomrw_Day,  99) tomrw_iDay

        HEADER(1)    = REAL(today_iYear); HEADER(7)     = REAL(tomrw_iYear)
        HEADER(2)    = REAL(today_iMon);  HEADER(8)     = REAL(tomrw_iMon)
        HEADER(3)    = REAL(today_iDay);  HEADER(9)     = REAL(tomrw_iDay)
        HEADER(4:6)  = 0.0;               HEADER(10:12) = 0.0

        HEADER(13)   = REAL(NLON_out);    HEADER(14)    = REAL(NLAT_out)
!---------------------------------------------------------------------------

!       Write out Reynolds fields for MERRA-x
!       SST, ICE:
        fileName_reynolds_SST  = 'Reynolds_sst_' // today //'.bin'
        fileName_reynolds_ICE  = 'Reynolds_ice_' // today //'.bin'

        OPEN (UNIT = 993, FILE = fileName_reynolds_SST, FORM = 'unformatted', STATUS = 'new')
        OPEN (UNIT = 994, FILE = fileName_reynolds_ICE, FORM = 'unformatted', STATUS = 'new')

        WRITE(993) HEADER
        WRITE(994) HEADER
        WRITE(993) reynolds_SST
        WRITE(994) reynolds_ICE
        CLOSE(993)
        CLOSE(994)
!---------------------------------------------------------------------------

        IF( iERR == 0) PRINT *, '...Finished!'
!---------------------------------------------------------------------------

 98     FORMAT(I4)
 99     FORMAT(I4)
!---------------------------------------------------------------------------
END PROGRAM proc_SST_FRACI_reynolds_quart
!
