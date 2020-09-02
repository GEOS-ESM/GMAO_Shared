!
PROGRAM lake_data_EIGTHdeg
!---------------------------------------------------------------------------
        IMPLICIT NONE

        INTEGER, PARAMETER      :: iDebug                 =  0

        REAL,    PARAMETER      :: myUNDEF                = 1.0e15
        REAL,    PARAMETER      :: TempLow                = 273.15d0 ! low sst (in deg K) below which there is ice
        REAL,    PARAMETER      :: Ice_thr                = 1.0e-4   ! threshold on ice concentration- related to TempLow

        INTEGER, PARAMETER      :: reynolds_NLAT          = 720
        INTEGER, PARAMETER      :: reynolds_NLON          = 1440

        INTEGER, PARAMETER      :: ostia_NLAT             = 3600
        INTEGER, PARAMETER      :: ostia_NLON             = 7200

        CHARACTER (LEN = 100)   :: inputBuffer, inputFile
        CHARACTER (LEN = 150)   :: fileNames(2)
        CHARACTER (LEN = 8)     :: today, tomrw

        CHARACTER (LEN = 150)   :: fileName_Reynolds,  fileName_Ostia
        CHARACTER (LEN = 40)    :: fileName_ostia_SST, fileName_ostia_ICE
        CHARACTER (LEN = 40)    :: fileName_mask

        INTEGER                 :: iERR, iMerra
        INTEGER                 :: NLAT_out
        INTEGER                 :: NLON_out
        INTEGER                 :: iLon, iLat
        INTEGER                 :: iAdjust_SST_SIC
        REAL                    :: SST_thr                           ! threshold on SST- above which we "expect" no sea-ice 

        REAL                    :: reynolds_LAT        (reynolds_NLAT), reynolds_LON(reynolds_NLON)
        REAL                    :: reynolds_SST_native (reynolds_NLON,  reynolds_NLAT)
        REAL                    :: reynolds_ICE_native (reynolds_NLON,  reynolds_NLAT)

        REAL, ALLOCATABLE       :: reynolds_SST_eigth(:,:), reynolds_ICE_eigth(:,:)
        REAL, ALLOCATABLE       :: ostia_SST_native  (:,:),   ostia_SST_eigth (:,:)
        REAL, ALLOCATABLE       :: ostia_ICE_native  (:,:),   ostia_ICE_eigth (:,:)
        REAL, ALLOCATABLE       :: mask(:,:)

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
        CALL read_input(inputFile, iDebug, today, tomrw, fileNames, NLAT_out, NLON_out, iMerra, iAdjust_SST_SIC, SST_Thr, iERR)
!---------------------------------------------------------------------------
        IF( iERR == 0) THEN
             PRINT *, 'Data over Great Lakes and Caspian Sea'
             PRINT *, 'Processing SST and ICE data from: ', today, '...To... ', tomrw
        ELSE
             PRINT *, 'User input is not in correct format- for this program to work!'
             PRINT *, 'SEE ABOVE LOG FOR DETAILS'
             STOP
        END IF

        fileName_Reynolds = fileNames(1)
        fileName_Ostia    = fileNames(2)
!------------------------------Read input files----------------------------
!       Read Reynolds 
!       SST               -> reynolds_SST_native 
!       Ice concentration -> reynolds_ICE_native
        CALL read_Reynolds(fileName_Reynolds, reynolds_NLAT, reynolds_NLON,                        &
                           reynolds_LAT, reynolds_LON,                                             &
                           reynolds_SST_native, reynolds_ICE_native, myUNDEF)

!       Read Ostia 
!       SST               -> ostia_SST_native
!       Ice concentration -> ostia_ICE_native
        ALLOCATE( ostia_SST_native(ostia_NLON, ostia_NLAT) )
        ALLOCATE( ostia_ICE_native(ostia_NLON, ostia_NLAT) )
        CALL read_Ostia( fileName_Ostia, "analysed_sst",     ostia_NLAT, ostia_NLON, ostia_SST_native, myUNDEF )
        CALL read_Ostia( fileName_Ostia, "sea_ice_fraction", ostia_NLAT, ostia_NLON, ostia_ICE_native, myUNDEF )

!------------------------------Process SST & ICE fields--------------------
!       reynolds ice has undef in open water as well. make that to 0.
        WHERE( (reynolds_SST_native .ne. myUNDEF) .and. (reynolds_ICE_native == myUNDEF))
             reynolds_ICE_native = 0.0d0
        END WHERE

!       Reynolds(SST, ICE): (1) flip, (2) interp to 1/8 deg
        ALLOCATE( reynolds_SST_eigth(NLON_out,      NLAT_out     ))
        ALLOCATE( reynolds_ICE_eigth(NLON_out,      NLAT_out     ))

        CALL hflip              ( reynolds_SST_native, reynolds_NLON, reynolds_NLAT )
        CALL hflip              ( reynolds_ICE_native, reynolds_NLON, reynolds_NLAT )

        CALL interp_to_eight_deg( reynolds_SST_native, reynolds_NLON, reynolds_NLAT,                &
                                  reynolds_SST_eigth, NLON_out, NLAT_out, myUNDEF)
        CALL interp_to_eight_deg( reynolds_ICE_native, reynolds_NLON, reynolds_NLAT,                &
                                  reynolds_ICE_eigth, NLON_out, NLAT_out, myUNDEF)

!---------------------------------------------------------------------------
        ALLOCATE( ostia_SST_eigth(NLON_out,  NLAT_out ))
        ALLOCATE( ostia_ICE_eigth(NLON_out,  NLAT_out ))
        ALLOCATE( mask            (NLON_out,  NLAT_out ))

!       Initialize
        mask = 0.

!       Ostia(SST, ICE): bin to 1/8 deg.
        CALL bin2bin( ostia_SST_native, ostia_NLON, ostia_NLAT, ostia_SST_eigth, NLON_out, NLAT_out, myUNDEF )
        CALL bin2bin( ostia_ICE_native, ostia_NLON, ostia_NLAT, ostia_ICE_eigth, NLON_out, NLAT_out, myUNDEF )
!---------------------------------------------------------------------------
!       Get Great Lakes SST and ICE from Reynolds into OSTIA (if needed)
        DO iLat = 1046, 1120
         DO iLon = 695, 842

            mask(iLon,iLat) = 1.0

            IF( (ostia_SST_eigth(iLon,iLat).eq.myUNDEF) .and. (reynolds_SST_eigth(iLon,iLat).ne.myUNDEF) ) THEN
                 ostia_SST_eigth(iLon,iLat) = reynolds_SST_eigth(iLon,iLat)
            END IF

            IF( (ostia_ICE_eigth(iLon,iLat).eq.myUNDEF) .and. (reynolds_ICE_eigth(iLon,iLat).ne.myUNDEF) ) THEN
                 ostia_ICE_eigth(iLon,iLat) = reynolds_ICE_eigth(iLon,iLat)    ! if OSTIA had no ice in Great Lakes, get data from Reynolds
            END IF
         END DO
        END DO
!---------------------------------------------------------------------------
!       Caspian Sea ice: there is no ice info in OSTIA ICE, when temp < freezing point, fix this problem.
        DO iLat = 1000, 1120
         DO iLon = 1800, 1890

            mask(iLon,iLat) = 1.0

            ! 1st. handle SST. if OSTIA SST = undef, then use Reynolds SST
            IF( (ostia_SST_eigth(iLon,iLat).eq.myUNDEF) .and. (reynolds_SST_eigth(iLon,iLat).ne.myUNDEF) ) THEN
                 ostia_SST_eigth(iLon,iLat) = reynolds_SST_eigth(iLon,iLat)
            END IF

            ! if sst < freezing temp
            IF( ostia_SST_eigth(iLon,iLat) <= 275.0d0) THEN
              ! if there is no ice data or ice ~ 0.0 
              IF( (ostia_ICE_eigth(iLon,iLat) .eq. myUNDEF) .or. (ostia_ICE_eigth(iLon,iLat) <= Ice_thr)) THEN
                  IF( (reynolds_ICE_eigth(iLon,iLat) .ne. myUNDEF) .and. (reynolds_ICE_eigth(iLon,iLat) > Ice_thr)) THEN
                     ostia_ICE_eigth(iLon,iLat)  = reynolds_ICE_eigth(iLon,iLat)                         ! use Reynolds ice
                  ELSE                                                                                   ! there is nothing useful from Reynolds, use empirical fit
                     ! this could be computed based on SST from OSTIA/Reynolds.
                     ostia_ICE_eigth(iLon,iLat)  = MIN( 1.0d0, MAX(-0.017451*((ostia_SST_eigth(iLon,iLat)- 271.38)/0.052747) + 0.96834, 0.0d0))
                  END IF
              END IF
            END IF  ! IF( ostia_SST_eigth(iLon,iLat) <= 275.0d0)
         END DO
        END DO
!---------------------------------------------------------------------------
        
        IF ( iAdjust_SST_SIC == 1) THEN
!          adjust SIC based on a threshold value of SST
           WHERE (ostia_SST_eigth > SST_Thr)
               ostia_ICE_eigth = 0.0d0
           END WHERE
        END IF
!---------------------------------------------------------------------------

!       get rid of the masked values
        WHERE( (ostia_SST_eigth <= 260.) .or. (ostia_SST_eigth >= 330.))
               ostia_SST_eigth = -999.
        ENDWHERE

        WHERE (ostia_SST_eigth == -999.)
               ostia_ICE_eigth = -999.
        ENDWHERE

        WHERE( (ostia_ICE_eigth < 0.0) .or. (ostia_ICE_eigth > 1.0))
               ostia_ICE_eigth = -999.
        ENDWHERE
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
!       Write out for OSTIA fields for FP & RPIT
!       SST, ICE:
        fileName_ostia_SST = 'Ostia_sst_'        // today //'.bin'
        fileName_ostia_ICE = 'Ostia_ice_'        // today //'.bin'
        fileName_mask      = 'mask'                       //'.bin'

        OPEN (UNIT = 991, FILE = fileName_ostia_SST, FORM = 'unformatted', STATUS = 'new')
        OPEN (UNIT = 992, FILE = fileName_ostia_ICE, FORM = 'unformatted', STATUS = 'new')
        OPEN (UNIT = 993, FILE = fileName_mask,      FORM = 'unformatted', STATUS = 'new')

        WRITE(991) HEADER
        WRITE(992) HEADER
        WRITE(993) HEADER
        WRITE(991) ostia_SST_eigth
        WRITE(992) ostia_ICE_eigth
        WRITE(993) mask
        CLOSE(991)
        CLOSE(992)
        CLOSE(993)
!---------------------------------------------------------------------------
        IF( iERR == 0) PRINT *, '...Finished!'
!---------------------------------------------------------------------------
 98     FORMAT(I4)
 99     FORMAT(I4)
!---------------------------------------------------------------------------
        DEALLOCATE(reynolds_SST_eigth)
        DEALLOCATE(reynolds_ICE_eigth)

        DEALLOCATE(ostia_SST_native)
        DEALLOCATE(ostia_ICE_native)

        DEALLOCATE(mask)

        DEALLOCATE(ostia_ICE_eigth)
        DEALLOCATE(ostia_SST_eigth)
!---------------------------------------------------------------------------
END PROGRAM lake_data_EIGTHdeg
!
