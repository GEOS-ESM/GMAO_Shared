#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

PROGRAM regrid_DYNVEG

  use MAPL_ConstantsMod,only: pi => MAPL_PI

  implicit none

  INCLUDE 'mpif.h'
  INCLUDE 'netcdf.inc'

  INTEGER, PARAMETER :: NUM_VEG = 4, NUM_ZONE = 3,  ntiles_cn = 686448, npft = 19, Noff_per_pft = 64, NFILES = 24
  character(len=300), parameter :: &
       InCNTilFile = '/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/' &
       //'DYNVEG/CF0180x6C_TM0720xTM0410/til/CF0180x6C_TM0720xTM0410-Pfafstetter.til', &
       InBCSDIR   = '/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/' &
       //'DYNVEG/CF0180x6C_TM0720xTM0410/' , &
       GLDASName  = 'CF180_DVG2_C4'

       !      InCNTilFile = '/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/' &
       !      //'DYNVEG/DE_02880x01440/til/DE_02880x01440_PE_2880x1440.til', &
       !      InBCSDIR   = '/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/' &
       !      //'DYNVEG/DE_02880x01440/'
  REAL, parameter :: fmin= 1.e-4, undef = -9999. ! ignore vegetation fractions at or below this value
  type :: tile_data_type
          
     integer :: tile_id    ! unique tile ID
     real    :: com_lon    ! center-of-mass longitude
     real    :: com_lat    ! center-of-mass latitude
     real    :: min_lon    ! minimum longitude (bounding box for tile)
     real    :: max_lon    ! maximum longitude (bounding box for tile)
     real    :: min_lat    ! minimum latitude (bounding box for tile)
     real    :: max_lat    ! maximum latitude (bounding box for tile)
     real, dimension (NUM_VEG) :: ITY     ! vegetation_type
     real, dimension (NUM_VEG) :: FVG     ! vegetation_fraction 

  end type tile_data_type

!  type :: input_data 
!     integer :: ncells
!     integer, allocatable, dimension (:) :: input_id
!     real,    allocatable, dimension (:) :: frac_area
!  end type input_data
!
!  type :: regrid_map
!     type (input_data)    :: map_pt1
!     type (input_data)    :: map_pt2
!     type (input_data)    :: map_st1
!     type (input_data)    :: map_st2
!  end type regrid_map
!  type(regrid_map), allocatable, dimension (:) :: regrid_global

  type(tile_data_type), pointer, contiguous    :: out_tile_data  (:)=>null()
  integer  :: myid=0, numprocs=1, STATUS, mpistatus(MPI_STATUS_SIZE)  
  logical  :: master_proc=.true., all_found = .true., ACTUAL_DV = .true., donot_regrid = .true.
  integer, allocatable, dimension (:)     :: low_ind, upp_ind, nt_local, tid_offl
  integer, allocatable, dimension (:)     :: sub_tid, sub_ityp1, sub_ityp2, sub_ityp3, sub_ityp4, CNT_TYP
  integer, allocatable, dimension (:,:)   :: t_count
  integer, allocatable, dimension (:,:)   :: CNT_TYP1, CNT_TYP2, CNT_TYP3, CNT_TYP4, INTPUT_TID
  integer, allocatable, dimension (:,:,:) :: INTPUT_TID1, INTPUT_TID2, INTPUT_TID3, INTPUT_TID4
  real   , allocatable, dimension (:)     :: lonc, latc, sub_fveg1, sub_fveg2, sub_fveg3, sub_fveg4, sub_area, tile_area
  real   , allocatable, dimension (:)     :: CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, CLMC_pt1, CLMC_pt2,CLMC_st1
  real   , allocatable, dimension (:)     :: CLMC_st2, LATT, LONN, min_lon, max_lon, min_lat, max_lat
  real   , allocatable, dimension (:,:)   :: AREA_TYP, fveg_offl,  ityp_offl
  real   , allocatable, dimension (:,:,:) :: AREA_TYP1, AREA_TYP2, AREA_TYP3, AREA_TYP4
  logical, allocatable, dimension(:)      :: mask

  INTEGER            :: YEARB, YEARE, YEAR, MONTH, DAY, NTILES, NTILES_G, i,n,nv, nplus, ityp_new, loc_size
  CHARACTER*8        :: YYYYMMDD
  CHARACTER*4        :: YYYY
  CHARACTER*2        :: MM
  CHARACTER*300      :: MODELING_SYSTEM, BCSPATH, TILFILE, EXPDIR, ACTORCLIM, BCSDIR
  real               :: dw, this_min_lon, this_max_lon, this_min_lat, this_max_lat, fveg_new
  integer            :: file_units (NFILES), NCInID(2), YR1, MN1, DY1 
  logical, parameter :: clim_sai = .true.

  CHARACTER*7, dimension (NFILES), PARAMETER :: OUT_VARS = (/     &
        'CNSAI11',  'CNSAI12',  'CNSAI13', &  
        'CNSAI21',  'CNSAI22',  'CNSAI23', &   
        'CNSAI31',  'CNSAI32',  'CNSAI33', &   
        'CNSAI41',  'CNSAI42',  'CNSAI43', &   
        'CNLAI11',  'CNLAI12',  'CNLAI13', &    
        'CNLAI21',  'CNLAI22',  'CNLAI23', &   
        'CNLAI31',  'CNLAI32',  'CNLAI33', &   
        'CNLAI41',  'CNLAI42',  'CNLAI43'/)   

  call init_MPI()

  ! READ ENVIROENMENT VARIABLES AND Out TileData
  ! ============================================
  
  call getenv ("MODELING_SYSTEM"  , MODELING_SYSTEM)
  call getenv ("BCSDIR"           , BCSDIR         )
  call getenv ("BCSPATH"          , BCSPATH        )
  call getenv ("TILFILE"          , TILFILE        )
  call getenv ("EXPDIR"           , EXPDIR         )
  call getenv ("ACTUAL_DV"        , ACTORCLIM      )
  
  call getenv ("YEARB"  , YYYY)
  read (YYYY,'(i4.4)') YEARB
  call getenv ("YEARE"  , YYYY)
  read (YYYY,'(i4.4)') YEARE
  if(trim(ACTORCLIM) == 'FALSE') ACTUAL_DV = .false.

  if (master_proc)  then

     print *,trim(YYYY), trim(BCSDIR)

     ! Check whether BCSGRID is CF180
     OPEN (11, FILE = TRIM(BCSDIR)//'/clsm/catchment.def',      FORM = 'FORMATTED', STATUS = 'OLD', ACTION = 'READ')     
     read (11, *) NTILES_G
     CLOSE (11, STATUS = 'KEEP')

  endif

  call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
  call MPI_BCAST(NTILES_G,     1,MPI_INTEGER,0,MPI_COMM_WORLD,STATUS)

  if(NTILES_G /= ntiles_cn) donot_regrid = .false.
  
  if (master_proc) call read_OUT_TileData (NTILES, out_tile_data, MODELING_SYSTEM, BCSPATH, TILFILE, EXPDIR)     

  call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
  call MPI_BCAST(NTILES,     1,MPI_INTEGER,0,MPI_COMM_WORLD,STATUS)

  DO_REGRID : if(.not. donot_regrid) then

     allocate(low_ind (numprocs))
     allocate(upp_ind (numprocs))
     allocate(nt_local(numprocs))
     
     low_ind (:)    = 1
     upp_ind (:)    = NTILES       
     nt_local(:)    = NTILES 
     
     ! Domain decomposition
     ! --------------------
     
     if (numprocs > 1) then      
        do i = 1, numprocs - 1
           upp_ind(i)   = low_ind(i) + (ntiles/numprocs) - 1 
           low_ind(i+1) = upp_ind(i) + 1
           nt_local(i)  = upp_ind(i) - low_ind(i) + 1
        end do
        nt_local(numprocs) = upp_ind(numprocs) - low_ind(numprocs) + 1
     endif
     
     allocate (lonn    (nt_local (myid + 1)))
     allocate (latt    (nt_local (myid + 1)))
     allocate (CLMC_pf1(nt_local (myid + 1)))
     allocate (CLMC_pf2(nt_local (myid + 1)))
     allocate (CLMC_sf1(nt_local (myid + 1)))
     allocate (CLMC_sf2(nt_local (myid + 1)))
     allocate (CLMC_pt1(nt_local (myid + 1)))
     allocate (CLMC_pt2(nt_local (myid + 1)))
     allocate (CLMC_st1(nt_local (myid + 1)))
     allocate (CLMC_st2(nt_local (myid + 1)))
     allocate (min_lon (nt_local (myid + 1)))
     allocate (max_lon (nt_local (myid + 1)))
     allocate (min_lat (nt_local (myid + 1)))
     allocate (max_lat (nt_local (myid + 1)))
     
     loc_size = size (lonn)
     
     call MPI_SCATTERV (out_tile_data%com_lon,nt_local,low_ind-1,MPI_real, lonn    ,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS) 
     call MPI_SCATTERV (out_tile_data%com_lat,nt_local,low_ind-1,MPI_real, latt    ,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS) 
     call MPI_SCATTERV (out_tile_data%min_lon,nt_local,low_ind-1,MPI_real, min_lon ,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS) 
     call MPI_SCATTERV (out_tile_data%max_lon,nt_local,low_ind-1,MPI_real, max_lon ,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)  
     call MPI_SCATTERV (out_tile_data%min_lat,nt_local,low_ind-1,MPI_real, min_lat ,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)  
     call MPI_SCATTERV (out_tile_data%max_lat,nt_local,low_ind-1,MPI_real, max_lat ,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)  
     call MPI_SCATTERV (out_tile_data%FVG(1) ,nt_local,low_ind-1,MPI_real, CLMC_pf1,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)  
     call MPI_SCATTERV (out_tile_data%FVG(2) ,nt_local,low_ind-1,MPI_real, CLMC_pf2,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)  
     call MPI_SCATTERV (out_tile_data%FVG(3) ,nt_local,low_ind-1,MPI_real, CLMC_sf1,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)   
     call MPI_SCATTERV (out_tile_data%FVG(4) ,nt_local,low_ind-1,MPI_real, CLMC_sf2,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)  
     call MPI_SCATTERV (out_tile_data%ITY(1) ,nt_local,low_ind-1,MPI_real, CLMC_pt1,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)  
     call MPI_SCATTERV (out_tile_data%ITY(2) ,nt_local,low_ind-1,MPI_real, CLMC_pt2,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)  
     call MPI_SCATTERV (out_tile_data%ITY(3) ,nt_local,low_ind-1,MPI_real, CLMC_st1,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)  
     call MPI_SCATTERV (out_tile_data%ITY(4) ,nt_local,low_ind-1,MPI_real, CLMC_st2,loc_size,MPI_real,0,MPI_COMM_WORLD, STATUS) ; VERIFY_(STATUS)   
     
     call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
               
     ! CONSTRUCT WEIGHT METRIX FOR REGRIDDING 
     ! ======================================
     
     allocate (tid_offl  (ntiles_cn))
     allocate (mask      (ntiles_cn))
     allocate (ityp_offl (ntiles_cn,num_veg))
     allocate (fveg_offl (ntiles_cn,num_veg))
     allocate (lonc      (1:ntiles_cn))
     allocate (latc      (1:ntiles_cn))
     allocate (tile_area (1:ntiles_cn))
     
     if (master_proc) then
        
        call ReadCNTilFile  (InCNTilFile,ntiles_cn,lonc,latc, tile_area = tile_area) 
        call ReadCLMvegFile (TRIM(InBCSDIR)//'/clsm/CLM_veg_typs_fracs', NTILES_CN,  &
             ityp_offl(:,1), ityp_offl(:,2), ityp_offl(:,3), ityp_offl(:,4),         &
             fveg_offl(:,1), fveg_offl(:,2), fveg_offl(:,3), fveg_offl(:,4))
        
        DO N = 1, ntiles_cn
           
           tid_offl(N) = n
           
           do nv = 1,num_veg
              if(ityp_offl(n,nv)<0 .or. ityp_offl(n,nv)>npft)    stop 'ityp'
              if((fveg_offl(n,nv)<0..or. fveg_offl(n,nv)>1.00001)) stop 'fveg'
           end do
           
           if((ityp_offl(n,3) == 0).and.(ityp_offl(n,4) == 0)) then
              if(ityp_offl(n,1) /= 0) then
                 ityp_offl(n,3) = ityp_offl(n,1)
              else
                 ityp_offl(n,3) = ityp_offl(n,2)
              endif
           endif
           
           if((ityp_offl(n,1) == 0).and.(ityp_offl(n,2) /= 0)) ityp_offl(n,1) = ityp_offl(n,2)
           if((ityp_offl(n,2) == 0).and.(ityp_offl(n,1) /= 0)) ityp_offl(n,2) = ityp_offl(n,1)
           if((ityp_offl(n,3) == 0).and.(ityp_offl(n,4) /= 0)) ityp_offl(n,3) = ityp_offl(n,4)
           if((ityp_offl(n,4) == 0).and.(ityp_offl(n,3) /= 0)) ityp_offl(n,4) = ityp_offl(n,3)
           
        END DO
     endif
     
     call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
     call MPI_BCAST(lonc,     ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,STATUS)
     call MPI_BCAST(latc,     ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,STATUS)
     call MPI_BCAST(tile_area,ntiles_cn,MPI_REAL,0,MPI_COMM_WORLD,STATUS)
     
     call MPI_BCAST(tid_offl ,size(tid_offl ),MPI_INTEGER,0,MPI_COMM_WORLD,STATUS)
     call MPI_BCAST(ityp_offl,size(ityp_offl),MPI_REAL   ,0,MPI_COMM_WORLD,STATUS)
     call MPI_BCAST(fveg_offl,size(fveg_offl),MPI_REAL   ,0,MPI_COMM_WORLD,STATUS)    
     
     allocate (t_count  (nt_local (myid + 1),4))
     allocate (CNT_TYP1 (nt_local (myid + 1),4))
     allocate (CNT_TYP2 (nt_local (myid + 1),4))
     allocate (CNT_TYP3 (nt_local (myid + 1),4))
     allocate (CNT_TYP4 (nt_local (myid + 1),4))
     allocate (AREA_TYP1(nt_local (myid + 1),4,Noff_per_pft))
     allocate (AREA_TYP2(nt_local (myid + 1),4,Noff_per_pft))
     allocate (AREA_TYP3(nt_local (myid + 1),4,Noff_per_pft))
     allocate (AREA_TYP4(nt_local (myid + 1),4,Noff_per_pft))
     allocate (INTPUT_TID1 (nt_local (myid + 1),4,Noff_per_pft))
     allocate (INTPUT_TID2 (nt_local (myid + 1),4,Noff_per_pft))
     allocate (INTPUT_TID3 (nt_local (myid + 1),4,Noff_per_pft))
     allocate (INTPUT_TID4 (nt_local (myid + 1),4,Noff_per_pft))
     allocate (CNT_TYP    (4))
     allocate (AREA_TYP   (4,Noff_per_pft))
     allocate (INTPUT_TID (4,Noff_per_pft))
     
     t_count  = -9999
     CNT_TYP1 = -9999
     CNT_TYP2 = -9999
     CNT_TYP3 = -9999
     CNT_TYP4 = -9999
     AREA_TYP1= -9999.
     AREA_TYP2= -9999.
     AREA_TYP3= -9999.
     AREA_TYP4= -9999.
     INTPUT_TID1 = -9999
     INTPUT_TID2 = -9999
     INTPUT_TID3 = -9999
     INTPUT_TID4 = -9999
     
     TILES : do n = 1, nt_local (myid + 1)
        
        ! First try the rectangular window encompassing the tile
        
        this_min_lon = min_lon(n)
        this_max_lon = max_lon(n)
        this_min_lat = min_lat(n)
        this_max_lat = max_lat(n)
        mask = .false.
        mask =  ((latc >= this_min_lat .and. latc <= this_max_lat).and.(lonc >= this_min_lon .and. lonc <= this_max_lon))
        nplus =  count(mask = mask)
        
        dw = 0.5 ! Start with a 1x1 window, then zoom out by increasing the size by 2-deg until 4 similar tiles are found for 4 PFT types        
        
        ZOOMOUT : do
           
           if(nplus < 1) then
              this_min_lon = MAX(lonn (n) - dw, -180.)
              this_max_lon = MIN(lonn (n) + dw,  180.)
              this_min_lat = MAX(latt (n) - dw,  -90.)
              this_max_lat = MIN(latt (n) + dw,   90.) 
              mask = .false.
              mask =  ((latc >= this_min_lat .and. latc <= this_max_lat).and.(lonc >= this_min_lon .and. lonc <= this_max_lon))
              nplus =  count(mask = mask)
              dw = dw + 1.0
              CYCLE
           endif
           
           allocate (sub_tid   (1:nplus))
           allocate (sub_ityp1 (1:nplus))
           allocate (sub_ityp2 (1:nplus))
           allocate (sub_ityp3 (1:nplus))
           allocate (sub_ityp4 (1:nplus))
           allocate (sub_fveg1 (1:nplus))
           allocate (sub_fveg2 (1:nplus))
           allocate (sub_fveg3 (1:nplus))
           allocate (sub_fveg4 (1:nplus))
           allocate (sub_area  (1:nplus))
           
           sub_tid   = PACK (tid_offl, mask= mask) 
           sub_ityp1 = NINT(ityp_offl (sub_tid,1))
           sub_ityp2 = NINT(ityp_offl (sub_tid,2))
           sub_ityp3 = NINT(ityp_offl (sub_tid,3))
           sub_ityp4 = NINT(ityp_offl (sub_tid,4))
           sub_fveg1 = fveg_offl (sub_tid,1)
           sub_fveg2 = fveg_offl (sub_tid,2)
           sub_fveg3 = fveg_offl (sub_tid,3)
           sub_fveg4 = fveg_offl (sub_tid,4)
           sub_area  = tile_area (sub_tid)
           
           NVLOOP : do nv = 1, num_veg
              
              if (nv == 1) ityp_new = CLMC_pt1(n)
              if (nv == 1) fveg_new = CLMC_pf1(n)
              if (nv == 2) ityp_new = CLMC_pt2(n)
              if (nv == 2) fveg_new = CLMC_pf2(n)
              if (nv == 3) ityp_new = CLMC_st1(n)
              if (nv == 3) fveg_new = CLMC_sf1(n)
              if (nv == 4) ityp_new = CLMC_st2(n) 
              if (nv == 4) fveg_new = CLMC_sf2(n)
              
              SEEK : if((t_count (n,nv) < 0).and.(fveg_new > fmin)) then

                 if (nv == 1) then ; CNT_TYP(:) = CNT_TYP1(n,:); AREA_TYP = AREA_TYP1(n,:,:); INTPUT_TID = INTPUT_TID1(n,:,:); endif
                 if (nv == 2) then ; CNT_TYP(:) = CNT_TYP2(n,:); AREA_TYP = AREA_TYP2(n,:,:); INTPUT_TID = INTPUT_TID2(n,:,:); endif
                 if (nv == 3) then ; CNT_TYP(:) = CNT_TYP3(n,:); AREA_TYP = AREA_TYP3(n,:,:); INTPUT_TID = INTPUT_TID3(n,:,:); endif
                 if (nv == 4) then ; CNT_TYP(:) = CNT_TYP4(n,:); AREA_TYP = AREA_TYP4(n,:,:); INTPUT_TID = INTPUT_TID4(n,:,:); endif
                 
                 ! Check INPUT PT1
                 call look4_input_cells (ityp_new, nplus, sub_tid, sub_ityp1, sub_fveg1, sub_area, &
                      t_count (n,nv), CNT_TYP(1), AREA_TYP(1,:), INTPUT_TID(1,:))
                 ! Check INPUT PT2
                 call look4_input_cells (ityp_new, nplus, sub_tid, sub_ityp2, sub_fveg2, sub_area, &
                      t_count (n,nv), CNT_TYP(2), AREA_TYP(2,:), INTPUT_TID(2,:))
                 ! Check INPUT ST1
                 call look4_input_cells (ityp_new, nplus, sub_tid, sub_ityp3, sub_fveg3, sub_area, &
                      t_count (n,nv), CNT_TYP(3), AREA_TYP(3,:), INTPUT_TID(3,:))
                 ! Check INPUT ST2
                 call look4_input_cells (ityp_new, nplus, sub_tid, sub_ityp4, sub_fveg4, sub_area, &
                      t_count (n,nv), CNT_TYP(4), AREA_TYP(4,:), INTPUT_TID(4,:))
                 
                 if (nv == 1) then ; CNT_TYP1(n,:) = CNT_TYP(:); AREA_TYP1(n,:,:) = AREA_TYP; INTPUT_TID1(n,:,:) = INTPUT_TID; endif
                 if (nv == 2) then ; CNT_TYP2(n,:) = CNT_TYP(:); AREA_TYP2(n,:,:) = AREA_TYP; INTPUT_TID2(n,:,:) = INTPUT_TID; endif
                 if (nv == 3) then ; CNT_TYP3(n,:) = CNT_TYP(:); AREA_TYP3(n,:,:) = AREA_TYP; INTPUT_TID3(n,:,:) = INTPUT_TID; endif
                 if (nv == 4) then ; CNT_TYP4(n,:) = CNT_TYP(:); AREA_TYP4(n,:,:) = AREA_TYP; INTPUT_TID4(n,:,:) = INTPUT_TID; endif

              endif SEEK
           END DO NVLOOP
                                   
           DEALLOCATE (sub_tid, sub_ityp1, sub_ityp2, sub_ityp3, sub_ityp4)
           DEALLOCATE (sub_fveg1, sub_fveg2, sub_fveg3, sub_fveg4, sub_area)
           
           all_found = .true.
           
           if((all_found).and.((CLMC_pf1(n) > fmin).and.(t_count(n,1) < 0))) all_found = .false.
           if((all_found).and.((CLMC_pf2(n) > fmin).and.(t_count(n,2) < 0))) all_found = .false.
           if((all_found).and.((CLMC_sf1(n) > fmin).and.(t_count(n,3) < 0))) all_found = .false.
           if((all_found).and.((CLMC_sf2(n) > fmin).and.(t_count(n,4) < 0))) all_found = .false.
           
           if(all_found) GO TO 100
           
           ! if not increase the window size
           nplus = 0
           
        END DO ZOOMOUT
        
100     continue !  if(mod (n,1000) == 0) print *, myid +1, n, Id_loc(local_id,:)
        !     if(mod (n,100) == 0) print *, myid +1, n, t_count (n,:)
        
     END DO TILES
     
     deallocate (min_lon, max_lon, min_lat, max_lat)
     deallocate (CLMC_pt1, CLMC_pt2, CLMC_st1, CLMC_st2, CNT_TYP, AREA_TYP, INTPUT_TID)
     deallocate (lonn, latt, tid_offl, mask, ityp_offl, fveg_offl, lonc, latc, tile_area)
     
     ! CREATE ANNUAL FILES (OR CLIM FILES ) DAILY DYNVEG DATA
     ! TOTALLY 24 FILES (LAI, SAI) x NUM_VEG (4) x NUM_ZONES (3)
     ! ---------------------------------------------------------
     
  endif DO_REGRID

  IF(ACTUAL_DV)  then

     Y_LOOP : DO YEAR = YEARB, YEARE 

        write (YYYY,'(i4.4)') YEAR
        ! Open the annual file and write header
        ! -------------------------------------

        IF (master_proc) THEN
           call create_output_files (file_units, YYYY, trim(EXPDIR))
           DO N = 1, NFILES; WRITE (file_units(N)) float((/YEAR-1,12,31,0,0,0,YEAR,1,1,0,0,0,NTILES,1/)); END DO
        ENDIF

        ! Read, process and gather CF180 LAI/SAI data
        ! -------------------------------------------

        write (YYYY,'(i4.4)') YEAR -1
        NCInID = 0
        STATUS = NF_OPEN (trim(GLDASName)//'.tavg1_tile_veg.'//YYYY//'1231_1200z.nc4', NF_NOWRITE, NCInID(1)) ; VERIFY_(STATUS) 

        if(donot_regrid) then
           if (master_proc) call READ_NOREGRID_PUT (NTILES, out_tile_data,NCInID, file_units)
           call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
        else
           DO NV = 1, num_veg
              if (nv == 1) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                   CNT_TYP1, AREA_TYP1, INTPUT_TID1, CLMC_pf1)
              if (nv == 2) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                   CNT_TYP2, AREA_TYP2, INTPUT_TID2, CLMC_pf2)
              if (nv == 3) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                   CNT_TYP3, AREA_TYP3, INTPUT_TID3, CLMC_sf1)
              if (nv == 4) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                   CNT_TYP4, AREA_TYP4, INTPUT_TID4, CLMC_sf2)
           END DO
        endif

        STATUS = NF_CLOSE (NCInID(1))
        write (YYYY,'(i4.4)') YEAR
        YR1 = YEAR

        M_LOOP : DO MONTH = 1,12

           MN1 = MONTH

           D_LOOP: DO DAY = 1, days_in_month(year, month)
      
              ! write the header
              ! ----------------

              DY1 = DAY + 1
              IF(DAY == days_in_month(year, month)) THEN; DY1  = 1; MN1  = MONTH + 1 ; ENDIF
              IF(MN1 == 13) then ; MN1 = 1; YR1 = YEAR + 1; ENDIF
              IF (master_proc) THEN; DO N = 1, NFILES; WRITE (file_units(N)) float((/YEAR,MONTH,DAY,0,0,0,YR1,MN1,DY1,0,0,0,NTILES,1/)); END DO; ENDIF

              ! Read, process and gather CF180 LAI/SAI data
              ! -------------------------------------------
              WRITE (YYYYMMDD,'(i4.4,i2.2,i2.2)') YEAR, MONTH, DAY
              NCInID = 0
              STATUS = NF_OPEN (trim(GLDASName)//'.tavg1_tile_veg.'//YYYYMMDD//'_1200z.nc4', NF_NOWRITE, NCInID(1)) ; VERIFY_(STATUS) 

              if(donot_regrid) then
                 if (master_proc) call READ_NOREGRID_PUT (NTILES, out_tile_data,NCInID, file_units)
                 call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
              else
              
                 DO NV = 1, num_veg
                    if (nv == 1) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                         CNT_TYP1, AREA_TYP1, INTPUT_TID1, CLMC_pf1)
                    if (nv == 2) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                         CNT_TYP2, AREA_TYP2, INTPUT_TID2, CLMC_pf2)
                    if (nv == 3) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                         CNT_TYP3, AREA_TYP3, INTPUT_TID3, CLMC_sf1)
                    if (nv == 4) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                         CNT_TYP4, AREA_TYP4, INTPUT_TID4, CLMC_sf2)
                 END DO
              endif

              STATUS = NF_CLOSE (NCInID(1))
              
           END DO D_LOOP
        END DO M_LOOP

        IF (master_proc) THEN; DO N = 1, NFILES; WRITE (file_units(N)) float((/YEAR+1,1,1,0,0,0,YEAR+1,1,2,0,0,0,NTILES,1/)); END DO; ENDIF

        write (YYYY,'(i4.4)') YEAR + 1
        NCInID = 0
        STATUS = NF_OPEN (trim(GLDASName)//'.tavg1_tile_veg.'//YYYY//'0101_1200z.nc4', NF_NOWRITE, NCInID(1)) ; VERIFY_(STATUS) 

        if(donot_regrid) then
           if (master_proc) call READ_NOREGRID_PUT (NTILES, out_tile_data,NCInID, file_units)
           call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
        else
           DO NV = 1, num_veg
              if (nv == 1) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                   CNT_TYP1, AREA_TYP1, INTPUT_TID1, CLMC_pf1)
              if (nv == 2) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                   CNT_TYP2, AREA_TYP2, INTPUT_TID2, CLMC_pf2)
              if (nv == 3) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                   CNT_TYP3, AREA_TYP3, INTPUT_TID3, CLMC_sf1)
              if (nv == 4) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                   CNT_TYP4, AREA_TYP4, INTPUT_TID4, CLMC_sf2)
           END DO
        endif
        
        STATUS = NF_CLOSE (NCInID(1))
        IF (master_proc) call create_output_files (file_units, YYYY, trim(EXPDIR), close_files = .true.)

     END DO Y_LOOP

  ELSE

     YYYY = 'CLIM'

     ! Open the annual file and write header
     ! -------------------------------------
     
     IF (master_proc) THEN
        call create_output_files (file_units, YYYY, trim(EXPDIR))
        DO N = 1, NFILES
           if (clim_sai) then
              if(N <= 12) then
                 WRITE (file_units(N)) float((/0,12,1,0,0,0,1,1,1,0,0,0,NTILES,1/))
              else
                 WRITE (file_units(N)) float((/0,12,31,0,0,0,1,1,1,0,0,0,NTILES,1/))
              endif
           else
              WRITE (file_units(N)) float((/0,12,31,0,0,0,1,1,1,0,0,0,NTILES,1/))
           endif
        END DO
     ENDIF

     ! Read, process and gather CF180 LAI/SAI data
     ! -------------------------------------------
     
     YYYY = 'YYYY'
     NCInID = 0
     STATUS = NF_OPEN (trim(GLDASName)//'.tavg1_tile_veg.'//YYYY//'1231_1200z.nc' , NF_NOWRITE, NCInID(1)) ; VERIFY_(STATUS) 
     if (clim_sai) STATUS = NF_OPEN (trim(GLDASName)//'.tavg1_tile_veg_monthly.'//YYYY//'12.nc4', NF_NOWRITE, NCInID(2)) ; VERIFY_(STATUS) 

     if(donot_regrid) then
        if (master_proc) call READ_NOREGRID_PUT (NTILES, out_tile_data,NCInID, file_units, 1)
        call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
     else
        DO NV = 1, num_veg
           if (nv == 1) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                CNT_TYP1, AREA_TYP1, INTPUT_TID1, CLMC_pf1,1)
           if (nv == 2) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                CNT_TYP2, AREA_TYP2, INTPUT_TID2, CLMC_pf2,1)
           if (nv == 3) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                CNT_TYP3, AREA_TYP3, INTPUT_TID3, CLMC_sf1,1)
           if (nv == 4) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                CNT_TYP4, AREA_TYP4, INTPUT_TID4, CLMC_sf2,1)
        END DO
     endif
     
     STATUS = NF_CLOSE (NCInID(1))
     if (clim_sai) STATUS = NF_CLOSE (NCInID(2))

     YR1  = 1
     YEAR = 1
     DO MONTH = 1,12

        MN1 = MONTH
        IF (master_proc) THEN       
           if (clim_sai) then
              MN1 = MONTH +1
              IF(MN1 == 13) then ; MN1 = 1; YR1 = 2; ENDIF
                 DO N = 1, NFILES
                    IF(N <= 12) WRITE (file_units(N)) float((/YEAR,MONTH,1,0,0,0,YR1,MN1,1,0,0,0,NTILES,1/))
                 END DO
                 YR1 = 1
              endif
           ENDIF

        MN1 = MONTH

        DO DAY = 1, days_in_month(1997, month)
      
           ! write the header
           ! ----------------

           DY1 = DAY + 1
           IF(DAY == days_in_month(1997, month)) THEN; DY1  = 1; MN1  = MONTH + 1 ; ENDIF
           IF(MN1 == 13) then ; MN1 = 1; YR1 = 2; ENDIF
           IF (master_proc) THEN              
              DO N = 1, NFILES
                 if (clim_sai) then
                    if(N > 12) WRITE (file_units(N)) float((/YEAR,MONTH,DAY,0,0,0,YR1,MN1,DY1,0,0,0,NTILES,1/))
                 else
                    WRITE (file_units(N)) float((/YEAR,MONTH,DAY,0,0,0,YR1,MN1,DY1,0,0,0,NTILES,1/))
                 endif
              END DO
           ENDIF
           
           ! Read, process and gather CF180 LAI/SAI data
           ! -------------------------------------------
           WRITE (YYYYMMDD,'(a4,i2.2,i2.2)') 'YYYY', MONTH, DAY
           WRITE (MM, '(i2.2)') MONTH
           NCInID = 0
           STATUS = NF_OPEN (trim(GLDASName)//'.tavg1_tile_veg.'//YYYYMMDD//'_1200z.nc', NF_NOWRITE, NCInID(1)) ; VERIFY_(STATUS) 
           if (clim_sai) STATUS = NF_OPEN (trim(GLDASName)//'.tavg1_tile_veg_monthly.'//YYYY//MM//'.nc4', NF_NOWRITE, NCInID(2)) ; VERIFY_(STATUS) 
 
           if(donot_regrid) then
              if (master_proc) then
                 if (day ==1) then
                    call READ_NOREGRID_PUT (NTILES, out_tile_data,NCInID, file_units,1)
                 else
                    call READ_NOREGRID_PUT (NTILES, out_tile_data,NCInID, file_units)
                 endif
              endif
              call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
           else
              DO NV = 1, num_veg
                 if (day ==1) then
                    if (nv == 1) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                         CNT_TYP1, AREA_TYP1, INTPUT_TID1, CLMC_pf1,1)
                    if (nv == 2) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                         CNT_TYP2, AREA_TYP2, INTPUT_TID2, CLMC_pf2,1)
                    if (nv == 3) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                         CNT_TYP3, AREA_TYP3, INTPUT_TID3, CLMC_sf1,1)
                    if (nv == 4) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                         CNT_TYP4, AREA_TYP4, INTPUT_TID4, CLMC_sf2,1)
                 else
                    if (nv == 1) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                         CNT_TYP1, AREA_TYP1, INTPUT_TID1, CLMC_pf1)
                    if (nv == 2) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                         CNT_TYP2, AREA_TYP2, INTPUT_TID2, CLMC_pf2)
                    if (nv == 3) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                         CNT_TYP3, AREA_TYP3, INTPUT_TID3, CLMC_sf1)
                    if (nv == 4) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                         CNT_TYP4, AREA_TYP4, INTPUT_TID4, CLMC_sf2)
                 endif
              END DO
           endif

           STATUS = NF_CLOSE (NCInID(1))
           if (clim_sai) STATUS = NF_CLOSE (NCInID(2))
        END DO 
        END DO 

        IF (master_proc) THEN
           DO N = 1, NFILES
              if (clim_sai) then
                 if(N <= 12) then
                    WRITE (file_units(N)) float((/2,1,1,0,0,0,2,2,1,0,0,0,NTILES,1/))
                 else
                    WRITE (file_units(N)) float((/2,1,1,0,0,0,2,1,2,0,0,0,NTILES,1/))
                 endif
              else
                 WRITE (file_units(N))  float((/2,1,1,0,0,0,2,1,2,0,0,0,NTILES,1/))
              endif
           END DO
        ENDIF
        NCInID = 0
        STATUS = NF_OPEN (trim(GLDASName)//'.tavg1_tile_veg.YYYY0101_1200z.nc', NF_NOWRITE, NCInID(1)) ; VERIFY_(STATUS) 
        if (clim_sai) STATUS = NF_OPEN (trim(GLDASName)//'.tavg1_tile_veg_monthly.'//YYYY//'01.nc4', NF_NOWRITE, NCInID(2)) ; VERIFY_(STATUS) 

        if(donot_regrid) then
           if (master_proc) call READ_NOREGRID_PUT (NTILES, out_tile_data,NCInID, file_units,1)
           call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
        else        
           DO NV = 1, num_veg
              if (nv == 1) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                   CNT_TYP1, AREA_TYP1, INTPUT_TID1, CLMC_pf1,1)
              if (nv == 2) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                   CNT_TYP2, AREA_TYP2, INTPUT_TID2, CLMC_pf2,1)
              if (nv == 3) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, & 
                   CNT_TYP3, AREA_TYP3, INTPUT_TID3, CLMC_sf1,1)
              if (nv == 4) call read_and_put (NTILES, loc_size, nt_local, low_ind, NCInID, file_units, NV, &
                   CNT_TYP4, AREA_TYP4, INTPUT_TID4, CLMC_sf2,1)
           END DO
        endif
        
        STATUS = NF_CLOSE (NCInID(1))
        if (clim_sai) STATUS = NF_CLOSE (NCInID(2))
        IF (master_proc) call create_output_files (file_units, YYYY, trim(EXPDIR), close_files = .true.)
  ENDIF

  call MPI_BARRIER( MPI_COMM_WORLD,STATUS)
  call MPI_Finalize(STATUS)

contains

  
  ! -----------------------------------------------------------------------------------------

  SUBROUTINE read_and_put (NTILES, NT, nt_local, low_ind, NCInID, file_units, NV, &
       CNT_TYP, AREA_TYP, INTPUT_TID, CLMC_fr, day)
    
    implicit none
    
    integer, intent (in) :: NTILES, NT, nt_local(numprocs), low_ind(numprocs), file_units (NFILES), &
         CNT_TYP (NT, NUM_VEG), INTPUT_TID(NT, NUM_VEG,Noff_per_pft), NV
    real   , intent (in) :: AREA_TYP (NT, NUM_VEG,Noff_per_pft), CLMC_fr (NT)
    INTEGER              :: TV, NZ, Var, file_no, STATUS
    CHARACTER*7          :: OUTVAR,   INVAR(4)
    REAL, allocatable, dimension (:,:) :: off_var
    integer, intent (in), optional     :: day
    integer, intent (in), dimension (:):: NCInID
    
    ALLOCATE (off_var   (1: NTILES_CN, 4))
    ! Must do         
    ! 'CNSAI[nv]1',  'CNSAI[nv]2',  'CNSAI[nv]3', &  
    ! 'CNLAI[nv]1',  'CNLAI[nv]2',  'CNLAI[nv]3', &    
    
    DO NZ = 1,NUM_ZONE
       DO Var = 1,2
          
          if (Var == 1) then
             write (OUTVAR  , '(a5,i1,i1)') 'CNSAI', nv,nz
             write (INVAR(1), '(a5,i1,i1)') 'CNSAI', 1,nz
             write (INVAR(2), '(a5,i1,i1)') 'CNSAI', 2,nz
             write (INVAR(3), '(a5,i1,i1)') 'CNSAI', 3,nz
             write (INVAR(4), '(a5,i1,i1)') 'CNSAI', 4,nz
          endif
          if (Var == 2) then
             write (OUTVAR  , '(a5,i1,i1)') 'CNLAI', nv,nz
             write (INVAR(1), '(a5,i1,i1)') 'CNLAI', 1,nz
             write (INVAR(2), '(a5,i1,i1)') 'CNLAI', 2,nz
             write (INVAR(3), '(a5,i1,i1)') 'CNLAI', 3,nz
             write (INVAR(4), '(a5,i1,i1)') 'CNLAI', 4,nz
          endif
              
          ! file_no = 6*nv -6 + 3*(var -1) + nz
          file_no = 12*(var -1) + 3*nv -3 + nz
          if (clim_sai) then  
             if (Var == 1) then
                if(present(day))then
                   DO TV = 1,NUM_VEG
                      STATUS = NF_GET_VARA_REAL (NCInID(2), VarID(NCInID(2),INVAR(TV)),  (/1,1/), (/NTILES_CN,1/), off_var(:,TV)) ; VERIFY_(STATUS)
                   END DO
                   call WRITE_OUTPUT (file_no, off_var, NTILES, NT, nt_local, low_ind,CNT_TYP, AREA_TYP, INTPUT_TID, CLMC_fr)
                endif
             else
                DO TV = 1,NUM_VEG
                   STATUS = NF_GET_VARA_REAL (NCInID(1), VarID(NCInID(1),INVAR(TV)),  (/1,1/), (/NTILES_CN,1/), off_var(:,TV)) ; VERIFY_(STATUS)
                END DO
                call WRITE_OUTPUT (file_no, off_var, NTILES, NT, nt_local, low_ind,CNT_TYP, AREA_TYP, INTPUT_TID, CLMC_fr)
             endif
          else
             
             ! Should check IN_VAR1[NZ] IN_VAR2[NZ] IN_VAR3[NZ] IN_VAR4[NZ]
             DO TV = 1,NUM_VEG
                STATUS = NF_GET_VARA_REAL (NCInID(1), VarID(NCInID(1),INVAR(TV)),  (/1,1/), (/NTILES_CN,1/), off_var(:,TV)) ; VERIFY_(STATUS)
             END DO
             call WRITE_OUTPUT (file_no, off_var, NTILES, NT, nt_local, low_ind,CNT_TYP, AREA_TYP, INTPUT_TID, CLMC_fr)
          endif
       END DO
    END DO
        
    deallocate (off_var)
  END SUBROUTINE read_and_put

  ! -------------------------------------------------------------------------
     
  SUBROUTINE WRITE_OUTPUT (FILE_UNIT, off_var, &
       NTILES, NT, nt_local, low_ind,CNT_TYP, AREA_TYP, INTPUT_TID, CLMC_fr)
      
      implicit none
      integer, intent (in) :: FILE_UNIT,  NTILES, NT, nt_local(numprocs), low_ind(numprocs), &
           CNT_TYP (NT, NUM_VEG), INTPUT_TID(NT, NUM_VEG,Noff_per_pft)  
      REAL,    intent (in), dimension (NTILES_CN, 4) :: off_var
      real   , intent (in) :: AREA_TYP (NT, NUM_VEG,Noff_per_pft), CLMC_fr (NT)
      REAL, allocatable, dimension (:)               :: out_local, out_global
      INTEGER              :: N, TV, NZ, I
      REAL                 :: sum_a
      
      IF (master_proc) ALLOCATE (out_global (1: NTILES))
      ALLOCATE (out_local (1: NT))
      
      out_local = undef
      
      LOCAL_LOOP : DO N = 1, NT
         if(CLMC_fr (n) > fmin) then
            ! this fraction is active
            sum_a         = 0.
            out_local (n) = 0.
            do TV = 1,  num_veg
               if(CNT_TYP (N, TV) > 0) then
                  ! loop through all available fractions of ssimilar out type
                  DO i = 1, CNT_TYP (N, TV) 
                     out_local (n) = out_local (n) + off_var(INTPUT_TID(n,tv,i),TV) * AREA_TYP (n,tv,i) 
                  END DO
                  sum_a = sum_a + SUM (AREA_TYP (n,tv,1:CNT_TYP (N, TV)))
               endif
            end do
            out_local (n) = out_local (n) / sum_a 
         endif
      END DO LOCAL_LOOP
      
      call MPI_Barrier(MPI_COMM_WORLD, STATUS)
      call MPI_GATHERV( out_local, nt_local(myid+1)   , MPI_real, &
           out_global, nt_local,low_ind-1, MPI_real, &
           0, MPI_COMM_WORLD, STATUS)       
      
      IF (master_proc) WRITE (file_unit) out_global
      IF (master_proc) deallocate (out_global)
      deallocate (out_local)
      
    END SUBROUTINE WRITE_OUTPUT
               
    ! -----------------------------------------------------------------------------------------
     
    SUBROUTINE READ_NOREGRID_PUT (NTILES, out_tile_data,NCInID, file_units, day)
       
      implicit none
      integer, intent (in) :: NTILES
      type(tile_data_type), INTENT (IN), DIMENSION(NTILES) :: out_tile_data 
      integer, intent (in) :: file_units (NFILES)
      integer, intent (in), dimension (:):: NCInID
      integer              :: n, STATUS
      integer, intent (in), optional     :: day
      REAL, allocatable, dimension (:)   :: out_global
      
      ALLOCATE (out_global (1: NTILES))
      do n = 1, NFILES       
         if (clim_sai) then       
            if(N <= 12) then
               if(present(day))then
                  STATUS = NF_GET_VARA_REAL (NCInID(2), VarID(NCInID(2),OUT_VARS(n)),  (/1,1/), (/NTILES,1/), out_global(:)) ; VERIFY_(STATUS)
                  WRITE (file_units (n)) out_global (out_tile_data%tile_id)
               endif
            else
               STATUS = NF_GET_VARA_REAL (NCInID(1), VarID(NCInID(1),OUT_VARS(n)),  (/1,1/), (/NTILES,1/), out_global(:)) ; VERIFY_(STATUS)    
               WRITE (file_units (n)) out_global (out_tile_data%tile_id)
            endif
         else
            STATUS = NF_GET_VARA_REAL (NCInID(1), VarID(NCInID(1),OUT_VARS(n)),  (/1,1/), (/NTILES,1/), out_global(:)) ; VERIFY_(STATUS)    
            WRITE (file_units (n)) out_global (out_tile_data%tile_id)
         endif
         
      end do
      
      STATUS = NF_CLOSE (NCInID(1))
      IF(NCInID(2) /= 0) STATUS = NF_CLOSE (NCInID(2))
      deallocate (out_global)
      
    END SUBROUTINE READ_NOREGRID_PUT

    ! -----------------------------------------------------------------------------------------

    SUBROUTINE look4_input_cells (OUT_TYPE, nplus, sub_tid, sub_ityp, sub_fveg, sub_area, &
                   t_count, CNT_TYP, AREA_TYP, INTPUT_TID)
      implicit none
      
      integer, intent (in)               :: OUT_TYPE, nplus
      integer, dimension(:), intent (in) :: sub_tid, sub_ityp
      real,    dimension(:), intent (in) :: sub_fveg, sub_area
      integer, intent (inout)            :: t_count, CNT_TYP
      integer, dimension(:), intent (inout) :: INTPUT_TID
      real,    dimension(:), intent (inout) :: AREA_TYP
      integer :: i
      
      WINDOW_LOOP : do i = 1, nplus
         
         if((sub_fveg (i) > fmin).and.(OUT_TYPE == sub_ityp(i)).and.(CNT_TYP < Noff_per_pft)) then
            if(t_count < 0)  t_count = 0 ! initialize if need be
            if(CNT_TYP < 0)  CNT_TYP = 0 ! initialize if need be
            t_count = t_count + 1
            CNT_TYP = CNT_TYP + 1
            INTPUT_TID(CNT_TYP) = sub_tid(i)
            AREA_TYP  (CNT_TYP) = sub_fveg(i) * sub_area (i)
            !          if( CNT_TYP == Noff_per_pft) then ; print *,'COUNT EXCEEDED Noff_per_pft' ; exit ; endif
         endif
         
      END DO WINDOW_LOOP
      
    END SUBROUTINE look4_input_cells

    ! -----------------------------------------------------------------------------------------
    
    SUBROUTINE read_OUT_TileData (NTILES, out_tile_data, MODELING_SYSTEM, BCSPATH, TILFILE, EXPDIR)
      
      implicit none
      
      CHARACTER(*), INTENT (IN) :: MODELING_SYSTEM, BCSPATH, TILFILE, EXPDIR
      INTEGER, INTENT (INOUT)   :: NTILES
      type(tile_data_type), pointer, INTENT (INOUT), DIMENSION(:) :: out_tile_data
      INTEGER, ALLOCATABLE, DIMENSION (:) :: IDUM
      REAL,    ALLOCATABLE, DIMENSION (:) :: RDUM
      REAL                                :: AREA
      INTEGER                             :: N, STATUS, NCID, PFAF
      CHARACTER*10                        :: TMPSTRING
      
      SYSTEM : IF (TRIM(MODELING_SYSTEM) == "GEOSldas") THEN
         
         OPEN (10, FILE = TRIM(TILFILE), FORM = 'UNFORMATTED', STATUS = 'OLD', ACTION = 'READ')
         
         read (10) NTILES
         
         ALLOCATE (OUT_TILE_DATA (1:NTILES))
         ALLOCATE (IDUM       (1:NTILES))
         ALLOCATE (RDUM       (1:NTILES))
         
         read (10) (out_tile_data(n)%tile_id,   n=1,NTILES)
         read (10) (IDUM(n),                 n=1,NTILES)
         read (10) (IDUM(n),                 n=1,NTILES)
         read (10) (out_tile_data(n)%com_lon,   n=1,NTILES)
         read (10) (out_tile_data(n)%com_lat,   n=1,NTILES)
         read (10) (out_tile_data(n)%min_lon,   n=1,NTILES)
         read (10) (out_tile_data(n)%max_lon,   n=1,NTILES)
         read (10) (out_tile_data(n)%min_lat,   n=1,NTILES)
         read (10) (out_tile_data(n)%max_lat,   n=1,NTILES)
         close (10, status = 'keep')
         
         STATUS = NF_OPEN (trim (EXPDIR)//'/input/restart/catchcn_internal_rst', NF_NOWRITE,NCID) ; VERIFY_(STATUS)
         
         DO N = 1, NUM_VEG
            STATUS = NF_GET_VARA_REAL (NCID, VarID(NCID,'ITY'),  (/1,N/), (/NTILES,1/), RDUM) ; VERIFY_(STATUS)
            out_tile_data(:)%ITY(N) = INT (RDUM)
            STATUS = NF_GET_VARA_REAL (NCID, VarID(NCID,'FVG'),  (/1,N/), (/NTILES,1/), out_tile_data(:)%FVG(N)) ; VERIFY_(STATUS) 
         END DO
         DEALLOCATE (IDUM, RDUM)
         STATUS = NF_CLOSE (NCID )
         
      ELSE
         
         OPEN (11, FILE = TRIM(BCSPATH)//'/clsm/catchment.def',      FORM = 'FORMATTED', STATUS = 'OLD', ACTION = 'READ') 
         
         read (11, *) NTILES
         ALLOCATE (OUT_TILE_DATA (1:NTILES))
         
         DO N = 1, NTILES
            read (11,*) out_tile_data(n)%tile_id, PFAF, out_tile_data(n)%min_lon, out_tile_data(n)%max_lon, out_tile_data(n)%min_lat, out_tile_data(n)%max_lat
         END DO
         
         close (11, status = 'keep')
         
         call ReadCNTilFile  (TRIM(BCSPATH)//'/'//trim(TILFILE), NTILES, out_tile_data(:)%com_lon, out_tile_data(:)%com_lat)
         call ReadCLMvegFile (TRIM(BCSPATH)//'/clsm/CLM_veg_typs_fracs', NTILES,                                   &
              out_tile_data(:)%ITY(1),  out_tile_data(:)%ITY(2), out_tile_data(:)%ITY(3), out_tile_data(:)%ITY(4), &
              out_tile_data(:)%FVG(1),  out_tile_data(:)%FVG(2), out_tile_data(:)%FVG(3), out_tile_data(:)%FVG(4))          
         
      ENDIF SYSTEM
      
    END SUBROUTINE read_OUT_TileData

    ! ------------------------------------------------------------------

    integer function days_in_month(year, month)
      
      ! return the number of days in a given month
      
      implicit none
      
      integer :: year, month
      
      integer, dimension(12), parameter :: days_in_month_leap = &
           (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
      
      integer, dimension(12), parameter :: days_in_month_nonleap = &
           (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
      
      if (is_leap_year(year)) then
         days_in_month = days_in_month_leap(month) 
      else
         days_in_month = days_in_month_nonleap(month) 
      end if
      
    end function days_in_month
    
  ! ------------------------------------------------------------------
    
    logical function is_leap_year(year)
      
      implicit none
      
      integer :: year
      
      if (mod(year,4) /= 0) then 
         is_leap_year = .false.
      else if (mod(year,400) == 0) then
         is_leap_year = .true.
      else if (mod(year,100) == 0) then 
         is_leap_year = .false.
      else
         is_leap_year = .true.
      end if
      
    end function is_leap_year

  ! ------------------------------------------------------------------
    
    integer function VarID (NCFID, VNAME) 
      
      integer, intent (in)      :: NCFID
      character(*), intent (in) :: VNAME
      integer                   :: status
      
      STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,VarID)
      IF (STATUS .NE. NF_NOERR) &
           CALL HANDLE_ERR(STATUS, trim(VNAME))  
      
    end function VarID
    
  ! -----------------------------------------------------------------------
    
  SUBROUTINE HANDLE_ERR(STATUS, Line)
    
    INTEGER,      INTENT (IN) :: STATUS
    CHARACTER(*), INTENT (IN) :: Line
    
    IF (STATUS .NE. NF_NOERR) THEN
       PRINT *, trim(Line),': ',STATUS, NF_STRERROR(STATUS)
       STOP 'Stopped'
    ENDIF
    
  END SUBROUTINE HANDLE_ERR

  ! *****************************************************************************
  
  subroutine init_MPI()
    
    ! initialize MPI
    
    call MPI_INIT(STATUS)
    
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, STATUS )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, STATUS )

    if (myid .ne. 0)  master_proc = .false.
    
!    write (*,*) "MPI process ", myid, " of ", numprocs, " is alive"    
!    write (*,*) "MPI process ", myid, ": master_proc=", master_proc

  end subroutine init_MPI

  ! ***************************************************************************** 

   subroutine ReadCNTilFile (InCNTileFile, nt, xlon, xlat, tile_area)
     
     implicit none
     character(*), intent (in) ::  InCNTileFile
     integer , intent (in) :: nt
     real, dimension (nt), intent(inout) :: xlon, xlat
     real, dimension (nt), intent(out), optional   :: tile_area
     integer :: n,icnt,ityp, ik,jk, ierr
     real    :: xval,yval, pf, fr
     real (kind =8) :: dxy, d2r, lats

     dxy = 360._8/2880.     
     d2r = PI/180._8
    
     open(11,file=InCNTileFile, &
          form='formatted',action='read',status='old')
     
     do n = 1,8 ! skip header
        read(11,*)
     end do
   
     icnt = 0
     ityp = 100
     ierr = 0

     do while (ierr == 0) ! loop over land tiles
        read(11,*,IOSTAT=ierr) ityp,pf,xval,yval,ik,jk, fr
        if((ierr ==0).and.(ityp == 100)) then
           icnt = icnt + 1
           xlon(icnt) = xval
           xlat(icnt) = yval

           if (present (tile_area)) then
              lats = DBLE (yval)
              tile_area (icnt) = fr * REAL((sin(d2r*(lats+0.5*dxy)) -sin(d2r*(lats-0.5*dxy)))*(dxy*d2r))
           endif
        endif
     end do
     
     close(11)
    

   end subroutine ReadCNTilFile

  ! ***************************************************************************** 

   subroutine ReadCLMvegFile (VegFile, nt, vt1,vt2,vt3,vt4,vf1,vf2,vf3,vf4)  
     
     implicit none
     character(*), intent (in) ::  VegFile
     integer , intent (in) :: nt
     real    , dimension (nt), intent(inout) :: vf1,vf2,vf3,vf4
     real    , dimension (nt), intent(inout) :: vt1,vt2,vt3,vt4
     integer :: n,tid, pfaf

     
     open(11,file=VegFile, &
          form='formatted',action='read',status='old')
     
     DO N = 1, NT
        read (11,*) TID, PFAF,vt1(n),vt2(n),vt3(n),vt4(n), & 
             vf1(n),vf2(n),vf3(n),vf4(n)
        vf1(n) =  vf1(n)/100.
        vf2(n) =  vf2(n)/100.
        vf3(n) =  vf3(n)/100.
        vf4(n) =  vf4(n)/100.
     END DO
     
     close(11)
     

   end subroutine ReadCLMvegFile

 ! ***************************************************************************** 

   SUBROUTINE create_output_files (file_units, YYYY, EXPDIR, close_files)

     implicit none

     integer, dimension (:), intent (inout) :: file_units
     character(*), intent (in)              :: YYYY, EXPDIR
     logical, intent(in), optional          :: close_files
     integer                                :: n

     if(present (close_files)) then        

        do n = 1, NFILES
           close (file_units (n), status = 'keep')
        end do       

     else  
     
        do n = 1, NFILES
          file_units (n) = 10 + n 
          open (file_units ( n), file = trim(EXPDIR)//'/VEGDATA/'//OUT_VARS(n)//'_'//YYYY//'.data', form = 'unformatted', action = 'write', status='new')
        end do

     endif
     

   END SUBROUTINE create_output_files

END PROGRAM regrid_DYNVEG
