module sst_ice_helpers

! !USES:
  USE netcdf

implicit none
private

public bin2bin
public fill_Land
public hflip
public interp_to_eight_deg
public read_Reynolds
public read_OSTIA
public read_Ostia_quart
public read_input
public read_input_quart
public write_bin
public check

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine is SAME as Larry's bin2bin.F90 [.../pre/OSTIA/bin2bin.F90] and comes 
! from https://github.com/GEOS-ESM/GMAO_Shared/tree/main/GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/
!
      SUBROUTINE bin2bin ( qi,imi,jmi,qo,imo,jmo,undef )
!***********************************************************************
!
!  PURPOSE:
!  ========
!    Bin an input field, qi(imi,jmi), to an output array qo(imo,jmo)
!
!  INPUT:
!  ======
!    qi ......... Input array(imi,jmi)
!
!  OUTPUT:
!  =======
!    qo ......... Output array(imo,jmo)
!
!  NOTES:
!  ======
!    Input and Output arrays are assumed Dateline_Edge and Pole_Edge.
!             Each box is referenced by the latitude and longitude of
!             its southwest corner, not its center point.  Thus,
!             the value associated with a coordinate actually
!             represents the value centered to the northeast of that point.
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************

      IMPLICIT NONE

      INTEGER,    INTENT(IN)   :: imi, jmi
      INTEGER,    INTENT(IN)   :: imo, jmo

      REAL,       INTENT(IN)   :: undef
      REAL,       INTENT(IN)   :: qi(imi,jmi)
      REAL,       INTENT(OUT)  :: qo(imo,jmo)

! Local variables
      INTEGER                  :: i,j,ibeg,iend,jbeg,jend
      INTEGER                  :: ii,jj,itmp
      real*8   dlam, dphi
      real*8    sum1,sum2
      real*8    zlat,zlon
      real*8    lon1,lon2,wx
      real*8    lat1,lat2,wy
      real*8    lonbeg,lonend,lat,coslat
      real*8    latbeg,latend
      real*8    pi,dz
      real*8    lon_out(imo)
      real*8    lat_out(jmo)

      pi   = 4*DATAN(1.0D0)
      dz   =   pi/jmi
      dlam = 2*pi/imo
      dphi =   pi/jmo

! Compute Output Lons and Lats
! ----------------------------
      lon_out(1) = -pi+0.5D0*dlam
      DO i = 2, imo
        lon_out(i) = lon_out(i-1) + dlam
      ENDDO

      lat_out(1) = -pi*0.5D0+0.5D0*dphi
      DO j = 2, jmo-1
        lat_out(j) = lat_out(j-1) + dphi
      ENDDO
      lat_out(jmo) =  pi*0.5D0-0.5D0*dphi

! Bin Input Array to Output Array
! -------------------------------
      DO j=1,jmo
         DO i=1,imo

           zlat = lat_out(j); zlon = lon_out(i)

           latbeg = zlat-dphi*0.5D0; latend = zlat+dphi*0.5D0
           lonbeg = zlon-dlam*0.5D0; lonend = zlon+dlam*0.5D0

            ibeg = 0.5D0+(lonbeg+pi)      /dz
            iend = 0.5D0+(lonend+pi)      /dz
            jbeg = 0.5D0+(latbeg+pi*0.5D0)/dz
            jend = 0.5D0+(latend+pi*0.5D0)/dz

! Check for Begin and End Errors due to Truncation
! ------------------------------------------------
      lon2 = -pi +  ibeg   *dz
      lon1 = -pi + (iend-1)*dz
      if( lon2.lt.lonbeg ) ibeg = ibeg + 1
      if( lon1.gt.lonend ) ibeg = ibeg - 1
      if( ibeg.gt.imi ) then
          print *, 'ERROR after Truncation Check!'
          print *, 'ibeg: ',ibeg
          stop
      endif
      if( ibeg.lt.0 ) then
          print *, 'ERROR after Truncation Check!'
          print *, 'ibeg: ',ibeg
          stop
      endif

      lat2 = -pi*0.5D0 +  jbeg   *dz
      lat1 = -pi*0.5D0 + (jend-1)*dz
      if( lat2.lt.latbeg ) jbeg = jbeg + 1
      if( lat1.gt.latend ) jbeg = jbeg - 1

            IF( jbeg.lt.0 .or. jend.gt.jmi ) THEN
              PRINT *, 'Bounding jbeg jend values Error for (i,j): ',i,j
              PRINT *, 'jbeg: ',jbeg,' latbeg: ',latbeg*180/pi
              PRINT *, 'jend: ',jend,' latend: ',latend*180/pi
              STOP
            ENDIF
            IF( jbeg.lt.1 ) jbeg = 1

            sum1 = 0.0D0 ; sum2 = 0.0D0

            DO jj=jbeg,jend
             lat2 = -pi*0.5D0 +  jj       *dz
             lat  = -pi*0.5D0 + (jj-0.5D0)*dz
             lat1 = -pi*0.5D0 + (jj-1)    *dz

             coslat = max( min(dcos(lat),1.0D0),0.0D0 )

             if( lat2.lt.latbeg .or. lat1.gt.latend ) then
                 print *
                 print *, ' ERROR!'
                 print *
                 print *, '     j: ',j
                 print *
                 print *, 'latend: ',latend*180/pi
                 print *, '  zlat: ',zlat  *180/pi
                 print *, 'latbeg: ',latbeg*180/pi
                 print *
                 print *, '    jj: ',jj
                 print *, '  jend: ',jend,' r_jend: ',0.5D0+(latend+pi*0.5D0)/dz
                 print *, '  jbeg: ',jbeg,' r_jbeg: ',0.5D0+(latbeg+pi*0.5D0)/dz
                 print *
                 print *, 'lat2  : ',lat2  *180/pi
                 print *, 'lat   : ',lat   *180/pi
                 print *, 'lat1  : ',lat1  *180/pi
                 print *
                 print *, 'lat2-latbeg: ',lat2-latbeg
                 print *, 'lat1-latend: ',lat1-latend
                 stop
             endif

                                  wy = 1.0D0
             IF( lat1.lt.latbeg ) wy = (lat2-latbeg)/dz
             IF( lat2.gt.latend ) wy = (latend-lat1)/dz

             IF( iend.gt.imi ) THEN
               PRINT *, 'Bounding iend Error for (i,j): ',i,j
               PRINT *, 'jbeg: ',jbeg,' latbeg: ',latbeg*180/pi
               PRINT *, 'jend: ',jend,' latend: ',latend*180/pi
               PRINT *, 'iend: ',iend,' lonend: ',lonend*180/pi
               PRINT *, 'iend: ',iend,' r_iend: ',0.5D0+(lonend+pi)/dz
               STOP
             ENDIF

             IF(ibeg.ge.1) THEN
               DO ii=ibeg,iend
                 lon1 = -pi  + (ii-1)*dz; lon2 = -pi  +  ii   *dz
                                      wx = 1.0D0
                 IF( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
                 IF( lon2.gt.lonend ) wx = (lonend-lon1)/dz
                 IF( wx.lt.0.0 .or. wy.lt.0.0 ) THEN
                     PRINT *
                     PRINT *, '     i: ',i
                     PRINT *
                     PRINT *, 'lonend: ',lonend*180/pi
                     PRINT *, '  zlon: ',zlon  *180/pi
                     PRINT *, 'lonbeg: ',lonbeg*180/pi
                     PRINT *
                     PRINT *, '    ii: ',ii
                     PRINT *, '  iend: ',iend,' r_iend: ',0.5D0+(latend+pi)/dz
                     PRINT *, '  ibeg: ',ibeg,' r_ibeg: ',0.5D0+(latbeg+pi)/dz
                     PRINT *
                     PRINT *, '  lon2: ',lon2  *180/pi
                     PRINT *, '  lon1: ',lon1  *180/pi
                     PRINT *, '    wx: ',wx
                     PRINT *
                     PRINT *, '     j: ',j
                     PRINT *
                     PRINT *, 'latend: ',latend*180/pi
                     PRINT *, '  zlat: ',zlat  *180/pi
                     PRINT *, 'latbeg: ',latbeg*180/pi
                     PRINT *
                     PRINT *, '    jj: ',jj
                     PRINT *, '  jend: ',jend,' r_jend: ',0.5D0+(latend+pi*0.5D0)/dz
                     PRINT *, '  jbeg: ',jbeg,' r_jbeg: ',0.5D0+(latbeg+pi*0.5D0)/dz
                     PRINT *
                     PRINT *, 'lat2  : ',lat2  *180/pi
                     PRINT *, 'lat   : ',lat   *180/pi
                     PRINT *, 'lat1  : ',lat1  *180/pi
                     PRINT *, '    wy: ',wy
                     STOP
                 ENDIF
                 if( qi(ii,jj).ne.undef ) then
                     sum1 = sum1 + qi(ii,jj)*coslat*wx*wy
                     sum2 = sum2 +           coslat*wx*wy
                 endif
               ENDDO
             ELSE
                 itmp = 0.5D0+(lonbeg+0.1D0*dz+3*pi)/dz
                 IF( itmp.gt.imi ) THEN
                   PRINT *, 'Bounding itmp Error for (i,j): ',i,j
                   PRINT *, 'jbeg: ',jbeg,' latbeg: ',latbeg*180/pi
                   PRINT *, 'jend: ',jend,' latend: ',latend*180/pi
                   PRINT *, 'itmp: ',itmp,' lontmp: ',lonbeg*180/pi
                   STOP
                 ENDIF
               DO ii=itmp,imi
                 lon1 = -pi  + (ii-1)*dz; lon2 = -pi  +  ii   *dz
                                           wx = 1.0D0
                 IF( lon1.lt.lonbeg+2*pi ) wx = (lon2-lonbeg-2*pi)/dz
                 IF( lon2.gt.lonend+2*pi ) wx = (2*pi+lonend-lon1)/dz
                 IF( wx.lt.0.0 .or. wy.lt.0.0 ) THEN
                     PRINT *
                     PRINT *, '     i: ',i
                     PRINT *
                     PRINT *, 'lonend: ',lonend*180/pi
                     PRINT *, '  zlon: ',zlon  *180/pi
                     PRINT *, 'lonbeg: ',lonbeg*180/pi
                     PRINT *
                     PRINT *, '    ii: ',ii
                     PRINT *, '  iend: ',iend,' r_iend: ',0.5D0+(latend+pi)/dz
                     PRINT *, '  ibeg: ',ibeg,' r_ibeg: ',0.5D0+(latbeg+pi)/dz
                     PRINT *
                     PRINT *, '  lon2: ',lon2  *180/pi
                     PRINT *, '  lon1: ',lon1  *180/pi
                     PRINT *, '    wx: ',wx
                     PRINT *
                     PRINT *, '     j: ',j
                     PRINT *
                     PRINT *, 'latend: ',latend*180/pi
                     PRINT *, '  zlat: ',zlat  *180/pi
                     PRINT *, 'latbeg: ',latbeg*180/pi
                     PRINT *
                     PRINT *, '    jj: ',jj
                     PRINT *, '  jend: ',jend,' r_jend: ',0.5D0+(latend+pi*0.5D0)/dz
                     PRINT *, '  jbeg: ',jbeg,' r_jbeg: ',0.5D0+(latbeg+pi*0.5D0)/dz
                     PRINT *
                     PRINT *, 'lat2  : ',lat2  *180/pi
                     PRINT *, 'lat   : ',lat   *180/pi
                     PRINT *, 'lat1  : ',lat1  *180/pi
                     PRINT *, '    wy: ',wy
                     STOP
                 ENDIF
                 if( qi(ii,jj).ne.undef ) then
                     sum1 = sum1 + qi(ii,jj)*coslat*wx*wy
                     sum2 = sum2 +           coslat*wx*wy
                 endif
               ENDDO
               DO ii=1,iend
                  lon1 = -pi  + (ii-1)*dz; lon2 = -pi  +  ii   *dz
                                       wx = 1.0D0
                  IF( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
                  IF( lon2.gt.lonend ) wx = (lonend-lon1)/dz
                  IF( wx.lt.0.0 .or. wy.lt.0.0 ) THEN
                      PRINT *
                      PRINT *, '     i: ',i
                      PRINT *
                      PRINT *, 'lonend: ',lonend*180/pi
                      PRINT *, '  zlon: ',zlon  *180/pi
                      PRINT *, 'lonbeg: ',lonbeg*180/pi
                      PRINT *
                      PRINT *, '    ii: ',ii
                      PRINT *, '  iend: ',iend,' r_iend: ',0.5D0+(latend+pi)/dz
                      PRINT *, '  ibeg: ',ibeg,' r_ibeg: ',0.5D0+(latbeg+pi)/dz
                      PRINT *
                      PRINT *, '  lon2: ',lon2  *180/pi
                      PRINT *, '  lon1: ',lon1  *180/pi
                      PRINT *, '    wx: ',wx
                      PRINT *
                      PRINT *, '     j: ',j
                      PRINT *
                      PRINT *, 'latend: ',latend*180/pi
                      PRINT *, '  zlat: ',zlat  *180/pi
                      PRINT *, 'latbeg: ',latbeg*180/pi
                      PRINT *
                      PRINT *, '    jj: ',jj
                      PRINT *, '  jend: ',jend,' r_jend: ',0.5D0+(latend+pi*0.5D0)/dz
                      PRINT *, '  jbeg: ',jbeg,' r_jbeg: ',0.5D0+(latbeg+pi*0.5D0)/dz
                      PRINT *
                      PRINT *, 'lat2  : ',lat2  *180/pi
                      PRINT *, 'lat   : ',lat   *180/pi
                      PRINT *, 'lat1  : ',lat1  *180/pi
                      PRINT *, '    wy: ',wy
                      STOP
                  ENDIF
                  if( qi(ii,jj).ne.undef ) then
                      sum1 = sum1 + qi(ii,jj)*coslat*wx*wy
                      sum2 = sum2 +           coslat*wx*wy
                  endif
               ENDDO
             ENDIF ! IF(ibeg.ge.1)
            ENDDO   ! DO jj=jbeg,jend

            if( sum2.ne.0.0D0 ) then
                qo(i,j) = sum1/sum2
            else
                qo(i,j) = undef
            endif

         ENDDO
      ENDDO

      RETURN
!--------------------------------------------------------------------------------
      END SUBROUTINE bin2bin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! From https://github.com/GEOS-ESM/GMAO_Shared/tree/main/GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/
! Purpose is to read a netcdf file from NOAA 0.25-deg OISST from 
! https://www.ncei.noaa.gov/products/optimum-interpolation-sst
!
      SUBROUTINE read_Reynolds(ncFileName, NLAT, NLON, LAT, LON, &
                               SST, ICE, myUNDEF)
!---------------------------------------------------------------------------
          IMPLICIT NONE

          CHARACTER (LEN = *),             INTENT(IN)    :: ncFileName
          INTEGER,                         INTENT(IN)    :: NLAT, NLON
          REAL,                            INTENT(IN)    :: myUNDEF
          REAL, DIMENSION(NLAT),           INTENT(OUT)   :: LAT
          REAL, DIMENSION(NLON),           INTENT(OUT)   :: LON
          REAL, DIMENSION(NLON,NLAT),      INTENT(OUT)   :: SST
          REAL, DIMENSION(NLON,NLAT),      INTENT(OUT)   :: ICE

! GET TO KNOW THESE BY ncdump -h
          REAL, PARAMETER :: sst_FillValue       = -999
          REAL, PARAMETER :: sst_offset          =  273.15              ! we need sst in K, Reynolds has it in deg C
          REAL, PARAMETER :: sst_scale_factor    =  0.01
          REAL, PARAMETER :: ice_FillValue       = -999
          REAL, PARAMETER :: ice_offset          =  0.0
          REAL, PARAMETER :: ice_scale_factor    =  0.01
          
! netCDF ID for the file and data variable.
          INTEGER :: ncid, varid1, varid2, varid3, varid4
!---------------------------------------------------------------------------

! Open the file.
          CALL check( nf90_open(ncFileName, nf90_nowrite, ncid))

! Get the varid of the data variable, based on its name.
! Get lat
          CALL check( nf90_inq_varid(ncid, "lat", varid1))
          CALL check( nf90_get_var(ncid, varid1,  LAT))
! Get lon
          CALL check( nf90_inq_varid(ncid, "lon", varid2))
          CALL check( nf90_get_var(ncid, varid2,  LON))
! Get SST
          CALL check( nf90_inq_varid(ncid, "sst", varid3))
          CALL check( nf90_get_var(ncid, varid3, SST))
! Get Ice Concentration
          CALL check( nf90_inq_varid(ncid, "ice", varid4))
          CALL check( nf90_get_var(ncid, varid4, ICE))
! Close nc file.
          CALL check( nf90_close(ncid))
! .....................................................................
! Use scale factor & offset
          WHERE( SST /= sst_FillValue)
              SST = sst_scale_factor * SST + sst_offset
          ENDWHERE
          WHERE( ICE /= ice_FillValue)
              ICE = ice_scale_factor * ICE + ice_offset
          ENDWHERE
! Unify undef
          WHERE( SST == sst_FillValue)
              SST = myUNDEF
          ENDWHERE
          WHERE( ICE == ice_FillValue)
              ICE = myUNDEF
          ENDWHERE
!---------------------------------------------------------------------------

      END SUBROUTINE read_Reynolds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! From https://github.com/GEOS-ESM/GMAO_Shared/tree/main/GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/
! Purpose is to replace `undef` values in input q(im,jm) with _nearest_ interpolated values.
! Output array: q(im,jm) has no undefined values.
      SUBROUTINE fill_Land (q,im,jm,undef)
!---------------------------------------------------------------------------
      implicit none

      integer  im,jm
      real     undef
      real     q(im,jm)
      integer  i,j,k,L,i0,nundef
      real     qz(im)
      real     dist,dq
!---------------------------------------------------------------------------

      do j=1,jm
         qz = q(:,j)
         nundef = count( qz.eq.undef )
         if( nundef.eq.im .or. nundef.eq.0 ) cycle

         do i0=1,im
         if( q(i0,j).ne.undef ) exit
         enddo

         do k=i0,im+i0-1
            L=k
            if(L.gt.im) L=L-im
            qz(k-i0+1) = q(L,j)
         enddo

         do i=2,im
             if( qz(i).ne.undef ) cycle
             do k=i+1,im
                if( qz(k).eq.undef ) cycle
                dist = k-i+1
                dq = ( qz(k)-qz(i-1) )/dist
                exit
             enddo
             if( k.eq.im+1) then
                dist = k-i+1
                dq = ( qz(1)-qz(i-1) )/dist
             endif
             do L=i,k-1
                qz(L) = qz(i-1) + (L-i+1)*dq
             enddo
         enddo

         do k=i0,im+i0-1
            L=k
            if(L.gt.im) L=L-im
            q(L,j) = qz(k-i0+1)
         enddo

      enddo ! j=1,jm

      return
!---------------------------------------------------------------------------
      END SUBROUTINE fill_Land
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! From https://github.com/GEOS-ESM/GMAO_Shared/tree/main/GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/
! flip lon: [0, 360] to [-180, 180]
      SUBROUTINE hflip( q, im, jm)
!---------------------------------------------------------------------------
          IMPLICIT NONE
          INTEGER  im,jm,i,j
          REAL   q(im,jm), dum(im)

            DO j=1,jm
              DO i=1,im/2
                 dum(i)      = q(i+im/2,j)
                 dum(i+im/2) = q(i,j)
              ENDDO
              q(:,j) = dum(:)
            ENDDO
!---------------------------------------------------------------------------
      END SUBROUTINE hflip
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! From https://github.com/GEOS-ESM/GMAO_Shared/tree/main/GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/
! This subroutine is SAME as Larry's interp_reynolds.F90 [.../pre/OSTIA/interp_reynolds.F90]
! Purpose is to interpolate input qi(imi,jmi) at a regular 0.25-deg to 0.125-deg qo(imo,jmo).
!
      subroutine interp_to_eight_deg (qi,imi,jmi,qo,imo,jmo,undef)
      implicit none

      integer,intent(in) :: imi,jmi
      integer,intent(in) :: imo,jmo
      real,   intent(in) :: undef
      real,   intent(in) :: qi(imi,jmi)
      real,   intent(out):: qo(imo,jmo)

      integer ib,jb
      integer ii1,ii2,ii3,ii4
      integer ji1,ji2,ji3,ji4
      integer io1,io2,io3,io4
      integer jo1,jo2,jo3,jo4

      real qz(imi,0:jmi+1)

      if( imi.ne.2*jmi .or. &
          imo.ne.2*jmo .or. &
          imo.ne.2*imi .or. &
          jmo.ne.2*jmi ) then
          PRINT *
          PRINT *, 'ERROR!       Output Resolution: ',imo,jmo
          PRINT *, 'must be twice Input Resolution: ',imi,jmi
          PRINT *
          stop 1
      endif

      qo = undef

      qz(:,0    ) = undef
      qz(:,1:jmi) = qi
      qz(:,jmi+1) = undef

      do jb=0,jmi
      do ib=1,imi

      ii1 = ib ; ii2 = ib+1 ; ii3 = ib+1 ; ii4 = ib
      ji1 = jb ; ji2 = jb   ; ji3 = jb+1 ; ji4 = jb+1

      io1 = 2*ib ; io2 = 2*ib+1 ; io3 = 2*ib+1 ; io4 = 2*ib
      jo1 = 2*jb ; jo2 = 2*jb   ; jo3 = 2*jb+1 ; jo4 = 2*jb+1

      if( ii2.gt.imi ) ii2 = ii2-imi
      if( ii3.gt.imi ) ii3 = ii3-imi
      if( io2.gt.imo ) io2 = io2-imo
      if( io3.gt.imo ) io3 = io3-imo

! Case (a)
! --------
      if( qz(ii1,ji1).ne.undef   .and. &
          qz(ii2,ji2).eq.undef   .and. &
          qz(ii3,ji3).eq.undef   .and. &
          qz(ii4,ji4).eq.undef ) then
          qo(io1,jo1) = qz(ii1,ji1)
      endif

! Case (b)
! --------
      if( qz(ii1,ji1).eq.undef   .and. &
          qz(ii2,ji2).ne.undef   .and. &
          qz(ii3,ji3).eq.undef   .and. &
          qz(ii4,ji4).eq.undef ) then
          qo(io2,jo2) = qz(ii2,ji2)
      endif

! Case (c)
! --------
      if( qz(ii1,ji1).eq.undef   .and. &
          qz(ii2,ji2).eq.undef   .and. &
          qz(ii3,ji3).ne.undef   .and. &
          qz(ii4,ji4).eq.undef ) then
          qo(io3,jo3) = qz(ii3,ji3)
      endif

! Case (d)
! --------
      if( qz(ii1,ji1).eq.undef   .and. &
          qz(ii2,ji2).eq.undef   .and. &
          qz(ii3,ji3).eq.undef   .and. &
          qz(ii4,ji4).ne.undef ) then
          qo(io4,jo4) = qz(ii4,ji4)
      endif

! Case (e)
! --------
      if( qz(ii1,ji1).ne.undef   .and. &
          qz(ii2,ji2).ne.undef   .and. &
          qz(ii3,ji3).eq.undef   .and. &
          qz(ii4,ji4).eq.undef ) then
          qo(io1,jo1) = 0.75*qz(ii1,ji1) + 0.25*qz(ii2,ji2)
          qo(io2,jo2) = 0.25*qz(ii1,ji1) + 0.75*qz(ii2,ji2)
      endif

! Case (f)
! --------
      if( qz(ii1,ji1).eq.undef   .and. &
          qz(ii2,ji2).ne.undef   .and. &
          qz(ii3,ji3).ne.undef   .and. &
          qz(ii4,ji4).eq.undef ) then
          qo(io2,jo2) = 0.75*qz(ii2,ji2) + 0.25*qz(ii3,ji3)
          qo(io3,jo3) = 0.25*qz(ii2,ji2) + 0.75*qz(ii3,ji3)
      endif

! Case (g)
! --------
      if( qz(ii1,ji1).eq.undef   .and. &
          qz(ii2,ji2).eq.undef   .and. &
          qz(ii3,ji3).ne.undef   .and. &
          qz(ii4,ji4).ne.undef ) then
          qo(io3,jo3) = 0.75*qz(ii3,ji3) + 0.25*qz(ii4,ji4)
          qo(io4,jo4) = 0.25*qz(ii3,ji3) + 0.75*qz(ii4,ji4)
      endif

! Case (h)
! --------
      if( qz(ii1,ji1).ne.undef   .and. &
          qz(ii2,ji2).eq.undef   .and. &
          qz(ii3,ji3).eq.undef   .and. &
          qz(ii4,ji4).ne.undef ) then
          qo(io1,jo1) = 0.75*qz(ii1,ji1) + 0.25*qz(ii4,ji4)
          qo(io4,jo4) = 0.25*qz(ii1,ji1) + 0.75*qz(ii4,ji4)
      endif

! Case (i)
! --------
      if( qz(ii1,ji1).ne.undef   .and. &
          qz(ii2,ji2).eq.undef   .and. &
          qz(ii3,ji3).ne.undef   .and. &
          qz(ii4,ji4).eq.undef ) then
          qo(io1,jo1) = 0.75*qz(ii1,ji1) + 0.25*qz(ii3,ji3)
          qo(io3,jo3) = 0.25*qz(ii1,ji1) + 0.75*qz(ii3,ji3)
      endif

! Case (j)
! --------
      if( qz(ii1,ji1).eq.undef   .and. &
          qz(ii2,ji2).ne.undef   .and. &
          qz(ii3,ji3).eq.undef   .and. &
          qz(ii4,ji4).ne.undef ) then
          qo(io2,jo2) = 0.75*qz(ii2,ji2) + 0.25*qz(ii4,ji4)
          qo(io4,jo4) = 0.25*qz(ii2,ji2) + 0.75*qz(ii4,ji4)
      endif

! Case (k)
! --------
      if( qz(ii1,ji1).ne.undef   .and. &
          qz(ii2,ji2).ne.undef   .and. &
          qz(ii3,ji3).ne.undef   .and. &
          qz(ii4,ji4).eq.undef ) then
          qo(io1,jo1) = 0.75*qz(ii1,ji1) + 0.125*(qz(ii2,ji2)+qz(ii3,ji3))
          qo(io3,jo3) = 0.75*qz(ii3,ji3) + 0.125*(qz(ii2,ji2)+qz(ii1,ji1))
          qo(io2,jo2) = 0.50*qz(ii2,ji2) + 0.250*(qz(ii1,ji1)+qz(ii3,ji3))
      endif

! Case (L)
! --------
      if( qz(ii1,ji1).eq.undef   .and. &
          qz(ii2,ji2).ne.undef   .and. &
          qz(ii3,ji3).ne.undef   .and. &
          qz(ii4,ji4).ne.undef ) then
          qo(io2,jo2) = 0.75*qz(ii2,ji2) + 0.125*(qz(ii4,ji4)+qz(ii3,ji3))
          qo(io4,jo4) = 0.75*qz(ii4,ji4) + 0.125*(qz(ii2,ji2)+qz(ii3,ji3))
          qo(io3,jo3) = 0.50*qz(ii3,ji3) + 0.250*(qz(ii2,ji2)+qz(ii4,ji4))
      endif

! Case (m)
! --------
      if( qz(ii1,ji1).ne.undef   .and. &
          qz(ii2,ji2).eq.undef   .and. &
          qz(ii3,ji3).ne.undef   .and. &
          qz(ii4,ji4).ne.undef ) then
          qo(io1,jo1) = 0.75*qz(ii1,ji1) + 0.125*(qz(ii4,ji4)+qz(ii3,ji3))
          qo(io3,jo3) = 0.75*qz(ii3,ji3) + 0.125*(qz(ii1,ji1)+qz(ii4,ji4))
          qo(io4,jo4) = 0.50*qz(ii4,ji4) + 0.250*(qz(ii1,ji1)+qz(ii3,ji3))
      endif

! Case (n)
! --------
      if( qz(ii1,ji1).ne.undef   .and. &
          qz(ii2,ji2).ne.undef   .and. &
          qz(ii3,ji3).eq.undef   .and. &
          qz(ii4,ji4).ne.undef ) then
          qo(io4,jo4) = 0.75*qz(ii4,ji4) + 0.125*(qz(ii1,ji1)+qz(ii2,ji2))
          qo(io2,jo2) = 0.75*qz(ii2,ji2) + 0.125*(qz(ii1,ji1)+qz(ii4,ji4))
          qo(io1,jo1) = 0.50*qz(ii1,ji1) + 0.250*(qz(ii2,ji2)+qz(ii4,ji4))
      endif

! Case (o)
! --------
      if( qz(ii1,ji1).ne.undef   .and. &
          qz(ii2,ji2).ne.undef   .and. &
          qz(ii3,ji3).ne.undef   .and. &
          qz(ii4,ji4).ne.undef ) then
          qo(io1,jo1) = ( 9*qz(ii1,ji1) + 3*qz(ii2,ji2) + 3*qz(ii4,ji4) + 1*qz(ii3,ji3) )/16.0
          qo(io2,jo2) = ( 9*qz(ii2,ji2) + 3*qz(ii1,ji1) + 3*qz(ii3,ji3) + 1*qz(ii4,ji4) )/16.0
          qo(io3,jo3) = ( 9*qz(ii3,ji3) + 3*qz(ii2,ji2) + 3*qz(ii4,ji4) + 1*qz(ii1,ji1) )/16.0
          qo(io4,jo4) = ( 9*qz(ii4,ji4) + 3*qz(ii1,ji1) + 3*qz(ii3,ji3) + 1*qz(ii2,ji2) )/16.0
      endif

      enddo
      enddo

      return
      end SUBROUTINE interp_to_eight_deg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! From https://github.com/GEOS-ESM/GMAO_Shared/tree/main/GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/
! Purpose is to read a netcdf file from UK MetOffice OSTIA SST system
! https://ghrsst-pp.metoffice.gov.uk/ostia-website/index.html
!
      SUBROUTINE read_OSTIA (ncFileName, VARNAME, NLAT, NLON, VAR, myUNDEF)
!---------------------------------------------------------------------------
          IMPLICIT NONE

          CHARACTER (LEN = *),               INTENT(IN)    :: ncFileName
          CHARACTER (LEN = *),               INTENT(IN)    :: VARNAME
          INTEGER,                           INTENT(IN)    :: NLAT, NLON
          REAL,                              INTENT(IN)    :: myUNDEF
          REAL, DIMENSION(NLON,NLAT),        INTENT(OUT)   :: VAR

! GET TO KNOW THESE BY ncdump -h
          REAL, PARAMETER :: sst_FillValue       = -32768
          REAL, PARAMETER :: sst_offset          =  273.15
          REAL, PARAMETER :: sst_scale_factor    =  0.01
          REAL, PARAMETER :: ice_FillValue       = -128
          REAL, PARAMETER :: ice_offset          =  0.0
          REAL, PARAMETER :: ice_scale_factor    =  0.01
          
! netCDF ID for the file and data variable.
          INTEGER :: ncid, varid3
          REAL    :: FillValue
          REAL    :: offset
          REAL    :: scale_factor
!---------------------------------------------------------------------------
! Open the file.
          CALL check( nf90_open(ncFileName, nf90_nowrite, ncid))

! Get lat
!         CALL check( nf90_inq_varid(ncid, "lat", varid1))
!         CALL check( nf90_get_var(ncid, varid1,  LAT))

! Get lon
!         CALL check( nf90_inq_varid(ncid, "lon", varid2))
!         CALL check( nf90_get_var(ncid, varid2,  LON))

! Get VAR
          CALL check( nf90_inq_varid(ncid, trim(VARNAME), varid3))
          CALL check( nf90_get_var(ncid, varid3, VAR))

! Close nc file.
          CALL check( nf90_close(ncid))
! .....................................................................
! Use scale factor & offset

          if( trim(VARNAME)=="analysed_sst" ) then
              FillValue    = sst_FillValue
              offset       = sst_offset
              scale_factor = sst_scale_factor
          endif
          if( trim(VARNAME)=="sea_ice_fraction" ) then
              FillValue    = ice_FillValue
              offset       = ice_offset
              scale_factor = ice_scale_factor
          endif

          WHERE( VAR /= FillValue)
              VAR = scale_factor * VAR + offset
          ENDWHERE

! Unify undef
          WHERE( VAR == FillValue)
                 VAR = myUNDEF
          ENDWHERE
!---------------------------------------------------------------------------
      END SUBROUTINE read_OSTIA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! From https://github.com/GEOS-ESM/GMAO_Shared/tree/main/GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/
! Purpose is to read a netcdf file that has been generated by grads/lats4d
      SUBROUTINE read_Ostia_quart (ncFileName, VARNAME, NLAT, NLON, VAR, myUNDEF)
!---------------------------------------------------------------------------
          IMPLICIT NONE

          CHARACTER (LEN = *),               INTENT(IN)    :: ncFileName
          CHARACTER (LEN = *),               INTENT(IN)    :: VARNAME
          INTEGER,                           INTENT(IN)    :: NLAT, NLON
          REAL,                              INTENT(IN)    :: myUNDEF
          REAL, DIMENSION(NLON,NLAT),        INTENT(OUT)   :: VAR

! GET TO KNOW THESE BY ncdump -h
          REAL, PARAMETER :: sst_FillValue       = -9.99e+33
          REAL, PARAMETER :: ice_FillValue       = -9.99e+33
          
! netCDF ID for the file and data variable.
          INTEGER :: ncid, varid3
          REAL    :: FillValue
!---------------------------------------------------------------------------
! Open the file.
          CALL check( nf90_open(ncFileName, nf90_nowrite, ncid))

! Get lat
!         CALL check( nf90_inq_varid(ncid, "lat", varid1))
!         CALL check( nf90_get_var(ncid, varid1,  LAT))

! Get lon
!         CALL check( nf90_inq_varid(ncid, "lon", varid2))
!         CALL check( nf90_get_var(ncid, varid2,  LON))

! Get VAR
          CALL check( nf90_inq_varid(ncid, trim(VARNAME), varid3))
          CALL check( nf90_get_var(ncid, varid3, VAR))

! Close nc file.
          CALL check( nf90_close(ncid))
! .....................................................................
! Use scale factor & offset

          if( trim(VARNAME)=="analysed_sst" ) then
              FillValue    = sst_FillValue
!             offset       = sst_offset
!             scale_factor = sst_scale_factor
          endif
          if( trim(VARNAME)=="sea_ice_fractio" ) then
              FillValue    = ice_FillValue
!             offset       = ice_offset
!             scale_factor = ice_scale_factor
          endif

!         WHERE( VAR /= FillValue)
!             VAR = scale_factor * VAR + offset
!         ENDWHERE

! Unify undef
          WHERE( VAR == FillValue)
                 VAR = myUNDEF
          ENDWHERE
!---------------------------------------------------------------------------
      END SUBROUTINE read_Ostia_quart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       NOTES:
!              1. REYNOLDS file name format changes-- cannot be hardcoded!
!      .......................................................................
!
      SUBROUTINE read_input(inputFile, iDebug, today, tomrw, fileName, NLAT, NLON, &
                            iAdjust_SST_SIC, SST_Thr, iERR)
!---------------------------------------------------------------------------
          IMPLICIT NONE

          CHARACTER (LEN = *),  INTENT(IN)    :: inputFile
          INTEGER,              INTENT(IN)    :: iDebug
          CHARACTER (LEN = *),  INTENT(OUT)   :: today, tomrw
          CHARACTER (LEN = *),  INTENT(OUT)   :: fileName(2)
          INTEGER,              INTENT(OUT)   :: NLAT, NLON
          INTEGER,              INTENT(OUT)   :: iERR
          INTEGER,              INTENT(OUT)   :: iAdjust_SST_SIC
          REAL,                 INTENT(OUT)   :: SST_Thr
!---------------------------------------------------------------------------

!       READ *, inputFileName
        OPEN (UNIT = 21, FILE = inputFile, STATUS = 'old')

!       Read multi-line input
        READ (21, '(A)') today
        READ (21, '(A)') tomrw
        READ (21, '(A)') fileName(1)                                 ! Reynolds file
        READ (21, '(A)') fileName(2)                                 ! OSTIA    file
        READ (21, '(I5)') NLAT
        READ (21, '(I5)') NLON
        READ (21, '(I5)') iAdjust_SST_SIC                             ! adjust SIC based on SST?
        READ (21, '(F5.2)') SST_Thr                                   ! adjust SIC based on what value of SST?
        CLOSE(21)
!      .......................................................................
!      CHECK USER INPUT. Die if not correct
!      All other checks must be done here.
!      .......................................................................
        iERR = 0
        IF( today == tomrw)  THEN
          iERR = 1
          PRINT *, 'Processing Start date: ', today
          PRINT *, 'is SAME as End date:   ', tomrw
          PRINT *, 'End date must be AFTER Start date'
        END IF
!      .......................................................................

        IF( iDebug /= 0 ) THEN
          PRINT *, '---------------------------------------'
          PRINT *, 'From read_input: '
          PRINT *, 'Today:         ', today
          PRINT *, 'Tomorrow:      ', tomrw
          PRINT *, 'Reynolds file: ', fileName(1)
          PRINT *, 'OSTIA    file: ', fileName(2)
          PRINT *, 'NLAT & NLON:   ', NLAT, NLON
          PRINT *, 'iAdjust_SST_SIC:', iAdjust_SST_SIC
          PRINT *, 'SST_Thr:        ', SST_Thr
          PRINT *, '---------------------------------------'
        END IF
!---------------------------------------------------------------------------
      END SUBROUTINE read_input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       NOTES:
!              1. REYNOLDS & NSIDC file name format changes-- cannot be hardcoded!
!      .......................................................................
!
      SUBROUTINE read_input_quart(inputFile, iDebug, today, tomrw, & 
                                  fileName, NLAT, NLON,            &
                                  iERR, max_diff_SST, max_diff_ICE)
!---------------------------------------------------------------------------
          IMPLICIT NONE

          CHARACTER (LEN = *),  INTENT(IN)    :: inputFile
          INTEGER,              INTENT(IN)    :: iDebug
          CHARACTER (LEN = *),  INTENT(OUT)   :: today, tomrw
          CHARACTER (LEN = *),  INTENT(OUT)   :: fileName(2)
          INTEGER,              INTENT(OUT)   :: NLAT, NLON
          INTEGER,              INTENT(OUT)   :: iERR

          REAL,                 INTENT(OUT)   :: max_diff_SST
          REAL,                 INTENT(OUT)   :: max_diff_ICE
!---------------------------------------------------------------------------

!       READ *, inputFileName
        OPEN (UNIT = 21, FILE = inputFile, STATUS = 'old')

!       Read multi-line input
        READ (21, '(A)') today
        READ (21, '(A)') tomrw
        READ (21, '(A)') fileName(1)             ! Reynolds file
        READ (21, '(A)') fileName(2)             ! OSTIA    file
        READ (21, '(I4)') NLAT
        READ (21, '(I4)') NLON
        READ (21, '(F8.4)') max_diff_SST            ! max allowed diff in SST between Reynolds and OSTIA
        READ (21, '(F8.4)') max_diff_ICE            ! max allowed diff in SIC between Reynolds and OSTIA
        CLOSE(21)
!      .......................................................................
!      CHECK USER INPUT. Die if not correct
!      All other checks must be done here.
!      .......................................................................
        iERR = 0
        IF( today == tomrw)  THEN
          iERR = 1
          PRINT *, 'Processing Start date: ', today
          PRINT *, 'is SAME as End date:   ', tomrw
          PRINT *, 'End date must be AFTER Start date'
        END IF
!      .......................................................................

       IF ( (max_diff_SST < 0.) .or. (max_diff_SST > 2.)) THEN
          PRINT *, 'Value of max_diff_SST is set to ', max_diff_SST
          PRINT *, 'Reynolds and OSTIA SSTs is out of bounds'
       END IF
       IF ( (max_diff_ICE < 1.e-6) .or. (max_diff_ICE > 1.)) THEN
          PRINT *, 'Value of max_diff_ICE is set to ', max_diff_ICE
          PRINT *, 'Reynolds and OSTIA Ice concentrations is out of bounds'
       END IF
!      .......................................................................
        IF( iDebug /= 0 ) THEN
          PRINT *, '---------------------------------------'
          PRINT *, 'From read_input: '
          PRINT *, 'Today:         ', today
          PRINT *, 'Tomorrow:      ', tomrw
          PRINT *, 'Reynolds file: ', fileName(1)
          PRINT *, 'OSTIA    file: ', fileName(2)
          PRINT *, 'NLAT & NLON:   ', NLAT, NLON
          PRINT *, '---------------------------------------'
        END IF
!---------------------------------------------------------------------------
      END SUBROUTINE read_input_quart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check(status)
        implicit none
        INTEGER, INTENT (IN) :: status
        IF (status /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(status))
         PRINT *, "ERROR in reading NetCDF file for SST and sea ice concentration."
         STOP
        END IF
      END SUBROUTINE check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Write out a binary file with `HEADER` and data, see below for its format.
!
    SUBROUTINE write_bin( today_year, today_mon, today_day, &
                          tomrw_year, tomrw_mon, tomrw_day, &
                          today, nlat, nlon, sst, ice)
!---------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN)     :: today_year, today_mon, today_day, &
                               tomrw_year, tomrw_mon, tomrw_day, &
                               nlat, nlon
    CHARACTER (LEN=*),INTENT(IN) :: today

    REAL, INTENT(IN)        :: sst(nlon, nlat), &
                               ice(nlon, nlat)

    REAL                    :: HEADER(14)
    CHARACTER (LEN = 40)    :: fileName_sst, fileName_ice
!---------------------------------------------------------------------------

    HEADER(1)    = REAL(today_year); HEADER(7)     = REAL(tomrw_year)
    HEADER(2)    = REAL(today_mon);  HEADER(8)     = REAL(tomrw_mon)
    HEADER(3)    = REAL(today_day);  HEADER(9)     = REAL(tomrw_day)
    HEADER(4:6)  = 0.0;              HEADER(10:12) = 0.0                  ! hours, min, sec are assumed to 00:00:00
    HEADER(13)   = REAL(nlon);       HEADER(14)    = REAL(nlat)

!---------------------------------------------------------------------------
!       Write out SST, ICE concentration data
!
    fileName_sst = 'Ostia_sst_' // today //'.bin'
    fileName_ice = 'Ostia_ice_' // today //'.bin'

    OPEN (UNIT = 991, FILE = fileName_sst, FORM = 'unformatted', STATUS = 'new')
    OPEN (UNIT = 992, FILE = fileName_ice, FORM = 'unformatted', STATUS = 'new')

    WRITE(991) HEADER
    WRITE(992) HEADER
    WRITE(991) sst
    WRITE(992) ice
    CLOSE(991)
    CLOSE(992)
!---------------------------------------------------------------------------

    END SUBROUTINE write_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sst_ice_helpers
