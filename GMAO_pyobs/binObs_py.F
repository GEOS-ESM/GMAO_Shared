      subroutine binObs2D(lon,lat,obs,nobs,gObs,im,jm,missing)
      implicit NONE
      integer, intent(in) :: im, jm, nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm)
!
! Bins obs. It assumes a global 2D GEOS-5 A-Grid:
!
!  Longitudes in [-180,10)
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, n
      real :: dLon, dLat, xLon
      real :: nsample(im,jm)

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gObs = 0.0
      nSample = 0.0
      do n = 1, nobs

         if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

         xLon = lon(n)
         if ( xLon >= 180. ) xLon = xLon - 360.
 
         i = 1 + nint((xlon + 180.  ) / dLon)
         j = 1 + nint((lat(n) + 90. ) / dLat)

         if ( i>im ) i = i - im
         if ( i<1  ) i = i + im
         
!        if (abs(obs(n)-missing) > 0.1) then
         if (abs(obs(n)-missing) > 0.01*abs(missing)) then
         gObs(i,j) = gObs(i,j) + obs(n)
         nSample(i,j) = nSample(i,j) + 1.0
         end if
      end do

      where ( nsample > 0 ) 
              gObs = gObs / nSample
      elsewhere
              gObs = missing
      end where

      end

!-----------

      subroutine binObs3D(lon,lat,obs,nobs,gObs,im,jm,km,missing)
      implicit NONE
      integer, intent(in) :: im, jm, km,nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs,km)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm,km)
!
! Bins obs. It assumes a global 3D GEOS-5 A-Grid:
!
!  Longitudes in [-180,10]
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, k, n
      real :: dLon, dLat, xLon
      real :: nSample(im,jm)

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gObs = 0.0

      do k = 1, km
         nSample = 0.0
         
         do n = 1,nobs
                 
            if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

            xLon = lon(n)
            if ( xLon >= 180. ) xLon = xLon - 360.
            
            i = 1 + nint((xlon + 180.  ) / dLon)
            j = 1 + nint((lat(n) + 90. ) / dLat)

            if ( i>im ) i = i - im
            if ( i<1  ) i = i + im
       
!           if (abs(obs(n,k)-missing) > 0.1) then
            if (abs(obs(n,k)-missing) > 0.01*abs(missing)) then
               gObs(i,j,k) = gObs(i,j,k) + obs(n,k)
               nSample(i,j) = nSample(i,j) + 1.0
            end if

         end do

         where ( nsample > 0 ) 
            gObs(:,:,k) = gObs(:,:,k) / nSample
         elsewhere
            gObs(:,:,k) = missing
         end where
      
      end do

      end
!-------------
      subroutine binRms2D(lon,lat,obs,nobs,gObs,im,jm,missing)
      implicit NONE
      integer, intent(in) :: im, jm, nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm)
!
! Bins RMS.  It assumes a global GEOS-5 A-Grid:
!
!  Longitudes in [-180,10)
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, n
      real :: dLon, dLat, xLon
      real :: nsample(im,jm)

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gObs = 0.0
      nSample = 0.0
      do n = 1, nobs

         if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

         xLon = lon(n)
         if ( xLon >= 180. ) xLon = xLon - 360.
 
         i = 1 + nint((xlon + 180.  ) / dLon)
         j = 1 + nint((lat(n) + 90. ) / dLat)

         if ( i>im ) i = i - im
         if ( i<1  ) i = i + im
         
!        if (abs(obs(n)-missing) > 0.1) then
         if (abs(obs(n)-missing) > 0.01*abs(missing)) then
         gObs(i,j) = gObs(i,j) + obs(n)*obs(n)
         nSample(i,j) = nSample(i,j) + 1.0
         end if
      end do

      where ( nsample > 0 ) 
              gObs = sqrt(abs(gObs / nSample))
      elsewhere
              gObs = missing
      end where

      end

!-----------

      subroutine binRms3D(lon,lat,obs,nobs,gObs,im,jm,km,missing)
      implicit NONE
      integer, intent(in) :: im, jm, km,nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs,km)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm,km)
!
! Bins RMS. It assumes a global 3D GEOS-5 A-Grid:
!
!  Longitudes in [-180,10]
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, k, n
      real :: dLon, dLat, xLon
      real :: nSample(im,jm)

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gObs = 0.0

      do k = 1, km
         nSample = 0.0
         
         do n = 1,nobs
                 
            if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

            xLon = lon(n)
            if ( xLon >= 180. ) xLon = xLon - 360.
            
            i = 1 + nint((xlon + 180.  ) / dLon)
            j = 1 + nint((lat(n) + 90. ) / dLat)

            if ( i>im ) i = i - im
            if ( i<1  ) i = i + im
       
!           if (abs(obs(n,k)-missing) > 0.1) then
            if (abs(obs(n,k)-missing) > 0.01*abs(missing)) then
               gObs(i,j,k) = gObs(i,j,k) + obs(n,k)*obs(n,k)
               nSample(i,j) = nSample(i,j) + 1.0
            end if

         end do

         where ( nsample > 0 ) 
            gObs(:,:,k) = sqrt(abs(gObs(:,:,k) / nSample))
         elsewhere
            gObs(:,:,k) = missing
         end where
      
      end do

      end
!--------------
      subroutine binObs3DP(lon,lat,obs,pobs,pe,nobs,nlev,
     & gObs,im,jm,km,missing)
      implicit NONE
      integer, intent(in) :: im, jm, km,nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs,nlev)
      real, intent(in)  :: pobs(nobs,nlev)
      real, intent(in)  :: pe(im,jm,km)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm,km-1)
!



!
! Bins obs. It assumes a global 3D GEOS-5 A-Grid:
!
!  Longitudes in [-180,10]
!  Latitudes  in [-90,90]
!
! 
      integer :: i, j, k, n, P, nlev, m
      real :: dLon, dLat, xLon
      real :: nSample(im,jm,km)
      
      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      
      gObs = 0.0
      m = 0
      do k = 1, km-1
         nSample = 0.0
         do n = 1,nobs
                    
            if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

            xLon = lon(n)
            if ( xLon >= 180. ) xLon = xLon - 360.
            
            i = 1 + nint((xlon + 180.  ) / dLon)
            j = 1 + nint((lat(n) + 90. ) / dLat)

            if ( i>im ) i = i - im
            if ( i<1  ) i = i + im
            
            do P = 1 , nlev

               if (k == km-1) then    ! if CALIPSO pressure > GEOS-5 P for edge 73 -> 
                                      ! put CALIPSO values in level 72 
                  if (pobs(n,P) > pe(i,j,k)) then
                   
                    if (abs(obs(n,P)-missing) > 0.01*abs(missing)) then
                                            
                      gObs(i,j,k) = gObs(i,j,k) + obs(n,P)
                      nSample(i,j,k) = nSample(i,j,k) + 1.0
                    end if
                  end if
               
               else if ( (pobs(n,P) > pe(i,j,k)).AND.
     &            (pobs(n,P) < pe(i,j,k+1) ) )then
                                      
                   if (abs(obs(n,P)-missing) > 0.01*abs(missing)) then

                      gObs(i,j,k) = gObs(i,j,k) + obs(n,P)
                      nSample(i,j,k) = nSample(i,j,k) + 1.0
                   end if
               end if
               
            end do
         end do
         
         where ( nsample(:,:,k) > 0 )            
            gObs(:,:,k) = gObs(:,:,k) / nSample(:,:,k)
         elsewhere
            gObs(:,:,k) = missing
         end where
      
      end do

      end

!--------------
      subroutine binAreas(lon,lat,area,nobs,gArea,im,jm)
      implicit NONE
      integer, intent(in) :: im, jm, nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: area(nobs)
      real, intent(out) :: gArea(im,jm)
!
! Bins areas of no-fire pixels. It assumes a global GEOS-5 A-Grid:
!
!  Longitudes in [-180,180)
!  Latitudes  in [-90,90]
!
!
      integer :: i, j, n
      real :: dLon, dLat, xLon

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gArea = 0.0
      do n = 1, nobs

         xLon = lon(n)
         if ( xLon >= 180. ) xLon = xLon - 360.
 
         i = 1 + nint((xlon + 180.  ) / dLon)
         j = 1 + nint((lat(n) + 90. ) / dLat)

         if ( i>im ) i = i - im
         if ( i<1  ) i = i + im

         gArea(i,j) = gArea(i,j) + area(n)

      end do

      end

!--------------
      subroutine binAreasNR(lon,lat,area,nobs,gArea,im,jm)
      implicit NONE
      integer, intent(in) :: im, jm, nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: area(nobs)
      real, intent(out) :: gArea(im,jm)
!
! Bins areas of no-fire pixels. It assumes a uniform DEPE grid.
!
!  Longitudes in (-180,180)
!  Latitudes  in (-90,90)
!
!
      integer :: i, j, n
      real :: dLon, dLat, xLon

      dLon = 360. / im
      dLat = 180. / jm

      gArea = 0.0
      do n = 1, nobs

         xLon = lon(n)
         if ( xLon >= 180. ) xLon = xLon - 360.
 
         i = nint((xLon   + 180.0) / dLon + 0.5)
         j = nint((lat(n) +  90.0) / dLat + 0.5)

         if ( i>im ) then
             i = i - im
!            print *, '[w] adjustment i>im', xLon, lat(n)
         end if
         if ( i<1  ) then
             i = i + im
!            print *, '[w] adjustment i<1', xLon, lat(n)
         end if

         if ( j>jm ) then
             j = jm
!            print *, '[w] adjustment j>jm', xLon, lat(n)
         end if
         if ( j<1  ) then
             j = 1
!            print *, '[w] adjustment j<1', xLon, lat(n)
         end if

         gArea(i,j) = gArea(i,j) + area(n)

      end do

      end
!..................................................................

      subroutine findBox(lon,lat,nobs,im,jm,gridbox)
      implicit NONE
      integer, intent(in)  :: im, jm, nobs
      real, intent(in)     :: lon(nobs)
      real, intent(in)     :: lat(nobs)
      integer, intent(out) :: gridbox(nobs)
!
!     Returns gridbox number for each lon,lat on input.
!
      print *, 'findBox: not implemented yet.'
      gridbox = 0

      end

!..................................................................

      subroutine decimateSwath(a, im, jm, b, ih, jh, idel, missing, rc )
      implicit NONE
      integer, intent(in) :: im     
      integer, intent(in) :: jm     
      integer, intent(in) :: ih     
      integer, intent(in) :: jh     
      integer, intent(in) :: idel   
                                    

      real, intent(in)  :: b(ih,jh) 
      real, intent(in)  :: missing  
      
      real, intent(out) :: a(im,jm) 
      integer, intent(out) :: rc    
      
!
!     Averages high res array *b* producing low res array *a*.
!

!                           ---

      integer :: idel_, jdel_, i, j, i_, j_, i1, i2, j1, j2, n
      real*8  :: accum

#define VALID(x) (abs(x-missing) > 0.01*missing)

      rc = 0

!     Consistency check
!     -----------------
      idel_ = nint(float(ih)/float(im))
      jdel_ = nint(float(jh)/float(jm))
      if ( idel_ /= idel .OR. jdel_ /= idel ) then
         rc = 1
         return
      end if
      if ( mod(idel,2) == 0 ) then
         rc = 2
         return
      end if
      
!     Average *good* values of high res data
!     --------------------------------------
      do j = 1, jm

       j1 = (j-1)*idel + 1
       j2 = min(jh,j1+idel-1)

       do i = 1, im

          i1 = (i-1)*idel + 1
          i2 = min(ih,i1+idel-1)

          accum = 0.0
          n = 0
          do j_ = j1, j2
             do i_ = i1, i2
                if ( VALID(b(i_,j_)) ) then
                   accum = accum + b(i_,j_)
                   n = n + 1
                end if
             end do
          end do

          if ( n == 0 ) then
             a(i,j) = missing
          else
             a(i,j) = accum / n
          end if

       end do
      end do

      return
      end 
