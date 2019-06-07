  subroutine gauss2grid(gfld,glon,glat,grad,glat2,rfld,rlon,rlat)

    use mod_param
    use m_const, only: undef
    implicit none 

    integer,              intent(in)  :: glon,glat,rlon,rlat,glat2
    real(kind=kind_evod), intent(in)  :: gfld(glon,glat),grad(glat)  
    real(kind=kind_evod), intent(out) :: rfld(rlon,rlat)  
    real(kind=kind_evod), allocatable :: lat2(:),dphi(:),dlam(:),lons(:),lats(:)
    real(kind=kind_evod)              :: pi,dl,dp,snp,ssp,lon,lat
    integer                           :: loc 
    integer                           :: i,j,k 


    pi = atan(1.0)*4.0         !constant

    ! create lats array with pole points

    allocate(lat2(glat+2))
    allocate(dphi(glat+2))
    allocate(dlam(glon))

    lat2(1) = -90.0
    do j=2,glat2+1
      lat2(glat2+3-j)= -grad(j-1)*180.0/pi 
      lat2(j+glat2)  =  grad(j-1)*180.0/pi
    enddo 
    lat2(glat+2) = 90.0

    ! compute dlam and dphi array
    dl = 2.0*pi/glon
    dlam(:) = dl

    do j=1,glat+1
      dphi(j) = ( lat2(j+1)-lat2(j) )*pi/180.0
    enddo
    dphi(glat+2) = undef

    ! create pole values 
    rfld = 0.
    snp=0.
    ssp=0.
    do i=1,glon
      snp=snp+gfld(i,glat)    
      ssp=ssp+gfld(i,1)
    enddo
    snp = snp/glon
    ssp = ssp/glon
    do i=1,rlon
     rfld(i,rlat) = snp
     rfld(i,1)    = ssp
    enddo

    ! create output lons and lats and interpolate

    dl = 2*pi/rlon
    dp = pi/(rlat-1)
    allocate(lons(rlon*rlat))
    allocate(lats(rlon*rlat))

    loc =0 
    do j=1,rlat
      do i=1,rlon
        loc=loc+1
        lon = -pi + (i-1)*dl
        lons(loc) = lon
      enddo
    enddo 

    loc =0 
    do j=1,rlat
      lat = -pi/2.0 + (j-1)*dp
      do i=1,rlon
        loc=loc+1
        lats(loc) = lat
      enddo
    enddo 

    call interp_h(gfld,glon,glat,1,dlam,dphi,0.0,90.0,0.0,&
                  rfld,rlon*rlat,lons,lats,1,3,.false.,undef)  

    deallocate(lons,lats)
    deallocate(lat2,dphi,dlam)

  end subroutine gauss2grid
