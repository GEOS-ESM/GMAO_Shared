subroutine g2g(datain,iin,jin,iout,jout,frac,dataout,undef,imin,jmin,imout,jmout,km,ntiles)
  implicit none

  integer, intent(in) :: imin,jmin,imout,jmout,km,ntiles
  real,    intent(in) :: datain(km,jmin,imin), frac(ntiles), undef 
  integer, intent(in) :: iin(ntiles), jin(ntiles),&
       iout(ntiles), jout(ntiles)
  real,    intent(out) :: dataout(km,jmout,imout)
  real, allocatable    :: fractmp1(:), fractmp3(:,:,:), tmp(:)
  integer :: nt

  !f2py intent(in) :: datain, iin, jin, iout, jout, frac, imout, jmout, undef
  !f2py intent(out) :: dataout
  !f2py intent(hide) :: imin, jmin, km, ntiles

  allocate(fractmp1(km), fractmp3(km,jmout,imout),tmp(km))
  dataout=0.0; fractmp1=0.0; fractmp3=0.0

  do nt=1,ntiles
     fractmp1=frac(nt)
     tmp=datain(:,jin(nt),iin(nt))
     where(tmp == undef)
        fractmp1=0.0
        tmp=0.0
     endwhere
     dataout(:,jout(nt),iout(nt))=dataout(:,jout(nt),iout(nt))+tmp*fractmp1
     fractmp3(:,jout(nt),iout(nt))=fractmp3(:,jout(nt),iout(nt))+fractmp1
  end do
  
  where(fractmp3 /= 0.0)
     dataout=dataout/fractmp3
  elsewhere
     dataout=undef
  end where
  
  deallocate(fractmp1, fractmp3, tmp)
end subroutine g2g
