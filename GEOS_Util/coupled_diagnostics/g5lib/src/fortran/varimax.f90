!!$subroutine varimax(X,Y,M,K,N)
!!$  
!!$  !  THIS SUBROUTINE ROTATES THE VECTORS X USING THE VARIMAX CRITERION. 
!!$  !  ON INPUT,  X CONTAINS THE N EIGENVECTORS, EACH OF LENGTH M, TO BE ROTATED. 
!!$  !  ON OUTPUT, X CONTAINS THE ROTATED SYSTEM.
!!$  
!!$  implicit none
!!$  
!!$  integer, intent(IN) :: M, K, N
!!$  real, intent(INOUT) :: X(M,N),Y(K,N)
!!$  
!!$  real, dimension(M) :: H,U,V
!!$  real, dimension(K) :: T
!!$  real    :: A, B, C, D, VVOLD, VVNEW
!!$  real    :: PHI, COSPHI, SINPHI
!!$  integer :: P, Q, L, MAXITER
!!$  
!!$  !f2py intent(in) :: M, K, N
!!$  !f2py intent(in,out) :: X, Y
!!$  !f2py intent(hide)  :: M, K, N
!!$
!!$  ! NORMALIZE COMPONENTS
!!$  
!!$  do P=1,M
!!$     H(P) =sqrt(sum(X(P,:)**2))
!!$  enddo
!!$  
!!$  do Q=1,N
!!$     where(H/=0.) X(:,Q) = X(:,Q)/H(:)
!!$  enddo
!!$  
!!$  ! VARIMAX ITERATION FOR ALL PAIRS OF VECTORS
!!$  
!!$  maxiter = 40
!!$  vvold =  M*sum(X(:,:)**4) - sum(sum(X**2,1)**2)
!!$  
!!$  do L=1,MAXITER
!!$     do P=1,N-1
!!$        do Q=P+1,N
!!$           U = X(:,P)**2 - X(:,Q)**2
!!$           V = 2.*X(:,P)*X(:,Q)
!!$           A = sum(U)
!!$           B = sum(V)
!!$           C = sum(U*U-V*V) - (A*A-B*B)/FLOAT(M)
!!$           D = 2.*(sum(U*V) - A*B/FLOAT(M))
!!$           PHI    = .25*atan2(D,C)
!!$           COSPHI = cos(PHI)
!!$           SINPHI = sin(PHI)
!!$           U      =  X(:,P)*COSPHI+X(:,Q)*SINPHI
!!$           X(:,Q) = -X(:,P)*SINPHI+X(:,Q)*COSPHI
!!$           X(:,P) = U
!!$           T      =  Y(:,P)*COSPHI+Y(:,Q)*SINPHI
!!$           Y(:,Q) = -Y(:,P)*SINPHI+Y(:,Q)*COSPHI
!!$           Y(:,P) = T
!!$        end do
!!$     enddo
!!$     
!!$     VVNEW =  M*sum(X(:,:)**4) - sum(sum(X**2,1)**2)
!!$     print *, vvold, vvnew, (VVNEW-VVOLD)/VVNEW 
!!$     if((VVNEW-VVOLD)/VVNEW < .0001) exit
!!$     
!!$     if(L==MAXITER) then
!!$        print *, 'eofs: VARIMAX rotation iteration did not converge'
!!$        call exit(10)
!!$     endif
!!$     
!!$     VVOLD = VVNEW
!!$     
!!$  enddo
!!$  
!!$  do Q=1,N
!!$     X(:,Q) = X(:,Q)*H(:)
!!$  end do
!!$  
!!$  return
!!$end subroutine VARIMAX


subroutine varimax(X,Y,M,K,N,NM)

  !  THIS SUBROUTINE ROTATES THE VECTORS X USING THE VARIMAX CRITERION. 
  !  ON INPUT,  X CONTAINS THE N EIGENVECTORS, EACH OF LENGTH M, TO BE ROTATED. 
  !  ON OUTPUT, X CONTAINS THE ROTATED SYSTEM.
  
  implicit none
  
  integer, intent(IN) :: M, K, N, NM
  real*8, intent(INOUT) :: X(M,N),Y(K,N)
  
  real*8, dimension(M) :: H,U,V
  real*8, dimension(K) :: T
  real*8    :: A, B, C, D, VVOLD, VVNEW
  real*8    :: PHI, COSPHI, SINPHI
  integer :: P, Q, L, MAXITER
  
  !f2py intent(in) :: M, K, N, NM
  !f2py intent(in,out) :: X, Y
  !f2py intent(hide)  :: M, K, N

  ! NORMALIZE COMPONENTS

  do P=1,M
     H(P) =sqrt(sum(X(P,:)**2))
  enddo
  
  do Q=1,N
     where(H/=0.) X(:,Q) = X(:,Q)/H(:)
  enddo
  
  ! VARIMAX ITERATION FOR ALL PAIRS OF VECTORS
  
  maxiter = 40
  vvold =  M*sum(X(:,:)**4) - sum(sum(X**2,1)**2)
  
  do L=1,MAXITER
     
     do P=1,NM-1 !N,N-NM+2,-1
        do Q=P+1,NM !P-1,N-NM+1,-1
           U = X(:,P)**2 - X(:,Q)**2
           V = 2.*X(:,P)*X(:,Q)
           A = sum(U)
           B = sum(V)
           C = sum(U*U-V*V) - (A*A-B*B)/FLOAT(M)
           D = 2.*(sum(U*V) - A*B/FLOAT(M))
           PHI    = .25*atan2(D,C)
           COSPHI = cos(PHI)
           SINPHI = sin(PHI)
           U      =  X(:,P)*COSPHI+X(:,Q)*SINPHI
           X(:,Q) = -X(:,P)*SINPHI+X(:,Q)*COSPHI
           X(:,P) = U
           T      =  Y(:,P)*COSPHI+Y(:,Q)*SINPHI
           Y(:,Q) = -Y(:,P)*SINPHI+Y(:,Q)*COSPHI
           Y(:,P) = T
        end do
     enddo
     
     VVNEW =  M*sum(X(:,:)**4) - sum(sum(X**2,1)**2)
     print *, vvold, vvnew, (VVNEW-VVOLD)/VVNEW 
     if((VVNEW-VVOLD)/VVNEW < .0001) exit
     
     if(L==MAXITER) then
        print *, 'eofs: VARIMAX rotation iteration did not converge'
        call exit(10)
     endif
     
     VVOLD = VVNEW
     
  enddo
  
  do Q=1,N
     X(:,Q) = X(:,Q)*H(:)
  end do
  
  return
end subroutine VARIMAX
