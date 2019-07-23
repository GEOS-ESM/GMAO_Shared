      subroutine varimax2(X,Y,M,K,N)

c     THIS SUBROUTINE ROTATES THE VECTORS X USING THE VARIMAX CRITERION. 
c     ON INPUT,  X CONTAINS THE N EIGENVECTORS, EACH OF LENGTH M, TO BE ROTATED. 
c     ON OUTPUT, X CONTAINS THE ROTATED SYSTEM.
      
c     implicit none
      
      integer M, K, L, N
      real    X(M,N), U(M), H(M)
      real    Y(K,N), V(M), T(K)

      real     VVOLD, COSPHI, SUM, A, B
      real     VVNEW, SINPHI, PHI, C, D
      integer  P, Q, I, MAXITER

c     f2py intent(in) :: M, K, N
c     f2py intent(in,out) :: X, Y
c     f2py intent(hide)  :: M, K, N

c     NORMALIZE COMPONENTS

      do P=1,M
         H(P) =0.0
         do Q=1,N
            H(P)=H(P)+X(P,Q)**2
         enddo
         H(P)=sqrt(H(P))
      enddo

      do P=1,M
      do Q=1,N
         IF(H(P) .NE. 0.0) X(P,Q) = X(P,Q)/H(P)
      enddo
      enddo

c     VARIMAX ITERATION FOR ALL PAIRS OF VECTORS

      maxiter=40

      vvold=0.0
      do Q=1,N
       sum=0.0
       do P=1,M
        sum=sum+X(P,Q)**2
       enddo
       vvold=vvold-sum*sum
       do P=1,M
        vvold=vvold+M*X(P,Q)**4
       enddo
      enddo

      do L=1,MAXITER

      do P=1,N-1
       do Q=P+1,N
        A=0.0
        B=0.0
        C=0.0
        D=0.0
        do I=1,M
         U(I)=X(I,P)**2-X(I,Q)**2
         V(I)=2.0*X(I,P)*X(I,Q)
         A=A+U(I)
         B=B+V(I)
         C=C+U(I)*U(I)-V(I)*V(I)
         D=D+U(I)*V(I)
        enddo
        C=C-(A*A-B*B)/FLOAT(M)
        D=2.0*(D-A*B /FLOAT(M))

        PHI=0.25*atan2(D,C)
        COSPHI=cos(PHI)
        SINPHI=sin(PHI)
        do I=1,M
         U(I)  = X(I,P)*COSPHI+X(I,Q)*SINPHI
         X(I,Q)=-X(I,P)*SINPHI+X(I,Q)*COSPHI
         X(I,P)= U(I)
        enddo
        Do I=1,K
         T(I)  = Y(I,P)*COSPHI+Y(I,Q)*SINPHI
         Y(I,Q)=-Y(I,P)*SINPHI+Y(I,Q)*COSPHI
         Y(I,P)= T(I)
        enddo
       enddo
      enddo

      vvnew=0.0
      do Q=1,N
       sum=0.0
       do P=1,M
        sum=sum+X(P,Q)**2
       enddo
       vvnew=vvnew-sum*sum
       do P=1,M
        vvnew=vvnew+M*X(P,Q)**4
       enddo
      enddo
c
      print *, vvold, vvnew, (VVNEW-VVOLD)/VVNEW
      if((VVNEW-VVOLD)/VVNEW < .0001) exit

      if(L==MAXITER) then
       print *, 'eofs: VARIMAX rotation iteration did not converge'
       stop
      endif

      VVOLD=VVNEW

      enddo

      do Q=1,N
      do P=1,M
       X(P,Q)=X(P,Q)*H(P)
      enddo
      enddo

      return
      end subroutine VARIMAX2
