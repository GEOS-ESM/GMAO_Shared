      subroutine sspg(n,x,m,eps,eps2,maxit,maxfc,output,f,cginfn,cgtwon,
     +                iter,fcnt,gcnt,flag)

C     Subroutine SPG implements the Spectral Projected Gradient
C     Method (Version 2: "feasible continuous projected path") to
C     find the local minimizers of a given function with convex
C     constraints, described in the paper "Nonmonotone spectral
C     projected gradient methods for convex sets", SIAM J. on
C     Optimization 10(4), 1196-1211, 2000, by E.G. Birgin, J.M.
C     Martinez and M. Raydan.
C
C     The user must supply the external subroutines spg_evalf, spg_evalg 
C     and spg_proj to evaluate the objective function and its gradient 
C     and to project an arbitrary point onto the feasible region.
C
C     This version 17 JAN 2000 by E.G.Birgin, J.M.Martinez and M.Raydan.
C     Reformatted 03 OCT 2000 by Tim Hopkins.
C     Final revision 04 MAY 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.
C     Nothing is every final! Revised on 01Sep2006, R. Todling:
C        - placed print out of summary in the end
C        - print out two-norm
C        - turned main arrays into allocatable ones
C
C     On Entry:
C
C     n     integer,
C           size of the problem
C
C     x     real x(n),
C           initial guess
C
C     m     integer,
C           number of previous function values to be considered 
C           in the nonmonotone line search.
C
C     eps   real,
C           stopping criterion: ||continuous grad||_inf < eps
C
C     eps2  real,
C           stopping criterion: ||continuous grad||_2 < eps2
C
C     maxit integer,
C           maximum number of iterations
C
C     maxfc integer,
C           maximum number of function evaluations
C
C     output logical,
C           true: print some information at each iteration,
C           false: no print.
C
C     On Return:
C
C     x     real X(N),
C           approximation to the local minimizer
C
C     f     real,
C           function value at the approximation to the local
C           minimizer
C
C     cginfn real,
C           ||continuous grad||_inf at the final iteration
C
C     cgtwon real,
C           ||continuous grad||_2^2 at the final iteration
C
C     iter  integer,
C           number of iterations
C
C     fcnt  integer,
C           number of function evaluations
C
C     gcnt  integer,
C           number of gradient evaluations
C
C     flag  integer,
C           termination parameter:
C           0= convergence with continuous gradient infinite-norm,
C           1= convergence with continuous gradient 2-norm,
C           2= too many iterations,
C           3= too many function evaluations,
C           4= error in spg_proj subroutine,
C           5= error in spg_evalf subroutine,
C           6= error in spg_evalg subroutine.
C           7= n larger the internal maximum.
C          98= cannot   allocate working space.
C          99= cannot deallocate working space.

C     PARAMETERS
      real lmin
      parameter (lmin=1.0e-30)
      real lmax
      parameter (lmax=1.0e+30)
      integer nmax
!     parameter (nmax=100000)
      integer mmax
      parameter (mmax=100)

C     SCALAR ARGUMENTS
      real cginfn,cgtwon,eps,eps2,f
      real summary(maxit,4)
      integer fcnt,flag,gcnt,iter,m,maxfc,maxit,n
      logical output

C     ARRAY ARGUMENTS
      real x(n)

C     LOCAL SCALARS
      real fbest,fnew,gtd,lambda,sts,sty
      integer i,lsflag,inform

C     LOCAL 
      
      real, allocatable :: cg(:),g(:),gnew(:),s(:),y(:),
     +                                 d(:),xbest(:),xnew(:)
      real lastfv(0:mmax-1)

C     EXTERNAL SUBROUTINES
      external sspg_ls,spg_evalf,spg_evalg,spg_proj

C     INTRINSIC FUNCTIONS
      intrinsic abs,max,min,mod

C     INITIALIZATION

      nmax = n
      allocate ( cg(nmax),g(nmax),gnew(nmax),s(nmax),y(nmax),
     +           d(nmax),xbest(nmax),xnew(nmax), stat=inform )
      if (inform .ne. 0) then
          flag = 98
          return
      endif

      iter = 0
      fcnt = 0
      gcnt = 0

      do i = 0,m - 1
          lastfv(i) = real(-1.0d+99)
      end do

C     PROJECT INITIAL GUESS

      call sspg_proj(n,x,inform)

      if (inform .ne. 0) then
        ! ERROR IN spg_PROJ SUBROUTINE, STOP
          flag = 4
          go to 200
      end if

C     INITIALIZE BEST SOLUTION

      do i = 1,n
          xbest(i) = x(i)
      end do

C     EVALUATE FUNCTION AND GRADIENT

      call sspg_evalf(n,x,f,inform)
      fcnt = fcnt + 1

      if (inform .ne. 0) then
        ! ERROR IN spg_EVALF SUBROUTINE, STOP
          flag = 5
          go to 200
      end if

      call sspg_evalg(n,x,g,inform)
      gcnt = gcnt + 1

      if (inform .ne. 0) then
        ! ERROR IN spg_EVALG SUBROUTINE, STOP
          flag = 6
          go to 200
      end if

C     STORE FUNCTION VALUE FOR THE NONMONOTONE LINE SEARCH

      lastfv(0) = f

C     INITIALIZE BEST FUNCTION VALUE

      fbest = f

C     COMPUTE PROJECTED CONTINUOUS GRADIENT (AND ITS NORMS)

      do i = 1,n
          cg(i) = x(i) - g(i)
      end do

      call sspg_proj(n,cg,inform)

      if (inform .ne. 0) then
        ! ERROR IN spg_PROJ SUBROUTINE, STOP
          flag = 4
          go to 200
      end if

      cgtwon = 0.0D0
      cginfn = 0.0D0
      do i = 1,n
          cg(i) = cg(i) - x(i)
          cgtwon = cgtwon + cg(i)**2
          cginfn = max(cginfn,abs(cg(i)))
      end do

C     PRINT ITERATION INFORMATION

      if (output) then
          summary(iter,1) = real(iter)
          summary(iter,2) = f
          summary(iter,3) = cginfn
          summary(iter,4) = cgtwon
      end if

C     DEFINE INTIAL SPECTRAL STEPLENGTH

      if (cginfn .ne. 0.0d0) then
          lambda =  min(lmax,max(lmin,1.0d0/cginfn))
      end if

C     MAIN LOOP

C     TEST STOPPING CRITERIA

 100  continue

      if (cginfn .le. eps) then
        ! GRADIENT INFINITE-NORM STOPPING CRITERION SATISFIED, STOP
          flag = 0
          go to 200
      end if

      if (cgtwon .le. eps2**2) then
        ! GRADIENT 2-NORM STOPPING CRITERION SATISFIED, STOP
          flag = 1
          go to 200
      end if

      if (iter .gt. maxit) then
        ! MAXIMUM NUMBER OF ITERATIONS EXCEEDED, STOP
          flag = 2
          go to 200
      end if

      if (fcnt .gt. maxfc) then
        ! MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED, STOP
          flag = 3
          go to 200
      end if

C     DO AN ITERATION

      iter = iter + 1

C     COMPUTE THE SPECTRAL PROJECTED GRADIENT DIRECTION
C     AND <G,D>

      do i = 1,n
          d(i) = x(i) - lambda*g(i)
      end do

      call sspg_proj(n,d,inform)

      if (inform .ne. 0) then
        ! ERROR IN spg_PROJ SUBROUTINE, STOP
          flag = 4
          go to 200
      end if

      gtd = 0.0d0
      do i = 1,n
          d(i) = d(i) - x(i)
          gtd = gtd + g(i)*d(i)
      end do

C     NONMONOTONE LINE SEARCH

      call sspg_ls(n,x,f,d,gtd,m,lastfv,maxfc,fcnt,fnew,xnew,lsflag)

      if (lsflag .eq. 3) then
        ! THE NUMBER OF FUNCTION EVALUATIONS WAS EXCEEDED 
        ! INSIDE  THE LINE SEARCH, STOP
          flag = 3
          go to 200
      end if

C     SET NEW FUNCTION VALUE AND SAVE IT FOR THE NONMONOTONE
C     LINE SEARCH

      f = fnew
      lastfv(mod(iter,m)) = f

C     COMPARE THE NEW FUNCTION VALUE AGAINST THE BEST FUNCTION 
C     VALUE  AND, IF SMALLER, ACTUALIZE THE BEST FUNCTION 
C     VALUE AND THE CORRESPONDING BEST POINT

      if (f .lt. fbest) then
          fbest = f
          do i = 1,n
              xbest(i) = xnew(i)
          end do
      end if

C     EVALUATE THE GRADIENT AT THE NEW ITERATE

      call sspg_evalg(n,xnew,gnew,inform)
      gcnt = gcnt + 1

      if (inform .ne. 0) then
        ! ERROR IN spg_EVALG SUBROUTINE, STOP
          flag = 6
          go to 200
      end if

C     COMPUTE S = XNEW - X, Y = GNEW - G, <S,S>, <S,Y>,
C     THE CONTINUOUS GRADIENT AND ITS NORMS

      sts = 0.0d0
      sty = 0.0d0
      do i = 1,n
          s(i) = xnew(i) - x(i)
          y(i) = gnew(i) - g(i)
          sts = sts + s(i)*s(i)
          sty = sty + s(i)*y(i)
          x(i) = xnew(i)
          g(i) = gnew(i)
          cg(i) = x(i) - g(i)
      end do

      call sspg_proj(n,cg,inform)

      if (inform .ne. 0) then
        ! ERROR IN spg_PROJ SUBROUTINE, STOP
          flag = 4
          go to 200
      end if

      cgtwon = 0.0D0
      cginfn = 0.0D0
      do i = 1,n
          cg(i) = cg(i) - x(i)
          cgtwon = cgtwon + cg(i)**2
          cginfn = max(cginfn,abs(cg(i)))
      end do

C     PRINT ITERATION INFORMATION

      if (output) then
          summary(iter,1) = real(iter)
          summary(iter,2) = f
          summary(iter,3) = cginfn
          summary(iter,4) = cgtwon
      end if

C     COMPUTE SPECTRAL STEPLENGTH

      if (sty .le. 0.0d0) then
          lambda = lmax

      else
          lambda = min(lmax,max(lmin,sts/sty))
      end if

C     FINISH THE ITERATION

      go to 100

C     STOP

 200  continue

C     SET X AND F WITH THE BEST SOLUTION AND ITS FUNCTION VALUE

      f = fbest
      do i = 1,n
          x(i) = xbest(i)
      end do

C     SUMMARIZE RESULTS
                                                                                                                           
      if (output) then
           do i = 1, iter
              write(*,fmt=9010) nint(summary(i,1)),(summary(i,j),j=2,4)
           end do
      endif

C     CLEAN UP

      deallocate ( cg,g,gnew,s,y,d,xbest,xnew, stat=inform )
      if (inform .ne. 0) then
          flag = 99
          return
      endif

 9010 format ('ITER= ',I10,' F= ',1P,D17.10,' CGINFNORM= ',1P,D17.10,
     .        ' CGTWONORM= ',D17.10)

      end


      subroutine sspg_ls(n,x,f,d,gtd,m,lastfv,maxfc,fcnt,fnew,xnew,flag)

C     Subroutine LS implements a nonmonotone line search with
C     safeguarded quadratic interpolation.
C
C     This version 17 JAN 2000 by E.G.Birgin, J.M.Martinez and M.Raydan.
C     Reformatted 03 OCT 2000 by Tim Hopkins.
C     Final revision 30 APR 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.
C
C     On Entry:
C
C     n     integer,
C           size of the problem
C
C     x     real x(n),
C           initial guess
C
C     f     real,
C           function value at the actual point
C
C     d     real d(n),
C           search direction
C
C     gtd   real,
C           internal product <g,d>, where g is the gradient at x,
C
C     m     integer,
C           number of previous function values to be considered 
C           in the nonmonotone line search,
C
C     lastfv real lastfv(m),
C           last m function values,
C
C     maxfc integer,
C           maximum number of function evaluations,
C
C     fcnt  integer,
C           actual number of fucntion evaluations.
C
C     On Return:
C
C     fcnt  integer,
C           actual number of fucntion evaluations,
C
C     fnew  real,
C           function value at the new point
C
C     xnew  real xnew(n),
C           new point,
C
C     flag  integer,
C           0= convergence with nonmonotone Armijo-like criterion,
C           3= too many function evaluations,
C           5= error in spg_evalf subroutine.

C     PARAMETERS
      real gamma
      parameter (gamma=1.0d-04)

C     SCALAR ARGUMENTS
      real f,fnew,gtd
      integer maxfc,fcnt,m,n,flag

C     ARRAY ARGUMENTS
      real d(n),lastfv(0:m-1),x(n),xnew(n)

C     LOCAL SCALARS
      real alpha,atemp,fmax
      integer i,inform

C     EXTERNAL SUBROUTINES
      external spg_evalf

C     INTRINSIC FUNCTIONS
      intrinsic max

C     INITIALIZATION

C     COMPUTE THE MAXIMUM FUNCTIONAL VALUE OF THE LAST M ITERATIONS

      fmax = lastfv(0)
      do i = 1,m - 1
          fmax = max(fmax,lastfv(i))
      end do

C     COMPUTE FIRST TRIAL

      alpha = 1.0d0

      do i = 1,n
          xnew(i) = x(i) + d(i)
      end do

      call sspg_evalf(n,xnew,fnew,inform)
      fcnt = fcnt + 1

      if (inform .ne. 0) then
        ! ERROR IN spg_EVALF SUBROUTINE, STOP
          flag = 5
          go to 200
      end if

C     MAIN LOOP

 100  continue

C     TEST STOPPING CRITERIA

      if (fnew .le. fmax + gamma*alpha*gtd) then
        ! NONMONOTONE ARMIJO-LIKE STOPPING CRITERION SATISFIED, STOP
          flag = 0
          go to 200
      end if

      if (fcnt .ge. maxfc) then
        ! MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED, STOP
          flag = 3
          go to 200
      end if 

C     DO AN ITERATION

C     SAFEGUARDED QUADRATIC INTERPOLATION

      if (alpha .le. 0.1d0) then
          alpha = alpha/2.0d0

      else
          atemp = (-gtd*alpha**2) / (2.0d0*(fnew-f-alpha*gtd))
          if (atemp .lt. 0.1d0 .or. atemp .gt. 0.9d0*alpha) then
              atemp = alpha/2.0d0
          end if
          alpha = atemp
      end if

C     COMPUTE TRIAL POINT 

      do i = 1,n
          xnew(i) = x(i) + alpha*d(i)
      end do

C     EVALUATE FUNCTION

      call sspg_evalf(n,xnew,fnew,inform)
      fcnt = fcnt + 1

      if (inform .ne. 0) then
        ! ERROR IN spg_EVALF SUBROUTINE, STOP
          flag = 5
          go to 200
      end if

C     ITERATE

      go to 100

C     STOP

 200  continue

      end

