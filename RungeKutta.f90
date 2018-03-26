



!******************************************************************
!			PART ONE: SOLVERS
!				
!******************************************************************
module DIFF_SOLVER

interface fwdEuler
module procedure fwdEuler, fwdEulerVec
end interface

interface rk4
	module procedure rk4, rk4Vec
end interface

contains

subroutine fwdEuler (t0, u0, dt, f, u)
implicit none
real ( kind = selected_real_kind(8) ) :: dt,t0,f0,u,u0

interface
subroutine f(t,u,uprime)
real ( kind = selected_real_kind(8) ) :: t,u,uprime
end subroutine f
end interface

call f(t0,u0,f0)
u=u0+dt*f0

end subroutine

subroutine fwdEulerVec (m,t0, u0, dt, fvec, u,beta, gamm)
implicit none
integer :: m
real ( kind = selected_real_kind(8) ) :: dt,t0,f0(m),u(m),u0(m), beta, gamm

interface
subroutine fvec(m,t,u,uprime,beta, gamm)
integer :: m
real (kind=selected_real_kind(8) ) :: u(m), uprime(m)
real (kind=selected_real_kind(8) ) :: t
real (kind=selected_real_kind(8) ) :: beta,gamm
end subroutine fvec
end interface

!call fvec(m,t0,u0,f0,beta,gamm)
  f0(1) = -1*beta*u0(1)*u0(2)
  f0(2) = beta*u0(1)*u0(2)-gamm*u0(2)
  f0(3) = gamm*u0(2)
u(1:m)=u0(1:m)+dt*f0(1:m)

end subroutine


subroutine rk4 ( t0, u0, dt, f, u )

  implicit none

  real ( kind = selected_real_kind(8) ) :: dt, f0,f1,f2,f3,t0,t1,t2,t3,u,u0,u1,u2,u3
!Interface the function describing the ODE
  interface
  subroutine f(t,u,uprime)
  real ( kind=selected_real_kind(8) )  :: t,u,uprime
  end subroutine f
  end interface

!Perform RK order 4
  call f ( t0, u0, f0 )
  
  t1 = t0 + 0.5*dt
  u1 = u0 + 0.5*dt * f0
  call f ( t1, u1, f1 )
  
  t2 = t0 + 0.5*dt
  u2 = u0 + 0.5*dt * f1
  call f ( t2, u2, f2 )
  
  t3 = t0 + dt
  u3 = u0 + dt * f2
  call f ( t3, u3, f3 ) 
  u = u0 + dt * ( f0 + 2.0 * f1 + 2.0 * f2 + f3 ) / 6.0
end subroutine


subroutine rk4Vec (m, t0, u0, dt, fvec, u,beta,gamm)

  implicit none

  integer :: m

  real ( kind=selected_real_kind(8) ) :: dt,f0(m),f1(m),f2(m),f3(m)
  real ( kind=selected_real_kind(8) ) :: t0,t1,t2,t3
  real ( kind=selected_real_kind(8) ) :: u(m),u0(m),u1(m),u2(m),u3(m),uprime(3)
  real ( kind=selected_real_kind(8) ) :: beta,gamm

  interface

  subroutine fvec(m,t,u,uprime, beta, gamm)
  integer:: m
  real ( kind=selected_real_kind(8) )  :: t,u(m),uprime(m),beta,gamm
  end subroutine fvec

  end interface  

  !f0(1) = -1*beta*u0(1)*u0(2)
  !f0(2) = beta*u0(1)*u0(2)-gamm*u0(2)
  !f0(3) = gamm*u0(2)
  call fvec (m, t0, u0, f0,beta,gamm )

  t1 = t0 + 0.5*dt
  u1(1:m) = u0(1:m) + 0.5*dt * f0(1:m)
  call fvec (m, t1, u1, f1 ,beta,gamm)
  !f1(1) = -1*beta*u1(1)*u1(2)
  !f1(2) = beta*u1(1)*u1(2)-gamm*u1(2)
  !f1(3) = gamm*u1(2)
  u2(1:m) = u0(1:m) + 0.5*dt * f1(1:m)
  call fvec (m, t1, u2, f2,beta,gamm )
  !f2(1) = -1*beta*u2(1)*u2(2)
  !f2(2) = beta*u2(1)*u2(2)-gamm*u2(2)
  !f2(3) = gamm*u2(2)

  t3 = t0 + dt
  u3(1:m) = u0(1:m) + dt * f2(1:m)
  call fvec (m, t3, u3, f3,beta,gamm )
  !f3(1) = -1*beta*u3(1)*u3(2)
  !f3(2) = beta*u3(1)*u3(2)-gamm*u3(2)
  !f3(3) = gamm*u3(2)
  u(1:m) = u0(1:m) + ( dt /6.0 ) * ( &
                 f0(1:m) &
     + 2.0* f1(1:m) &
     + 2.0* f2(1:m) &
     +           f3(1:m) )

  return
end subroutine
end module

!******************************************************************
!		    PART TWO: DIFFERENTIAL EQUATIONS
!				
!******************************************************************

module DIFF

implicit none

contains
!EQUATION TO BE SOLVED GOES HERE
subroutine f (t,u,uprime)
real (kind=selected_real_kind(8) ) u
real (kind=selected_real_kind(8) ) t
real (kind=selected_real_kind(8) ) uprime
uprime = 3.0*t**2
end subroutine f
!VECTORIAL VERSION OF DIFFERENTIAL EQUATION
subroutine fvec (m,t,u,uprime, beta, gamm)
integer :: m
real (kind=selected_real_kind(8) ) :: u(m), uprime(m)
real (kind=selected_real_kind(8) ) :: t
real (kind=selected_real_kind(8) ) :: beta,gamm
uprime(1) = -1*beta*u(1)*u(2)
uprime(2) = beta*u(1)*u(2)-gamm*u(2)
uprime(3) = gamm*u(2)
end subroutine fvec
end module
!******************************************************************
!		  PART THREE: RANDOM SAMPLERS
!				
!******************************************************************

module RANDOM

implicit none

contains
subroutine drawDirichlet(m,param,theta)
implicit none
!Subroutine that draws from m-dimensional Dirichlet 
!using gamma distribution Gamm(alpha,1)
integer, intent(in) :: m
integer	:: i
real (kind=selected_real_kind(8) ), intent(in) :: param(m)
real (kind=selected_real_kind(8) ), intent(out) ::theta(m)
real (kind=selected_real_kind(8) ) :: y,s
call init_random_seed()
s=0.0
do i=1,m
y= rand_gamma(real(1.0,kind=selected_real_kind(8)),param(i))
theta(i)=y
s=s+y
end do
theta=(1.0/s)*theta
end subroutine drawDirichlet

function log_gamma(shape,scale) result(ans)
real (kind=selected_real_kind(8) ) :: ans,r,lambda,w,h,eta,z,u2,u
real (kind=selected_real_kind(8) ) :: shape,scale
	  call init_random_seed
	  call random_number(u2)
	  h=-1.
	  eta=1.
	  do while((h/eta).le.u2)
		w=shape/(exp(1.0)*(1-shape))
		lambda=(1.0/shape)-1.0
		r=1./(1.+w)
		call init_random_seed()
		call random_number(u)
		if(u.le.r)then
		z=-LOG(u/r)
		else
		call init_random_seed()
		call random_number(u)
		z=log(u)/lambda
		end if
		h=(1./GAMMA(shape+1))*EXP(-z-EXP(-z/shape))
		if(z .ge. 0.) then
		eta=(1./GAMMA(shape+1))*EXP(-z)
		else
		eta=(1./GAMMA(shape+1))*w*lambda*EXP(lambda*z)
		end if
		call init_random_seed
	    call random_number(u2)
		end do
		ans=-z/shape
end function

RECURSIVE FUNCTION rand_gamma(shape, scale) RESULT(ans)
	implicit none
	real (kind=selected_real_kind(8) ) :: ans,w
	real (kind=selected_real_kind(8) ), intent(in) :: shape, scale !shape and scale
     	real (kind=selected_real_kind(8) ) :: u,d,c,x,xsq,g,v
     IF (shape <= 0.0) THEN

        WRITE(*,*) "Shape PARAMETER must be positive"
      END IF
      IF (scale <= 0.0) THEN
        WRITE(*,*) "Scale PARAMETER must be positive"
      END IF
      
      IF (shape >= 1.0) THEN
        d = shape - 1.0/3.0
        c = 1.0/(9.0*d)**0.5
        DO while (.true.)
            x = rand_normal(0.0_8, 1.0_8)
		!print*,x
            v = 1.0 + c*x
            DO while (v <= 0.0)
                x = rand_normal(0.0_8, 1.0_8)
                v = 1.0 + c*x
            END DO

            v = v*v*v
	    call init_random_seed()
            CALL RANDOM_NUMBER(u)
            xsq = x*x
            IF ((u < 1.0 -.0331*xsq*xsq) .OR.  &
              (log(u) < 0.5*xsq + d*(1.0 - v + log(v))) )then
                ans=scale*d*v
                RETURN
            END IF

        END DO
      ELSE
        g = rand_gamma(shape+1.0, 1.0_8)
	call init_random_seed()
        CALL RANDOM_NUMBER(w)
        ans=scale*g*(w)**(1.0/shape)
        RETURN
      END IF
END FUNCTION


SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock,clock2
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
            real	:: y
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
	    clock2=clock
            seed = clock + 9929 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
	    
	    do while(clock .eq. clock2)
		CALL SYSTEM_CLOCK(COUNT=clock2)
	    end do
END SUBROUTINE

FUNCTION rand_beta(a, b) RESULT(ans)
	implicit none
	real*8, intent(in) :: a,b
	real*8 :: ans
      	real*8 :: u,v
		call init_random_seed()
      IF ((a <= 0.0) .OR. (b <= 0.0)) THEN

        WRITE(*,*) "Beta PARAMETERs must be positive"
      END IF
      
       u = rand_gamma(a, 1.0_8)
       v = rand_gamma(b, 1.0_8)
       ans = u / (u + v)
END FUNCTION

function rand_normal ( a, b)
  implicit none

  real ( kind = selected_real_kind(8) ) :: a,b,r1,r2,rand_normal,two
  real ( kind = selected_real_kind(8) ), parameter :: pi = 4*ATAN(real(1.0,kind=selected_real_kind(8)))
  
  call init_random_seed()
  call random_number(r1)
  call random_number(r2)
  two=2.0
  rand_normal = sqrt ( - two * log ( r1 ) ) * cos ( two * pi * r2 )

  rand_normal = a + b * rand_normal

  return
end function

end module


!******************************************************************
!		    		MAIN
!				
!******************************************************************
!Feed this subroutine the estimated parameters (beta, gamma, kappa, lambda) and the 
!number of iterates K
module DBSSM_MOD


contains

subroutine DBSSM(K,dt, t0, beta,gamm,kappa,lambda,theta0,theta,y)

use RANDOM
use DIFF_SOLVER
use DIFF

real ( kind=selected_real_kind(8) ), intent (in)   :: dt,t0
real ( kind=selected_real_kind(8) ), intent (inout):: theta(3,K), y(K),theta0(3)
real ( kind=selected_real_kind(8) ), intent (in)   :: beta,gamm,kappa, lambda
integer 					   :: seed,clock,count_rate,count_max,i,K

call rk4vec(3, t0, theta0, dt, fvec, theta(1,1),beta,gamm)
 y(1)=rand_beta(lambda*theta(2,1), lambda*(1-theta(2,1)))
do i=2,K
 call rk4Vec(3, t0, theta(:,i-1), dt, fvec, theta(:,i),beta,gamm)
 call drawDirichlet(3, kappa*theta(:,i),theta(:,i))
 y(i)=rand_beta(lambda*theta(2,i), lambda*(1-theta(2,i)))
end do

end subroutine

function DBSSM_FUNC(K,dt, t0, beta,gamm,kappa,lambda,theta0) result(y_and_theta)
 use RANDOM
 use DIFF_SOLVER
 use DIFF

real ( kind=selected_real_kind(8) )   :: dt,t0
real ( kind=selected_real_kind(8) )   :: theta0(3)
real ( kind=selected_real_kind(8) )		   :: y_and_theta(4,K)
real ( kind=selected_real_kind(8) )  :: beta,gamm,kappa, lambda
integer 					   :: seed,clock,count_rate,count_max,i,K
 y_and_theta=0.0
 call rk4Vec(3, t0, theta0, dt, fvec, y_and_theta(1:3,1),beta,gamm)
 y_and_theta(4,1)=rand_beta(lambda*y_and_theta(2,1), lambda*(1.0-y_and_theta(2,1)))
 !print *, y_and_theta
do i=2,K
 call rk4Vec(3, t0, y_and_theta(1:3,i-1), dt, fvec, y_and_theta(1:3,i),beta,gamm)
 call drawDirichlet(3, kappa*y_and_theta(1:3,i),y_and_theta(1:3,i))
 y_and_theta(4,i)=rand_beta(lambda*y_and_theta(2,i), lambda*(1-y_and_theta(2,i)))
end do
i=1
end function

end module

!program main

!use DBSSM_MOD
!use random
!real ( kind=selected_real_kind(8) )   :: dt,t0,s
!real ( kind=selected_real_kind(8) )   :: theta0(3)
!real ( kind=selected_real_kind(8) )		   :: y_and_theta(4,10),x,theta(3),param(3)
!real ( kind=selected_real_kind(8) )   :: beta,gamm,kappa, lambda
!integer 	:: seed,clock,count_rate,count_max,i,K,j
!integer		:: i,j
!t0=0.0
!s=0.0
!dt=1.0
!beta=2.0
!gamm=1.4
!lambda=2000
!kappa=1.0
!theta0(1)=0.9
!theta0(2)=0.0002
!theta0(3)=0.0998
!param=(/1.0,1.0,1.0/)
!j=0
!x=0.0
!do i=1,100
!y_and_theta=DBSSM_FUNC(10,dt, t0, beta,gamm,kappa,lambda,theta0)
!print *,i, "Y IS: ",y_and_theta(4,:)
!if( sqrt(y_and_theta(2,10)*(1.0-y_and_theta(2,10))/(1.0+lambda))- abs(y_and_theta(4,10)-y_and_theta(2,10))<0) then
!j=j+1
!end if
!end do
!open(unit=10,file=Dirichlet1.txt)
!open(unit=11,file=Dirichlet2.txt)
!open(unit=12,file=Dirichlet3.txt)
!do i=1,5000
!call drawDirichlet(3,param, theta)
!write (10,*) theta(1)
!write (11,*) theta(2)
!write (12,*) theta(3)
!end do
!s=0.0
!do i=1,5000
!s=s+rand_normal(real(2.0,8),real(1.0,8))
!end do
!print *, s/5000
!x=x-s
!print *, "VARIANCE IS: ", sum(x**2)/49999.0

!print *, "OUTSIDE OF ST. DEV. ", j, " OUT OF 100 TIMES"
!do i=1,100
!print *, rand_beta(real(0.0003,8)*lambda,lambda*0.9997)
!end do
!print*, y_and_theta
!print *, "test performed"
!end program
