
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
implicit none
integer :: m
real (kind=selected_real_kind(8) ) :: u(m), uprime(m)
real (kind=selected_real_kind(8) ) :: t
real (kind=selected_real_kind(8) ) :: beta,gamm
end subroutine fvec
end interface

call fvec(m,t0,u0,f0,beta,gamm)
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
  real ( kind=selected_real_kind(8) ) :: u(m),u0(m),u1(m),u2(m),u3(m)
  real ( kind=selected_real_kind(8) ) :: beta,gamm

  interface

  subroutine fvec(m,t,u,uprime, beta, gamm)
  integer:: m
  real ( kind=selected_real_kind(8) )  :: t,u(m),uprime(m),beta,gamm
  end subroutine fvec

  end interface  

  call fvec (m, t0, u0, f0,beta,gamm )

  t1 = t0 + 0.5*dt
  u1(1:m) = u0(1:m) + 0.5*dt * f0(1:m)
  call fvec (m, t1, u1, f1 ,beta,gamm)
  u2(1:m) = u0(1:m) + 0.5*dt * f1(1:m)
  call fvec (m, t1, u2, f2,beta,gamm )

  t3 = t0 + dt
  u3(1:m) = u0(1:m) + dt * f2(1:m)
  call fvec (m, t3, u3, f3,beta,gamm )
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
implicit none
real (kind=selected_real_kind(8) ) u
real (kind=selected_real_kind(8) ) t
real (kind=selected_real_kind(8) ) uprime
uprime = 3.0*t**2
end subroutine f
!VECTORIAL VERSION OF DIFFERENTIAL EQUATION
subroutine fvec (m,t,u,uprime, beta, gamm)
implicit none
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


RECURSIVE FUNCTION rand_gamma(shape, scale) RESULT(ans)
	implicit none
	real (kind=selected_real_kind(8) ) :: ans
	real (kind=selected_real_kind(8) ), intent(in) :: shape, scale !shape and scale
     	real (kind=selected_real_kind(8) ) :: u,w,d,c,x,xsq,g,v
call init_random_seed()
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
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
	    
	    do while(clock .eq. clock2)
		CALL SYSTEM_CLOCK(COUNT=clock2)
	    end do
END SUBROUTINE

FUNCTION rand_beta(a, b) RESULT(ans)
	implicit none
	real(kind=selected_real_kind(8) ), intent(in) :: a,b
	real(kind=selected_real_kind(8) )			  :: ans
    real(kind=selected_real_kind(8) ) 			  :: u,v
	call init_random_seed()
      IF ((a <= 0.0) .OR. (b <= 0.0)) THEN

        WRITE(*,*) "Beta PARAMETERs must be positive"
      END IF
      
       u = rand_gamma(a, 1.0_8)
       v = rand_gamma(b, 1.0_8)
       ans = u / (u + v)
	   
END FUNCTION

FUNCTION rand_normal(mean,stdev) RESULT(c)
	implicit none
	real(kind=selected_real_kind(8) ), intent(in) :: mean, stdev
    real(kind=selected_real_kind(8) ) :: theta,r,c
	real(kind=selected_real_kind(8) ), dimension(2) :: temp
	real(kind=selected_real_kind(8) ),parameter :: PI_8 = 4*ATAN(1.0_8)
	

	call init_random_seed()
	
        CALL RANDOM_NUMBER(temp)
	!print*,temp
        r=(-2.0*log(temp(1)))**0.5
        theta = 2.0*PI_8*temp(2)
        c = mean+stdev*r*sin(theta)

      
END FUNCTION

end module


!******************************************************************
!		    		MAIN
!				
!******************************************************************
!Feed this subroutine the estimated parameters (beta, gamma, kappa, lambda) and the 
!number of iterates K
subroutine DBSSM(K,dt, t0, beta,gamm,kappa,lambda,theta0,theta,y)
use RANDOM
use DIFF_SOLVER
use DIFF

real ( kind=selected_real_kind(8) ), intent (in)   :: dt,t0
real ( kind=selected_real_kind(8) ), intent (inout):: theta(3,K), y(K),theta0(3)
real ( kind=selected_real_kind(8) ), intent (in)   :: beta,gamm,kappa, lambda
integer 					   :: seed,clock,count_rate,count_max,i,K

 call init_random_seed()
 call rk4vec(3, t0, theta0, dt, fvec, theta(1,1),beta,gamm)
 y(1)=rand_beta(lambda*theta(2,1), lambda*(1-theta(2,1)))
do i=2,K
 call init_random_seed()
 call rk4vec(3, t0, theta(:,i-1), dt, fvec, theta(:,i),beta,gamm)
 call drawDirichlet(3, kappa*theta(:,i),theta(:,i))
 y(i)=rand_beta(lambda*theta(2,i), lambda*(1-theta(2,i)))
end do

end subroutine





program main
use random
implicit none
real*8 :: S0,dt,I0,R0,kappa,lambda,beta,gamm,t0,rho,minmiss,currentmiss,mean
real*8, dimension(3) :: theta0
real*8, dimension(:,:), allocatable :: y_and_theta, y_and_theta_best
integer :: K,i,j
K = 30
t0 = 0 !voorlopig
allocate(y_and_theta(4,K))
allocate(y_and_theta_best(4,K))
open(unit=1,file='I0.txt')
open(unit=2,file='R0.txt')
open(unit=3,file='kappa.txt')
open(unit=4,file='lambda.txt')
open(unit=5,file='betas.txt')
open(unit=7,file='rho.txt') !gamm
minmiss = huge(0.0)

do j=1,5000
print*,j

read(1,*) I0
read(2,*) R0
S0 = 0.9
theta0(1) = S0

theta0(2) = I0
theta0(3) = R0
read(3,*) kappa
read(4,*) lambda
read(5,*) beta
read(7,*) rho
gamm = rho*beta


do i = 1,20
	y_and_theta = DBSSM_FUNC(K,dt, t0, beta,gamm,kappa,lambda,theta0)
	currentmiss = calc_likelihood2(y_and_theta)

	if (currentmiss .LT. minmiss) then
		minmiss = currentmiss
		y_and_theta_best(:,:) = y_and_theta(:,:)
	endif

enddo
print *,minmiss	

enddo
print*,'part1'
print*,y_and_theta_best(1,:)
print*,'part2'
print*,y_and_theta_best(2,:)
print*,'part3'
print*,y_and_theta_best(3,:)
print*,'part4'
print*,y_and_theta_best(4,:)

print *, minmiss

	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	close(7)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

function calc_likelihood(y_and_theta) result (miss)
	implicit none
	real*8, dimension(4,K) :: y_and_theta
	real*8 :: miss
	integer :: i
	real*8, dimension(17) :: y_week
	y_week = (/ 21.752,13.566,20.835,10.612,15.355,21.55,26.911,33.207,41.233,49.664,58.386,67.528,85.89,118.46,193.18,290.95,391.85 /)
	miss = 0.0
	do i=1,16 !enkel eerste 16 punten vergelijken
		miss = miss + (abs(y_and_theta(4,i) - y_week(i+1))/y_week(i+1))
	enddo
end function

!maximale fout minimaliseren, dus deze methode berekend maximale miss
function calc_likelihood2(y_and_theta) result (maxmiss)
	implicit none
	real*8, dimension(4,K) :: y_and_theta
	real*8 :: maxmiss,miss
	integer :: i
	real*8, dimension(17) :: y_week
	y_week = (/ 21.752,13.566,20.835,10.612,15.355,21.55,26.911,33.207,41.233,49.664,58.386,67.528,85.89,118.46,193.18,290.95,391.85 /)
	maxmiss = 0.0
	do i=1,16 !enkel eerste 16 punten vergelijken
		miss = abs(y_and_theta(4,i) - y_week(i+1))/y_week(i+1)
		if (miss > maxmiss) then
			maxmiss = miss
		endif
	enddo
end function

function DBSSM_FUNC(K,dt, t0, beta,gamm,kappa,lambda,theta0) result(y_and_theta)
 use RANDOM
 use DIFF_SOLVER
 use DIFF

implicit none
real ( kind=selected_real_kind(8) )  :: dt,t0
real ( kind=selected_real_kind(8) )  :: theta0(3)
real ( kind=selected_real_kind(8) )   :: y_and_theta(4,K)
real ( kind=selected_real_kind(8) )   :: beta,gamm,kappa, lambda
integer 					   :: seed,clock,count_rate,count_max,i,K

 call init_random_seed()
 call rk4vec(3, t0, theta0, dt, fvec, y_and_theta(1:3,1),beta,gamm)
 y_and_theta(4,1)=rand_beta(lambda*y_and_theta(2,1), lambda*(1-y_and_theta(2,1)))
do i=2,K
 call init_random_seed()
 call rk4vec(3, t0, y_and_theta(1:3,i-1), dt, fvec, y_and_theta(1:3,i),beta,gamm)
 call drawDirichlet(3, kappa*y_and_theta(1:3,i),y_and_theta(1:3,i))
 y_and_theta(4,i)=rand_beta(lambda*y_and_theta(2,i), lambda*(1-y_and_theta(2,i)))
!print *,y_and_theta(4,i)
end do
y_and_theta(4,:) = y_and_theta(4,:)*100000
end function


end program


