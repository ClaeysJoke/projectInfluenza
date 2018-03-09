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

subroutine fwdEulerVec (m,t0, u0, dt, f, u)
implicit none
integer :: m
real ( kind = selected_real_kind(8) ) :: dt,t0,f0(m),u(m),u0(m)

interface
subroutine fvec(m, t,u,uprime)
integer :: m
real ( kind = selected_real_kind(8) ) :: t,u,uprime
end subroutine fvec
end interface

call f(m,t0,u0,f0)
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
end


subroutine rk4Vec (m, t0, u0, dt, fvec, u,beta,gamm)

  implicit none

  integer :: m

  real ( kind=selected_real_kind(8) ) :: dt,f0(m),f1(m),f2(m),f3(m)
  real ( kind=selected_real_kind(8) ) :: t0,t1,t2,t3
  real ( kind=selected_real_kind(8) ) :: u(m),u0(m),u1(m),u2(m),u3(m)
  real ( kind=selected_real_kind(8) ) :: beta,gamm

  interface

  subroutine fvec(m,t,u,uprime,beta,gamm)
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
end
end module
!******************************************************************
!		    PART TWO: DIFFERENTIAL EQUATIONS
!				
!******************************************************************


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

!******************************************************************
!		  PART THREE: RANDOM SAMPLERS
!				
!******************************************************************


subroutine drawDirichlet(m,param,theta,seed)
!Subroutine that draws from Dirichlet using gamma distribution Gamm(alpha,1)
integer, intent(in) :: m
real (kind=selected_real_kind(8) ), intent(in) :: param(m)
real (kind=selected_real_kind(8) ), intent(out) ::theta(m)
real (kind=selected_real_kind(8) ) :: y,s
integer :: seed
s=0.0
do i=1,m
call drawGamma(param(i),y,seed)
theta(i)=y
s=s+y
end do
theta=(1.0/s)*theta
end subroutine drawDirichlet


subroutine drawGamma(alpha,X, seed)
!Draws from Gamm(alpha,1)
!Simple Wilson-Hilferty algorithm (see http://www.nrbook.com/devroye/Devroye_files/chapter_four.pdf )
implicit none
real (kind=selected_real_kind(8) ) :: alpha,sigma,sigma2, z0,X,Z,U,N,normal,m,s,t
integer :: seed,i

do i=1,10000
t=t+1
end do

Z=-0.5
i=0
m=0.0
s=1.0
sigma2=(9*alpha ) *(1.0- ( 1 / (3*alpha) ))**(.3333)
sigma=sigma2**(.5)
z0=( (3*alpha-1) / (3*alpha) )**(.3333)

i=0

do while((Z<0 .or. U*EXP(-0.5*N**2)>((Z/z0)**(3*alpha - 1))*EXP(-1*alpha*(Z**3-z0**3))) .and. i<10000)
i=i+1
print *,i
N=normal(m, s, seed)
U=rand()
Z=z0+sigma*N
end do
X=alpha*Z**3
end subroutine

function normal( mu, sigma2, seed )
! draws from normal distribution
! Using transformation
  implicit none

  real ( kind = 8 ) mu
  real ( kind = 8 ) sigma2
  real ( kind = 8 ) u1, u2
  real ( kind = 8 ) normal
  real ( kind = 8 ), parameter :: pi = 4*ATAN(1.0)
  real ( kind = 8 ) U
  integer ( kind = 4 ) seed
  real ( kind = 8 ) Z
  call srand(seed)
  u1 = rand()
  u2 = rand()
  Z = sqrt ( - 2.0D+00 * log ( u1 ) ) * cos ( 2.0D+00 * pi * u2 )

  normal = mu + sigma2 * Z

  return
end




!******************************************************************
!		    		MAIN
!				
!******************************************************************

program test_prog
real ( kind=selected_real_kind(8) ) u(3), ui(3)
real ( kind=selected_real_kind(8) ):: dt,ti
real ( kind=selected_real_kind(8) ):: beta,gamm,X
integer :: seed,clock,count_rate,count_max,i

external :: fvec
dt=1
ui(1)=0.9
ui(2)=0.0002
ui(3)=0.0998
ti=0
beta=3.5
gamm=0.5
!do i=0,40
!call rk4vec (3,ti+i*dt, ui, dt, fvec, u, beta, gamm)
!print*, ti+(i+1)*dt,u, u(1)+u(2)+u(3)
!ui=u
!end do


end program


