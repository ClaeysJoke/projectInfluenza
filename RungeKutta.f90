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


subroutine rk4vec (m, t0, u0, dt, fvec, u )

  implicit none

  integer :: m

  real ( kind=selected_real_kind(8) ) :: dt,f0(m),f1(m),f2(m),f3(m)
  real ( kind=selected_real_kind(8) ) :: t0,t1,t2,t3
  real ( kind=selected_real_kind(8) ) :: u(m),u0(m),u1(m),u2(m),u3(m)

  interface

  subroutine fvec(m,t,u,uprime)
  integer:: m
  real ( kind=selected_real_kind(8) )  :: t,u(m),uprime(m)
  end subroutine fvec

  end interface  

  call fvec (m, t0, u0, f0 )

  t1 = t0 + 0.5*dt
  u1(1:m) = u0(1:m) + 0.5*dt * f0(1:m)
  call fvec (m, t1, u1, f1 )
  u2(1:m) = u0(1:m) + 0.5*dt * f1(1:m)
  call fvec (m, t1, u2, f2 )

  t3 = t0 + dt
  u3(1:m) = u0(1:m) + dt * f2(1:m)
  call fvec (m, t3, u3, f3 )
  u(1:m) = u0(1:m) + ( dt /6.0 ) * ( &
                 f0(1:m) &
     + 2.0* f1(1:m) &
     + 2.0* f2(1:m) &
     +           f3(1:m) )

  return
end

!EQUATION TO BE SOLVED GOES HERE
subroutine f (t,u,uprime)
real (kind=selected_real_kind(8) ) u
real (kind=selected_real_kind(8) ) t
real (kind=selected_real_kind(8) ) uprime
uprime = 3.0*t**2
end subroutine f
!VECTORIAL VERSION OF DIFFERENTIAL EQUATION
subroutine fvec (m,t,u,uprime)
integer :: m
real (kind=selected_real_kind(8) ) :: u(m), uprime(m)
real (kind=selected_real_kind(8) ) :: t
uprime(1) = -2.0*u(1)*u(2)
uprime(2) = 2.0*u(1)*u(2)-1.4*u(2)
uprime(3) = 1.4*u(2)
end subroutine fvec

program test_prog
real ( kind=selected_real_kind(8) ) u(3), ui(3)
real ( kind=selected_real_kind(8) ):: dt,ti
external :: fvec
dt=1
ui(1)=0.9
ui(2)=0.0002
ui(3)=0.0998
ti=0
do i=0,40
call rk4vec (3,ti+i*dt, ui, dt, fvec, u)
print*, ti+(i+1)*dt,u, u(1)+u(2)+u(3)
ui=u
end do
end program


