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


subroutine rk4vec ( t0, m, u0, dt, f, u )

  implicit none

  integer ( kind = selected_int_kind(4) ) :: m

  real ( kind=selected_real_kind(8) ) :: dt,f0(m),f1(m),f2(m),f3(m)
  real ( kind=selected_real_kind(8) ) :: t0,t1,t2,t3
  real ( kind=selected_real_kind(8) ) :: u(m),u0(m),u1(m),u2(m),u3(m)

  interface
  subroutine fvec(t,m,u,uprime)
  real ( kind=selected_real_kind(8) )  :: t,u,uprime
  integer(kind = selected_int_kind(4)) :: m
  end subroutine fvec
  end interface  

  call f ( t0, m, u0, f0 )

  t1 = t0 + 0.5*dt
  u1(1:m) = u0(1:m) + 0.5*dt * f0(1:m)
  call f ( t1, m, u1, f1 )

  t2 = t0 + 0.5*dt
  u2(1:m) = u0(1:m) + 0.5*dt * f1(1:m)
  call f ( t2, m, u2, f2 )

  t3 = t0 + dt
  u3(1:m) = u0(1:m) + dt * f2(1:m)
  call f ( t3, m, u3, f3 )
!
!  Combine them to estimate the solution U at time T1.
!
  u(1:m) = u0(1:m) + ( dt /6.0 ) * ( &
                 f0(1:m) &
     + 2.0* f1(1:m) &
     + 2.0* f2(1:m) &
     +           f3(1:m) )

  return
end

!EQUATION TO BE SOLVED GOES HERE
subroutine f (t,u,uprime)
real (kind=8) u
real (kind=8) t
real (kind=8) uprime
uprime = 2.0*t
end subroutine f

subroutine fvec (t,u,uprime)
real (kind=8) u
real (kind=8) t
real (kind=8) uprime
integer(kind= selected_int_kind(4)) :: m
uprime = 2.0*t
end subroutine fvec

program test_prog
real ( kind = 8 ) u, ui
real ( kind =8) dt,t0
external :: f
dt=0.1
ui=1
t0=1
do i=0,10
call rk4 (t0+i*dt, ui, dt, f, u)
print*, t0+(i+1)*dt,u
ui=u
end do
end program



