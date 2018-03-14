module ODE

  use DIFF_SOLVER


contains
  function rkvec(theta,beta,gamma)
    real*8, dimension(3) :: theta
    real*8 :: beta,gamma
    integer :: m
    real*8 :: t0,dt

    real*8,dimension(3) :: rkvec

    m = 3
    t0 = 0.0D0
    dt = 1.0D0

    call rk4vec(m,t0,theta,dt,fvec,rkvec,beta,gamma)

  end function



end module
