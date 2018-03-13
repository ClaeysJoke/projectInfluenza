module ODE

contains
  function rkvec(theta,beta,gamma)
    real*8, dimension(3) :: theta
    real*8 :: beta,gamma

    real*8,dimension(3) :: rkvec

    ! DUMMY FUNCTION
    rkvec = theta

  end function



end module
