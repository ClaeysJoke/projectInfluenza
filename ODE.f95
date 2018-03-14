module ODE

  use

contains
  function rkvec(theta,beta,gamma)
    real*8, dimension(3) :: theta
    real*8 :: beta,gamma

    real*8,dimension(3) :: rkvec

    rkvec = theta



  end function



end module
