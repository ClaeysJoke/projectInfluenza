module Prior
  use Random
  use Sampling
  use Normal
  use ODE
  use GetModule

contains

  subroutine SamplePrior(X,y)
    real*8,intent(inout) :: X(:)
    real*8,intent(in) :: y(:)
    integer :: i
    real*8 :: shape,rate
    real*8 :: alpha,beta

    shape = 2.0D0
    rate = 0.0001D0
    alpha = 1.62D0
    beta = 7084.1D0
    open(unit=10,file='/Users/Nick/Documents/2eMaster/Project/Code/projectInfluenza/Prior.txt')


    ! Sample Kappa
    print *, "============ PRIOR ================="

    !X(1) = rand_gamma(shape,1.0D0/rate)
    ! Initialize around found good value
    X(1) = 500.0D0

    ! Sample Lambda


    !X(2) = rand_gamma(shape,1.0D0/rate)
    ! Initialize around found value
    X(2) = 90.0D0
    !X(2) = 90.0D0 + 200.0D0
    ! "Sample" S0

    X(3) = 0.9D0
    ! Sample I0

    ! randGamma is still giving shitty values :'(
    X(4) = 0.00023D0

    !X(4) = 0.0003D0


    ! "Sample" R0

    X(5) = 1.0D0 - X(3) - X(4)

    ! "Sample" PI and PT

    X(6) = 0.0144D0
    X(7) = 17.9D0


    ! "Sample" rho


    X(8) = getRho(X(3),X(4),X(6))
    !X(8) = 0.7D0
    ! "Sample" beta

    X(9) = getBeta(X)
    !X(9) = 2.1D0
    ! "Sample" thetas, by simulating them with the ODE solver from theta0


    do i = 1,size(y)
      select case(i)
      case(1)
        X(10:12) = rkvec(X(3:5),X(9),X(8)*X(9))
      case default
        X(10+3*(i-1):12+3*(i-1)) = rkvec(X(10+3*(i-2):12+3*(i-2)),X(9),X(8)*X(9))
      end select

    enddo
    write(10,*) X(8),X(9),X(10),X(11),X(12),X(13),X(14),X(15)
    !print *, "rho is",X(8)
  end subroutine





end module
