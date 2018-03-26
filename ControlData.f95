program ControlData
  use ODE
  use Random
  use Normal

  implicit none


  integer :: T = 30
  integer :: k
  real*8 :: lambda,kappa,rho,beta
  real*8 :: S,I,R

  real*8,dimension(3) :: theta
  real*8 :: y

  ! Setting parameters

  lambda = 90.0D0
  kappa = 1500.0D0

  S = 0.9D0
  I = 0.0003D0
  R = 1.0D0 - S - I

  rho = 0.7D0
  beta = 2.1D0


  open(unit=7,file='/Users/Nick/Documents/2eMaster/Project/Code/projectInfluenza/Control.txt')

  ! First update
  theta(1) = S
  theta(2) = I
  theta(3) = R

! call init_random_seed()
!  call rand_Dirichlet(3,kappa*rkvec(theta,beta,rho*beta),theta)
!
  do k = 1,T
!    call rand_Dirichlet(3,kappa*rkvec(theta,beta,rho*beta),theta)
!
!    y = rand_Beta(lambda*theta(2),lambda*(1.0D0-theta(2)))
!

    theta = rkvec(theta,beta,rho*beta)
    y = theta(2)
    write(7,*) y

  enddo









end program
