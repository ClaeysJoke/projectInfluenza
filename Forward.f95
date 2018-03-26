program ForwardSolver
  use ODE
  use Random

  implicit none

  integer :: nb_samples = 500
  integer :: i,j
  real*8 :: thetaI
  real*8,dimension(3) :: theta
  real*8,allocatable :: hiddenTheta(:)
  real*8,dimension(50) :: y
  integer :: cutoff = 10
  real*8 :: kappa,lambda,S0,I0,R0,PI,PT,rho,beta
  real*8 :: one = 1.0D0

  open(unit=10,file='/Users/Nick/Documents/2eMaster/Project/Code/projectInfluenza/Samples/Samples_10.txt')
  open(unit=11,file='/Users/Nick/Documents/2eMaster/Project/Code/projectInfluenza/Samples/ILI_10.txt')

  allocate(hiddenTheta(cutoff))

  do i = 1,nb_samples
    read(10,*) kappa,lambda,S0,I0,R0,PI,PT,rho,beta
    read(10,*) hiddenTheta(1),hiddenTheta(2),hiddenTheta(3),hiddenTheta(4),hiddenTheta(5),hiddenTheta(6)&
    ,hiddenTheta(7),hiddenTheta(8),hiddenTheta(9),theta(1),hiddenTheta(10),theta(3)!hiddenTheta(11),theta(1),hiddenTheta(12),theta(3)

    do j = 1,50
      if (j .le. cutoff) then
        y(j) = rand_beta(lambda*hiddenTheta(j),lambda*(one-hiddenTheta(j)))
      elseif ( j == cutoff+1 ) then
        theta(2) = hiddenTheta(cutoff)

        call rand_Dirichlet(3,kappa*rkvec(theta,beta,rho*beta),theta)
        y(j) = rand_beta(lambda*theta(2),lambda*(one-theta(2)))
      else
        call rand_Dirichlet(3,kappa*rkvec(theta,beta,rho*beta),theta)
        y(j) = rand_beta(lambda*theta(2),lambda*(one-theta(2)))

      endif

    enddo
    write(11,*) y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),y(12),y(13),y(14),y(15),y(16)&
    ,y(17),y(18),y(19),y(20),y(21),y(22),y(23),y(24),y(25),y(26)&
    ,y(27),y(28),y(29),y(30),y(31),y(32),y(33),y(34),y(35),y(36)&
    ,y(37),y(38),y(39),y(40),y(41),y(42),y(43),y(44),y(45),y(46)&
    ,y(47),y(48),y(49),y(50)

  enddo

  deallocate(hiddenTheta)
  close(10)
  close(11)

end program
