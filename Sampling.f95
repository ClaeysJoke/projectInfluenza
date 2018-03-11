module Sampling
  use Normal
  use Pdf
  implicit none
  real*8,parameter :: PI = 4*ATAN(1.0_8)

contains
  subroutine GenerateSample(X,gradX,pdfX,y,dt)
    real*8, intent(inout) :: X(:), gradX(:), pdfX
    real*8, intent(in) :: y(:)
    real*8, intent(in) :: dt
    real*8, allocatable :: trialX(:),gradTrialX(:)

    real*8 :: fromTrial, toTrial
    real*8 :: U
    real*8 :: pdfTrialX

    allocate(trialX(size(X)),gradTrialX(size(gradX)))

    do !Loop until acceptable sample is generated
    ! Generate trial
      call GenerateTrial(X,gradX,dt,trialX)
    ! Calculate trial gradient
      call Gradient(trialX,y,gradTrialX)

    ! Sample Uniform distribution for accept-reject
      call RANDOM_NUMBER(U)

    ! Call functions for accept reject
      call Posterior(trialX, y, pdfTrialX)
      fromTrial =  Transition(trialX,gradTrialX,dt,X)
      toTrial = Transition(X,gradX,dt,trialX)

    ! Accept reject
      if (U <= pdfTrialX*fromTrial/(pdfX*toTrial)) then
        X = trialX
        gradX = gradTrialX
        pdfX = pdfTrialX
        exit
      endif
    enddo

    deallocate(trialX,gradTrialX)
  end subroutine

  subroutine GenerateTrial(X,gradX,dt,trialX)
    real*8, intent(out) :: trialX(:)
    real*8, intent(in) :: X(:), gradX(:), dt
    integer :: i
    real*8, allocatable :: xi(:)

    allocate(xi(size(X)))

    do i = 1,size(xi)
      xi(i) = rand_normal(0.0D0,1.0D0)
    enddo

    trialX = X - gradX*dt + SQRT(2.0*dt)*xi

    ! Setting theta0 to passable values
    trialX(3) = 0.9D0
    trialX(5) = 1 - trialX(4) - trialX(3)

    ! Setting rho
    trialX(8) = getRho(X)
    ! Setting beta
    trialX(9) = getBeta(X)

    deallocate(xi)
  end subroutine

  subroutine Gradient(X,y,gradX)
    real*8, intent(in) :: X(:),y(:)
    real*8, intent(out) :: gradX(:)
    gradX = -X
  end subroutine


  real*8 function Transition(X, gradX, dt, destination)
    real*8 :: X(:), gradX(:), dt, destination(:)
    real*8, allocatable :: mu(:)

    allocate(mu(size(X)))
    mu = destination - (X-gradX*dt)

    Transition = EXP(-0.5*(sqrt(2*dt)**(-1)*DOT_PRODUCT(destination-mu,destination-mu))) &
    /(sqrt((2*PI)**size(X)*sqrt(2*dt)**size(X)))

    deallocate(mu)
  end function

  ! Calculates posterior pdf for a given parameter vector X
  ! Assumes fixed values of X to still be fixed
  subroutine Posterior(X,y,pdfX)
    real*8, intent(in) :: X(:),y(:)
    real*8, intent(out) :: pdfX
    real*8 :: prior, simulation, ILI
    integer :: i
    real*8 :: thetaI
    real*8,dimension(3) :: dirichletWeights

    ! Initialize to zero
    prior = 1
    simulation = 1
    ILI = 1

    ! Prior pdf
    prior = prior*pdfGamma(X(1))
    prior = prior*pdfGamma(X(2))
    prior = prior*pdfBeta(X(4),1.62D0,7084.1D0)
    prior = prior*pdfTruncatedNormal(X(6:7),X(4))

    ! ILI and simulation pdf
    do i = 1,size(y)
      thetaI = X(11 + 3*(i-1))
      ILI = ILI*pdfBeta(y(i),X(2)*thetaI,X(2)*(1.0D0-thetaI))
      select case (i)
      case(1)
        dirichletWeights = rkvec(X(3:5),X(9),X(8)*X(9))
        simulation = simulation*pdfDirichlet(X(10:12),X(1)*dirichletWeights)
      case default
        dirichletWeights = rkvec(X(10+3*(i-1):9+3*i),X(9),X(8)*X(9))
        simulation = simulation*pdfDirichlet(X(10+3*i:12+3*i),X(1)*dirichletWeights)
      end select
    enddo

    pdfX = prior*ILI*simulation

  end subroutine

  real*8 function getRho(X)
    real*8 :: X(:)
    real*8 :: PI,I0,S0,rho
    real*8 :: func, der_func
    S0 = X(3)
    I0 = X(4)
    PI = X(6)
    rho = X(8)

    ! Newton-Raphson loop
    do
      func = I0 + S0 - PI - rho*(log(S0)+1.0D0-log(rho))
      if (abs(func)<100*EPSILON(func)) then
        getRho = rho
        return
      endif
      der_func = -log(S0) + log(rho)
      rho = rho - func/der_func
    end do
  end function

  real*8 function getBeta(X)
    real*8 :: X(:)
    real*8 :: sigma = 0.0421D0
    real*8,dimension(17) :: tau,log_reg
    real*8 :: PT,I0,rho
    tau(1) = -49.7540D0
    tau(2) = -0.9577D0
    tau(3) = -0.0065D0
    tau(4) = -9.4896D0
    tau(5) = -0.3761D0
    tau(6) = -590.0001D0
    tau(7) = -2537.6102D0
    tau(8) = -4756.1828D0
    tau(9) = -3265.2458D0
    tau(10) = -102.2665D0
    tau(11) = -4.0162D0
    tau(12) = -430.9596D0
    tau(13) = -16.7104D0
    tau(14) = -798.3443D0
    tau(15) = -30.6638D0
    tau(16) = -543.8857D0
    tau(17) = -20.7459D0

    PT = X(7)
    I0 = X(4)
    rho = X(8)

    log_reg(1) = 1
    log_reg(2) = log(PT)
    log_reg(3) = log(PT)**2
    log_reg(4) = log(I0)
    log_reg(5) = log(I0)**2
    log_reg(6) = log(rho)
    log_reg(7) = log(rho)**2
    log_reg(8) = log(rho)**3
    log_reg(9) = log(rho)**4
    log_reg(10) = log(I0)*log(rho)
    log_reg(11) = log(I0)**2*log(rho)
    log_reg(12) = log(I0)*log(rho)**2
    log_reg(13) = log(I0)**2*log(rho)**2
    log_reg(14) = log(I0)*log(rho)**3
    log_reg(15) = log(I0)**2*log(rho)**3
    log_reg(16) = log(I0)*log(rho)**4
    log_reg(17) = log(I0)**2*log(rho)**4

    getBeta = exp(DOT_PRODUCT(log_reg,tau)+0.5D0*sigma**2)

  end function

end module
