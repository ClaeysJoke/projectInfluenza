module Sampling
  use Normal
  use Pdf
  use GradientModule
  use ODE
  use GetModule
  implicit none
  real*8,parameter :: PI = 4*ATAN(1.0_8)

contains
  subroutine GenerateSample(X,gradX,pdfX,y,dt,nbTries)
    real*8, intent(inout) :: X(:), gradX(:), pdfX
    real*8, intent(in) :: y(:)
    real*8, intent(in) :: dt
    integer,intent(out) :: nbTries
    real*8, allocatable :: trialX(:),gradTrialX(:)


    real*8 :: fromTrial, toTrial
    real*8 :: U
    real*8 :: pdfTrialX

    allocate(trialX(size(X)),gradTrialX(size(gradX)))
    nbTries = 0
    print *, "================= Generating Sample ====================="
    do !Loop until acceptable sample is generated
      print *, "================ Trial ", nbTries, " ============"
    ! Generate trial
      call GenerateTrial(X,gradX,dt,trialX)
    ! Calculate probability of the new trial, if p=0, skip entire loop.
      call LogPosterior(trialX, y, pdfTrialX)

      print *, "LogPosterior density", pdfTrialX

      if (pdfTrialX > 0) then

    ! Calculate trial gradient
        call Gradient(trialX,y,gradTrialX)

    ! Sample Uniform distribution for accept-reject
        call RANDOM_NUMBER(U)

    ! Call functions for accept reject

        fromTrial =  Transition(trialX,gradTrialX,dt,X)
        toTrial = Transition(X,gradX,dt,trialX)

    ! Accept reject
        if (U <= pdfTrialX*fromTrial/(pdfX*toTrial)) then
          X = trialX
          gradX = gradTrialX
          pdfX = pdfTrialX
          exit
        endif
      endif
      nbTries = nbTries + 1
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


  real*8 function Transition(X, gradX, dt, destination)
    real*8 :: X(:), gradX(:), dt, destination(:)
    real*8, allocatable :: mu(:)

    allocate(mu(size(X)))
    mu = destination - (X-gradX*dt)

    Transition = EXP(-0.5*(DOT_PRODUCT(destination-mu,destination-mu)/(2*dt))) &
    /(sqrt((2*PI)**size(X)*(2*dt)**size(X)))

    deallocate(mu)
  end function

  ! Calculates posterior pdf for a given parameter vector X
  ! Assumes fixed values of X to still be fixed
  subroutine LogPosterior(X,y,pdfX)
    real*8, intent(in) :: X(:),y(:)
    real*8, intent(out) :: pdfX
    real*8 :: prior, simulation, ILI
    integer :: i
    real*8 :: thetaI
    real*8,dimension(3) :: dirichletWeights

    ! Initialize to zero
    print *, "------------- Sampling Posterior --------------"
    prior = 0.0D0
    simulation = 0.0D0
    ILI = 0.0D0

    ! Prior pdf
    prior = prior + log(pdfGamma(X(1)))
    print *, "log p(kappa) ", prior
    prior = prior + log(pdfGamma(X(2)))
    print *, "log p(kappa)*p(lambda) ", prior
    prior = prior + log(pdfBeta(X(4),1.62D0,7084.1D0))
    print *, "log p(kappa)*p(lambda)*p(theta) ", prior
    prior = prior + log(pdfTruncatedNormal(X(6:7),X(4)))
    print *, "log p(kappa)*p(lambda)*p(theta)*p(z|theta) ", prior

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

end module
