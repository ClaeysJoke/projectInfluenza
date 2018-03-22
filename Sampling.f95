module Sampling
  use Normal
  use Pdf
  use Random
  use ODE
  use GetModule
  implicit none

  real*8 :: dt = 1
  real*8 :: h = 0.001D0

contains
  subroutine GenerateSample(X,pdfX,y,nbTries)
    real*8, intent(inout) :: X(:), pdfX
    real*8, intent(in) :: y(:)
    integer,intent(out) :: nbTries
    real*8, allocatable :: trialX(:)


    real*8 :: fromTrial, toTrial
    real*8 :: U
    real*8 :: pdfTrialX

    allocate(trialX(size(X)))
    nbTries = 0
    print *, "================= Generating Sample ====================="
    do !Loop until acceptable sample is generated
      print *, "================ Trial ", nbTries, " ============"
    ! Generate trial
      call GenerateTrial(X,y,pdfX,trialX)
    ! Calculate probability of the new trial, if p=0, skip entire loop.
      call LogPosterior(trialX, y, pdfTrialX)

      print *, "LogPosterior density", pdfTrialX

      if (pdfTrialX > -huge(pdfTrialX)) then

    ! Calculate trial gradient
        !call GradientFinite(trialX,y,pdfX,gradTrialX)

    ! Sample Uniform distribution for accept-reject
        call RANDOM_NUMBER(U)

    ! Call functions for accept reject

        fromTrial =  Transition(trialX,X,y)
        toTrial = Transition(X,trialX,y)

    ! Accept reject
        if (log(U) <= pdfTrialX + log(fromTrial) - pdfX -log(toTrial)) then
          X = trialX

          pdfX = pdfTrialX
          exit
        endif
      endif
      nbTries = nbTries + 1
    enddo

    deallocate(trialX)
  end subroutine

  subroutine GenerateTrial(X,y,pdfX,trialX)
    real*8, intent(out) :: trialX(:)
    real*8, intent(in) :: X(:)
    real*8, intent(in) :: y(:)
    real*8 :: pdfX
    integer :: i
    real*8, allocatable :: gradientStep(:)
    real*8 :: pdfPlus,pdfMinus
    real*8 :: gradient
    real*8,allocatable :: X_placeholder(:)
    real*8,dimension(3) :: theta_estimate
    real*8 :: variance
!    real*8, allocatable :: xi(:)
    allocate(gradientStep(size(X)))
    allocate(X_placeholder(size(X)))
    ! Updating kappa
    print *, "Updating Kappa"
      ! Grad V of kappa
      gradientStep = 0.0D0
      gradientStep(1) = 0.5*h
      call LogPosterior(X+gradientStep,y,pdfPlus)
      call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
      trialX(1) = X(1) - dt*(-pdfPlus+pdfMinus)/h + sqrt(2*dt)*rand_normal(0.0D0,1.0D0)
    ! Updating lambda
      ! Grad V of lambda
      print *, "Updating Lambda"
      gradientStep = 0.0D0
      gradientStep(2) = 0.5*h
      call LogPosterior(X+gradientStep,y,pdfPlus)
      call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
      trialX(2) = X(2) - dt*(-pdfPlus+pdfMinus)/h + sqrt(2*dt)*rand_normal(0.0D0,1.0D0)

    ! Updating S0
    print *, "Updating S0"
      trialX(3) = 0.9D0

    ! Updating I0
    print *, "Updating I0"
    gradientStep = 0.0D0
    gradientStep(4) = 0.5*h
    call LogPosterior(X+gradientStep,y,pdfPlus)
    call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
    gradient = (-pdfPlus+pdfMinus)/abs(-pdfPlus+pdfMinus)
    variance = 0.00018**2
    trialX(4) = X(4) - sqrt(variance)*gradient + sqrt(variance)*rand_normal(0.0D0,1.0D0)
    if (trialX(4) < tiny(0.0D0)) then
      trialX(4) = tiny(0.0D0)
    elseif (trialX(4) > 0.1D0) then
      trialX(4) = 0.1D0
    endif

    ! Updating R0
    print *, "Updating R0"
      trialX(5) = 1.0D0 - trialX(3) - trialX(4)

    ! Updating PI
    print *, "Updating PI"
    gradientStep = 0.0D0
    gradientStep(6) = 0.5*h
    call LogPosterior(X+gradientStep,y,pdfPlus)
    call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
    gradient = (-pdfPlus+pdfMinus)/abs(-pdfPlus+pdfMinus)
    variance = 0.000036D0
    trialX(6) = X(6) - sqrt(variance)*gradient + sqrt(variance)*rand_normal(0.0D0,1.0D0)
    if (trialX(6)<trialX(4)) then
      trialX(6) = trialX(4)
    elseif (trialX(6) > 1.0D0) then
      trialX(6) = 1.0D0
    endif

    ! Updating PT
    print *, "Updating PT"
    gradientStep = 0.0D0
    gradientStep(7) = 0.5*h
    call LogPosterior(X+gradientStep,y,pdfPlus)
    call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
    gradient = (-pdfPlus+pdfMinus)/abs(-pdfPlus+pdfMinus)
    variance = 16.09
    trialX(7) = X(7) - sqrt(variance)*gradient + sqrt(variance)*rand_normal(0.0D0,1.0D0)
    if (trialX(7)<1.0D0) then
      trialX(7) = 1.0D0
    elseif (trialX(7) > 35.0D0) then
      trialX(7) = 35.0D0
    endif

    ! Updating rho
    print *, "Updating Rho"
      trialX(8) = getRhoIter(trialX)
    ! Updating beta
    print *, "Updating Beta"
      trialX(9) = getBeta(trialX)

    ! Updating thetas
      do i = 1,size(y)
      select case(i)
      case(1)
        print *,"Updating theta 1"
        call rand_Dirichlet(3,trialX(1)*rkvec(trialX(3:5),trialX(9),trialX(8)*trialX(9)),theta_estimate)
        trialX(10:12) = theta_estimate
      case default
        print *, "Updating theta",i
        call rand_Dirichlet(3,trialX(1)*rkvec(trialX(10+3*(i-2):12+3*(i-2)),trialX(9),trialX(8)*trialX(9)),theta_estimate)
        trialX(10+3*(i-1):12+3*(i-1)) = theta_estimate
      end select
    enddo

    !allocate(xi(size(X)))

    !do i = 1,size(xi)
!      xi(i) = rand_normal(0.0D0,1.0D0)
!    enddo
!
!    trialX = X - gradX*dt + SQRT(2.0*dt)*xi
!
!    ! Setting theta0 to passable values
!    trialX(3) = 0.9D0
!    if (trialX(4) > 0.1D0) then
!      trialX(4) = 0.1D0
!    endif
!    trialX(5) = 1 - trialX(4) - trialX(3)
!
!    ! Setting rho
!    trialX(8) = getRho(X)
!    ! Setting beta
!    trialX(9) = getBeta(X)

!    deallocate(xi)
    deallocate(gradientStep)
    deallocate(X_placeholder)
  end subroutine

  ! Transition distribution for phi vector (for starters)
  real*8 function Transition(X,destination,y)
    real*8 :: X(:),destination(:),y(:)
    real*8 :: kappa
    real*8 :: lambda
    real*8 :: gradient
    real*8 :: mu,sigma2
    real*8 :: kappaTransition,lambdaTransition,I0Transition,PITransition,PTTransition
    real*8 :: pdfPlus,pdfMinus
    real*8,allocatable :: gradientStep(:)

    allocate(gradientStep(size(X)))

    ! Kappa
    gradientStep = 0.0D0
    gradientStep(1) = 0.5*h
    call LogPosterior(X+gradientStep,y,pdfPlus)
    call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
    mu =  X(1) - dt*(-pdfPlus+pdfMinus)/h
    sigma2 = 2*dt

    kappaTransition = exp(-(mu-destination(1))**2/(2*sigma2))/sqrt(2*PI*sigma2)

    ! Lambda
    gradientStep = 0.0D0
    gradientStep(2) = 0.5*h
    call LogPosterior(X+gradientStep,y,pdfPlus)
    call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
    mu =  X(2) - dt*(-pdfPlus+pdfMinus)/h
    sigma2 = 2*dt

    kappaTransition = exp(-(mu-destination(2))**2/(2*sigma2))/sqrt(2*PI*sigma2)

    ! I0
    gradientStep = 0.0D0
    gradientStep(4) = 0.5*h
    call LogPosterior(X+gradientStep,y,pdfPlus)
    call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
    gradient = (-pdfPlus+pdfMinus)/abs(-pdfPlus+pdfMinus)
    sigma2 = 0.00018**2
    mu = X(4) - sqrt(sigma2)*gradient

    I0Transition = exp(-(mu-destination(4))**2/(2*sigma2))/sqrt(2*PI*sigma2)

    ! PI
    gradientStep = 0.0D0
    gradientStep(6) = 0.5*h
    call LogPosterior(X+gradientStep,y,pdfPlus)
    call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
    gradient = (-pdfPlus+pdfMinus)/abs(-pdfPlus+pdfMinus)
    sigma2 = 0.000036D0
    mu = X(6) - sqrt(sigma2)*gradient

    PITransition = exp(-(mu-destination(6))**2/(2*sigma2))/sqrt(2*PI*sigma2)

    ! PT
    gradientStep = 0.0D0
    gradientStep(7) = 0.5*h
    call LogPosterior(X+gradientStep,y,pdfPlus)
    call LogPosterior(X-2.0D0*gradientStep,y,pdfMinus)
    gradient = (-pdfPlus+pdfMinus)/abs(-pdfPlus+pdfMinus)
    sigma2 = 16.09
    mu = X(7) - sqrt(sigma2)*gradient

    PTTransition = exp(-(mu-destination(7))**2/(2*sigma2))/sqrt(2*PI*sigma2)

    deallocate(gradientStep)
  end function

end module
